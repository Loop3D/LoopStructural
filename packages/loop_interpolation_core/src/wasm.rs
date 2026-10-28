//! WebAssembly bindings (wasm-bindgen). Compiled only with the feature `wasm`.
//!
//! # JavaScript signatures (used by `web/worker.js`)
//!
//! ```text
//! memoryBytes() -> number                 // size of the wasm memory in bytes
//!
//! new Problem(origin: Float64Array(3), counts: Uint32Array(3), step: Float64Array(3))
//! problem.setRegularisation(w: object)    // keys, all optional: dxx dyy dzz dxy dxz dyz dx dy dz
//!                                         // (missing key = 1.0; dx/dy/dz are the border weights)
//! problem.setRidgeFactor(ridgeFactor: number)          // 0 removes the ridge term; default 1e-8
//! problem.setInitialGuess(x0: Float64Array(nodes))     // start value; default 0
//! problem.addValue(points: Float64Array(3n), values: Float64Array(n), weights: Float64Array(n))
//! problem.addGradient(points, vectors: Float64Array(3n), weights)
//! problem.addNormal(points, vectors, weights)
//! problem.addTangent(points, vectors, weights)
//! problem.addInterface(points, ids: Float64Array(n), weight: number)
//! problem.addValueBound(points, lower, upper, weights)
//! problem.addPairBound(pointsA, pointsB, lower, upper, weights)
//! problem.estimate(settings: object|undefined) -> Estimate { total, vectors, constraints, working }
//! problem.start(settings: object|undefined) -> Session   // consumes the problem
//!
//! settings keys, all optional:
//!   precision: 'f64' | 'f32', tol, maxIter, budgetBytes, rho,
//!   innerMaxIter, innerTol, admmAbsTol, admmRelTol
//!
//! session.step(maxIters: number) -> StepResult
//!   { done: boolean, stage: 'setup'|'solve'|'check'|'done', iteration: number,
//!     residual: number, violation: number,
//!     converged: boolean, peakBytes: number, warnings: string[] }
//!   (converged, peakBytes and warnings are valid when done is true.
//!    After done, step keeps returning done.)
//! session.field() -> Float64Array | Float32Array   // one copy out
//! session.evaluate(points: Float64Array(3m)) -> Float64Array(m)
//! session.evaluateGradient(points: Float64Array(3m)) -> Float64Array(3m)
//! session.free()                                    // drops the session
//! ```
//!
//! # Errors
//!
//! Every method that can fail throws a `JsError`. Its `message` is JSON text
//! with a `kind` key, for example
//! `{"kind":"OverBudget","estimate":123,"budget":100,"maxNodesPerAxis":40}`.
//! The kinds are `InvalidInput` (`constraint`, `reason`), `GridTooSmall`
//! (`axis`, `nodes`), `GridTooLarge` (`nodes`, `maxNodes`), `NoConstraints`,
//! `OverBudget`, `OutOfMemory` (`requestedBytes`), `Infeasible` (`worst`: list
//! of `{index, kind, violation}`), `PointOutsideGrid` (`index`),
//! `Unsupported` (`setting`) and `SolveRunning`. A non-finite number is `null`.

use wasm_bindgen::prelude::*;

use crate::error::{BoundKind, CoreError};
use crate::field::{Field, FieldData};
use crate::interpolator::{BorderRegularisationWeights, RegularisationWeights};
use crate::problem::Problem;
use crate::report::{SolveReport, Stage, Step};
use crate::session::SolveSession;
use crate::settings::{Precision, Settings};

// ---------------------------------------------------------------------------
// JSON text for errors (no serde)
// ---------------------------------------------------------------------------

/// Write `s` as a JSON string with quotes and escapes.
fn json_str(out: &mut String, s: &str) {
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
}

/// Write a number as JSON. A non-finite number is `null`.
fn json_num(out: &mut String, x: f64) {
    if x.is_finite() {
        out.push_str(&format!("{x}"));
    } else {
        out.push_str("null");
    }
}

/// Write `,"key":` (the key is a fixed ASCII name).
fn json_key(out: &mut String, key: &str) {
    out.push_str(",\"");
    out.push_str(key);
    out.push_str("\":");
}

/// Make the JSON text of an error.
fn error_json(e: &CoreError) -> String {
    let mut o = String::from("{\"kind\":");
    match e {
        CoreError::InvalidInput { constraint, reason } => {
            json_str(&mut o, "InvalidInput");
            json_key(&mut o, "constraint");
            json_str(&mut o, constraint);
            json_key(&mut o, "reason");
            json_str(&mut o, reason);
        }
        CoreError::GridTooSmall { axis, nodes } => {
            json_str(&mut o, "GridTooSmall");
            json_key(&mut o, "axis");
            json_num(&mut o, *axis as f64);
            json_key(&mut o, "nodes");
            json_num(&mut o, *nodes as f64);
        }
        CoreError::GridTooLarge { nodes, max_nodes } => {
            json_str(&mut o, "GridTooLarge");
            json_key(&mut o, "nodes");
            json_num(&mut o, *nodes as f64);
            json_key(&mut o, "maxNodes");
            json_num(&mut o, *max_nodes as f64);
        }
        CoreError::NoConstraints => json_str(&mut o, "NoConstraints"),
        CoreError::OverBudget {
            estimate_bytes,
            budget_bytes,
            max_nodes_per_axis,
        } => {
            json_str(&mut o, "OverBudget");
            json_key(&mut o, "estimate");
            json_num(&mut o, *estimate_bytes as f64);
            json_key(&mut o, "budget");
            json_num(&mut o, *budget_bytes as f64);
            json_key(&mut o, "maxNodesPerAxis");
            json_num(&mut o, *max_nodes_per_axis as f64);
        }
        CoreError::OutOfMemory { requested_bytes } => {
            json_str(&mut o, "OutOfMemory");
            json_key(&mut o, "requestedBytes");
            json_num(&mut o, *requested_bytes as f64);
        }
        CoreError::Infeasible { worst } => {
            json_str(&mut o, "Infeasible");
            json_key(&mut o, "worst");
            o.push('[');
            for (i, w) in worst.iter().enumerate() {
                if i > 0 {
                    o.push(',');
                }
                o.push_str("{\"index\":");
                json_num(&mut o, w.index as f64);
                json_key(&mut o, "kind");
                json_str(
                    &mut o,
                    match w.kind {
                        BoundKind::Value => "value",
                        BoundKind::Pair => "pair",
                    },
                );
                json_key(&mut o, "violation");
                json_num(&mut o, w.violation);
                o.push('}');
            }
            o.push(']');
        }
        CoreError::PointOutsideGrid { index } => {
            json_str(&mut o, "PointOutsideGrid");
            json_key(&mut o, "index");
            json_num(&mut o, *index as f64);
        }
        CoreError::Unsupported { setting } => {
            json_str(&mut o, "Unsupported");
            json_key(&mut o, "setting");
            json_str(&mut o, setting);
        }
        CoreError::SolveRunning => json_str(&mut o, "SolveRunning"),
    }
    o.push('}');
    o
}

/// Make a thrown JavaScript error from a core error.
fn js_err(e: CoreError) -> JsError {
    JsError::new(&error_json(&e))
}

/// Make a thrown error for a wrong argument of the binding.
fn arg_err(constraint: &str, reason: &str) -> JsError {
    js_err(CoreError::InvalidInput {
        constraint: constraint.into(),
        reason: reason.into(),
    })
}

// ---------------------------------------------------------------------------
// Plain JavaScript objects (settings and weights), read with getters
// ---------------------------------------------------------------------------

#[wasm_bindgen]
extern "C" {
    /// A plain JavaScript object with the settings keys.
    pub type SettingsObj;
    #[wasm_bindgen(method, getter)]
    fn precision(this: &SettingsObj) -> Option<String>;
    #[wasm_bindgen(method, getter)]
    fn tol(this: &SettingsObj) -> Option<f64>;
    #[wasm_bindgen(method, getter, js_name = maxIter)]
    fn max_iter(this: &SettingsObj) -> Option<f64>;
    #[wasm_bindgen(method, getter, js_name = budgetBytes)]
    fn budget_bytes(this: &SettingsObj) -> Option<f64>;
    #[wasm_bindgen(method, getter)]
    fn rho(this: &SettingsObj) -> Option<f64>;
    #[wasm_bindgen(method, getter, js_name = innerMaxIter)]
    fn inner_max_iter(this: &SettingsObj) -> Option<f64>;
    #[wasm_bindgen(method, getter, js_name = innerTol)]
    fn inner_tol(this: &SettingsObj) -> Option<f64>;
    #[wasm_bindgen(method, getter, js_name = admmAbsTol)]
    fn admm_abs_tol(this: &SettingsObj) -> Option<f64>;
    #[wasm_bindgen(method, getter, js_name = admmRelTol)]
    fn admm_rel_tol(this: &SettingsObj) -> Option<f64>;

    /// A plain JavaScript object with the regularisation weight keys.
    pub type WeightsObj;
    #[wasm_bindgen(method, getter)]
    fn dxx(this: &WeightsObj) -> Option<f64>;
    #[wasm_bindgen(method, getter)]
    fn dyy(this: &WeightsObj) -> Option<f64>;
    #[wasm_bindgen(method, getter)]
    fn dzz(this: &WeightsObj) -> Option<f64>;
    #[wasm_bindgen(method, getter)]
    fn dxy(this: &WeightsObj) -> Option<f64>;
    #[wasm_bindgen(method, getter)]
    fn dxz(this: &WeightsObj) -> Option<f64>;
    #[wasm_bindgen(method, getter)]
    fn dyz(this: &WeightsObj) -> Option<f64>;
    #[wasm_bindgen(method, getter)]
    fn dx(this: &WeightsObj) -> Option<f64>;
    #[wasm_bindgen(method, getter)]
    fn dy(this: &WeightsObj) -> Option<f64>;
    #[wasm_bindgen(method, getter)]
    fn dz(this: &WeightsObj) -> Option<f64>;
}

/// Convert a JavaScript number to a count. Reject values that are not whole.
fn to_count(name: &str, x: f64) -> Result<usize, JsError> {
    if x.is_finite() && x >= 0.0 && x.fract() == 0.0 && x <= 9.0e15 {
        Ok(x as usize)
    } else {
        Err(arg_err(name, "must be a whole number that is 0 or more"))
    }
}

/// Read the settings from a JavaScript object. `undefined` gives the defaults.
fn parse_settings(obj: Option<SettingsObj>) -> Result<Settings, JsError> {
    let mut s = Settings::default();
    let Some(o) = obj else {
        return Ok(s);
    };
    if let Some(p) = o.precision() {
        s.precision = match p.as_str() {
            "f64" => Precision::F64,
            "f32" => Precision::F32,
            _ => return Err(arg_err("settings.precision", "must be 'f64' or 'f32'")),
        };
    }
    if let Some(v) = o.tol() {
        s.tol = v;
    }
    if let Some(v) = o.max_iter() {
        s.max_iter = to_count("settings.maxIter", v)?;
    }
    if let Some(v) = o.budget_bytes() {
        s.budget_bytes = Some(to_count("settings.budgetBytes", v)? as u64);
    }
    if let Some(v) = o.rho() {
        s.rho = v;
    }
    if let Some(v) = o.inner_max_iter() {
        s.inner_max_iter = to_count("settings.innerMaxIter", v)?;
    }
    if let Some(v) = o.inner_tol() {
        s.inner_tol = v;
    }
    if let Some(v) = o.admm_abs_tol() {
        s.admm_abs_tol = v;
    }
    if let Some(v) = o.admm_rel_tol() {
        s.admm_rel_tol = v;
    }
    Ok(s)
}

/// Convert a typed array of 3 values to an array.
fn triple(name: &str, v: &[f64]) -> Result<[f64; 3], JsError> {
    <[f64; 3]>::try_from(v).map_err(|_| arg_err(name, "must have 3 values"))
}

// ---------------------------------------------------------------------------
// Exported types
// ---------------------------------------------------------------------------

/// Memory estimate, in bytes. All values are JavaScript numbers.
#[wasm_bindgen(js_name = Estimate)]
pub struct WasmEstimate {
    pub total: f64,
    pub vectors: f64,
    pub constraints: f64,
    pub working: f64,
}

/// Result of one `step` call.
/// Only the `String` and `Vec` fields need a clone in their getter.
#[wasm_bindgen(js_name = StepResult)]
pub struct WasmStepResult {
    pub done: bool,
    #[wasm_bindgen(getter_with_clone)]
    pub stage: String,
    pub iteration: f64,
    pub residual: f64,
    pub violation: f64,
    pub converged: bool,
    #[wasm_bindgen(js_name = peakBytes)]
    pub peak_bytes: f64,
    #[wasm_bindgen(getter_with_clone)]
    pub warnings: Vec<String>,
}

/// An interpolation problem (JavaScript name `Problem`).
#[wasm_bindgen(js_name = Problem)]
pub struct WasmProblem {
    inner: Problem,
}

#[wasm_bindgen(js_class = Problem)]
impl WasmProblem {
    #[wasm_bindgen(constructor)]
    pub fn new(origin: &[f64], counts: &[u32], step: &[f64]) -> Result<WasmProblem, JsError> {
        let origin = triple("origin", origin)?;
        let step = triple("step", step)?;
        let counts = <[u32; 3]>::try_from(counts)
            .map_err(|_| arg_err("counts", "must have 3 values"))?
            .map(|c| c as usize);
        let inner = Problem::new(origin, counts, step).map_err(js_err)?;
        Ok(WasmProblem { inner })
    }

    #[wasm_bindgen(js_name = setRegularisation)]
    pub fn set_regularisation(&mut self, weights: WeightsObj) {
        let d = RegularisationWeights::default();
        self.inner.set_regularisation(RegularisationWeights {
            dxx: weights.dxx().unwrap_or(d.dxx),
            dyy: weights.dyy().unwrap_or(d.dyy),
            dzz: weights.dzz().unwrap_or(d.dzz),
            dxy: weights.dxy().unwrap_or(d.dxy),
            dxz: weights.dxz().unwrap_or(d.dxz),
            dyz: weights.dyz().unwrap_or(d.dyz),
        });
        let (dx, dy, dz) = (weights.dx(), weights.dy(), weights.dz());
        if dx.is_some() || dy.is_some() || dz.is_some() {
            let b = BorderRegularisationWeights::default();
            self.inner
                .set_border_regularisation(BorderRegularisationWeights {
                    dx: dx.unwrap_or(b.dx),
                    dy: dy.unwrap_or(b.dy),
                    dz: dz.unwrap_or(b.dz),
                });
        }
    }

    #[wasm_bindgen(js_name = setRidgeFactor)]
    pub fn set_ridge_factor(&mut self, ridge_factor: f64) -> Result<(), JsError> {
        self.inner.set_ridge_factor(ridge_factor).map_err(js_err)
    }

    #[wasm_bindgen(js_name = setInitialGuess)]
    pub fn set_initial_guess(&mut self, x0: &[f64]) -> Result<(), JsError> {
        self.inner.set_initial_guess(x0).map_err(js_err)
    }

    #[wasm_bindgen(js_name = addValue)]
    pub fn add_value(
        &mut self,
        points: &[f64],
        values: &[f64],
        weights: &[f64],
    ) -> Result<(), JsError> {
        self.inner
            .add_value(points, values, weights)
            .map_err(js_err)
    }

    #[wasm_bindgen(js_name = addGradient)]
    pub fn add_gradient(
        &mut self,
        points: &[f64],
        vectors: &[f64],
        weights: &[f64],
    ) -> Result<(), JsError> {
        self.inner
            .add_gradient(points, vectors, weights)
            .map_err(js_err)
    }

    #[wasm_bindgen(js_name = addNormal)]
    pub fn add_normal(
        &mut self,
        points: &[f64],
        vectors: &[f64],
        weights: &[f64],
    ) -> Result<(), JsError> {
        self.inner
            .add_normal(points, vectors, weights)
            .map_err(js_err)
    }

    #[wasm_bindgen(js_name = addTangent)]
    pub fn add_tangent(
        &mut self,
        points: &[f64],
        vectors: &[f64],
        weights: &[f64],
    ) -> Result<(), JsError> {
        self.inner
            .add_tangent(points, vectors, weights)
            .map_err(js_err)
    }

    #[wasm_bindgen(js_name = addInterface)]
    pub fn add_interface(
        &mut self,
        points: &[f64],
        ids: &[f64],
        weight: f64,
    ) -> Result<(), JsError> {
        self.inner
            .add_interface(points, ids, weight)
            .map_err(js_err)
    }

    #[wasm_bindgen(js_name = addValueBound)]
    pub fn add_value_bound(
        &mut self,
        points: &[f64],
        lower: &[f64],
        upper: &[f64],
        weights: &[f64],
    ) -> Result<(), JsError> {
        self.inner
            .add_value_bound(points, lower, upper, weights)
            .map_err(js_err)
    }

    #[wasm_bindgen(js_name = addPairBound)]
    pub fn add_pair_bound(
        &mut self,
        points_a: &[f64],
        points_b: &[f64],
        lower: &[f64],
        upper: &[f64],
        weights: &[f64],
    ) -> Result<(), JsError> {
        self.inner
            .add_pair_bound(points_a, points_b, lower, upper, weights)
            .map_err(js_err)
    }

    pub fn estimate(&self, settings: Option<SettingsObj>) -> Result<WasmEstimate, JsError> {
        let s = parse_settings(settings)?;
        let e = self.inner.estimate_memory(&s).map_err(js_err)?;
        Ok(WasmEstimate {
            total: e.total as f64,
            vectors: e.vectors as f64,
            constraints: e.constraints as f64,
            working: e.working as f64,
        })
    }

    /// Start a solve. This consumes the problem.
    pub fn start(self, settings: Option<SettingsObj>) -> Result<WasmSession, JsError> {
        let s = parse_settings(settings)?;
        let session = self.inner.start(s).map_err(js_err)?;
        Ok(WasmSession {
            session: Some(session),
            field: None,
            report: None,
        })
    }
}

/// A running solve (JavaScript name `Session`).
#[wasm_bindgen(js_name = Session)]
pub struct WasmSession {
    /// `Some` until the solve is done.
    session: Option<SolveSession>,
    /// `Some` after the solve is done.
    field: Option<Field>,
    /// `Some` after the solve is done.
    report: Option<SolveReport>,
}

/// Make the step result for a finished solve.
fn done_result(r: &SolveReport) -> WasmStepResult {
    WasmStepResult {
        done: true,
        stage: "done".into(),
        iteration: r.iterations as f64,
        residual: r.residual,
        violation: r.max_violation,
        converged: r.converged,
        peak_bytes: r.peak_bytes as f64,
        warnings: r.warnings.clone(),
    }
}

impl WasmSession {
    /// Run `f` on the field. Before the end of the solve, the field is the
    /// best result so far.
    fn with_field<R>(&self, f: impl FnOnce(&Field) -> Result<R, CoreError>) -> Result<R, JsError> {
        if let Some(field) = &self.field {
            return f(field).map_err(js_err);
        }
        match &self.session {
            Some(s) => {
                let (field, _) = s.best_so_far().map_err(js_err)?;
                f(&field).map_err(js_err)
            }
            None => Err(arg_err("session", "the session has no field")),
        }
    }
}

#[wasm_bindgen(js_class = Session)]
impl WasmSession {
    /// Run at most `max_iters` iterations.
    pub fn step(&mut self, max_iters: f64) -> Result<WasmStepResult, JsError> {
        if let Some(report) = &self.report {
            return Ok(done_result(report));
        }
        let max_iters = to_count("maxIters", max_iters)?;
        let Some(session) = self.session.as_mut() else {
            return Err(arg_err("session", "the session has no field"));
        };
        match session.step(max_iters).map_err(js_err)? {
            Step::Running(p) => Ok(WasmStepResult {
                done: false,
                stage: match p.stage {
                    Stage::Setup => "setup",
                    Stage::Solve => "solve",
                    Stage::CheckConstraints => "check",
                }
                .into(),
                iteration: p.iteration as f64,
                residual: p.residual,
                violation: p.max_violation,
                converged: false,
                peak_bytes: 0.0,
                warnings: Vec::new(),
            }),
            Step::Done(report) => {
                // Keep the field. `finish` moves it out of the session.
                if let Some(session) = self.session.take() {
                    let (field, _) = session.finish().map_err(js_err)?;
                    self.field = Some(field);
                }
                let result = done_result(&report);
                self.report = Some(report);
                Ok(result)
            }
        }
    }

    /// Copy the field out: `Float64Array` or `Float32Array`.
    pub fn field(&self) -> Result<JsValue, JsError> {
        self.with_field(|field| {
            Ok(match &field.data {
                FieldData::F64(v) => JsValue::from(v.clone()),
                FieldData::F32(v) => JsValue::from(v.clone()),
            })
        })
    }

    pub fn evaluate(&self, points: &[f64]) -> Result<Vec<f64>, JsError> {
        self.with_field(|f| f.evaluate(points))
    }

    #[wasm_bindgen(js_name = evaluateGradient)]
    pub fn evaluate_gradient(&self, points: &[f64]) -> Result<Vec<f64>, JsError> {
        self.with_field(|f| f.evaluate_gradient(points))
    }
}

/// Size of the wasm linear memory in bytes. It is 0 on other targets.
/// The memory never shrinks, so after a solve this is the peak.
#[wasm_bindgen(js_name = memoryBytes)]
pub fn memory_bytes() -> f64 {
    #[cfg(target_arch = "wasm32")]
    {
        (core::arch::wasm32::memory_size::<0>() * 65536) as f64
    }
    #[cfg(not(target_arch = "wasm32"))]
    {
        0.0
    }
}
