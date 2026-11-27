//! Python bindings (PyO3). Compiled only with the feature `python`.
//!
//! The module `loop_interpolation_core` has one class, `Problem`.
//! Arrays come in as contiguous `float64` numpy arrays.
//! `solve` releases the GIL and consumes the problem.

use numpy::{PyArray1, PyArrayMethods, PyReadonlyArray1, PyReadonlyArray2, PyUntypedArrayMethods};
use pyo3::exceptions::{PyMemoryError, PyNotImplementedError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use crate::error::{BoundKind, CoreError};
use crate::fault::{
    BlockSettings, Extension, Fault, FaultChain, GridField, Profile, Profile1, Region, Side,
};
use crate::field::FieldData;
use crate::interpolator::{BorderRegularisationWeights, RegularisationWeights};
use crate::problem::Problem;
use crate::report::SolveReport;
use crate::settings::{Precision, Settings};

/// Make a Python error from a core error. The text of the error stays.
fn map_err(e: CoreError) -> PyErr {
    let text = e.to_string();
    match e {
        CoreError::InvalidInput { .. }
        | CoreError::GridTooSmall { .. }
        | CoreError::GridTooLarge { .. }
        | CoreError::NoConstraints
        | CoreError::PointOutsideGrid { .. } => PyValueError::new_err(text),
        CoreError::OutOfMemory { .. } | CoreError::OverBudget { .. } => {
            PyMemoryError::new_err(text)
        }
        CoreError::Infeasible { worst } => {
            let rows: Vec<String> = worst
                .iter()
                .map(|w| {
                    let kind = match w.kind {
                        BoundKind::Value => "value",
                        BoundKind::Pair => "pair",
                    };
                    format!("{kind}[{}] violation {:e}", w.index, w.violation)
                })
                .collect();
            PyRuntimeError::new_err(format!("{text}: {}", rows.join("; ")))
        }
        CoreError::SolveRunning => PyRuntimeError::new_err(text),
        CoreError::Unsupported { .. } => PyNotImplementedError::new_err(text),
    }
}

/// Get a contiguous slice of an array, or raise `ValueError`.
fn slice<'a>(name: &str, a: &'a PyReadonlyArray1<'_, f64>) -> PyResult<&'a [f64]> {
    a.as_slice()
        .map_err(|_| PyValueError::new_err(format!("array '{name}' must be contiguous")))
}

/// Read a number from a dict value.
fn num(key: &str, v: &Bound<'_, PyAny>) -> PyResult<f64> {
    v.extract::<f64>()
        .map_err(|_| PyValueError::new_err(format!("value of '{key}' must be a number")))
}

/// Read an unsigned integer from a dict value.
fn uint(key: &str, v: &Bound<'_, PyAny>) -> PyResult<u64> {
    v.extract::<u64>().map_err(|_| {
        PyValueError::new_err(format!("value of '{key}' must be a non-negative integer"))
    })
}

/// Read a settings dict. An unknown key raises `ValueError`.
fn parse_settings(d: &Bound<'_, PyDict>) -> PyResult<Settings> {
    let mut s = Settings::default();
    for (k, v) in d.iter() {
        let key: String = k
            .extract()
            .map_err(|_| PyValueError::new_err("settings keys must be text"))?;
        match key.as_str() {
            "precision" => {
                let p: String = v
                    .extract()
                    .map_err(|_| PyValueError::new_err("'precision' must be \"f64\" or \"f32\""))?;
                s.precision = match p.as_str() {
                    "f64" => Precision::F64,
                    "f32" => Precision::F32,
                    other => {
                        return Err(PyValueError::new_err(format!(
                            "unknown precision '{other}'; use \"f64\" or \"f32\""
                        )));
                    }
                };
            }
            "tol" => s.tol = num(&key, &v)?,
            "max_iter" => s.max_iter = uint(&key, &v)? as usize,
            "budget_bytes" => {
                s.budget_bytes = if v.is_none() {
                    None
                } else {
                    Some(uint(&key, &v)?)
                }
            }
            "rho" => s.rho = num(&key, &v)?,
            "inner_max_iter" => s.inner_max_iter = uint(&key, &v)? as usize,
            "inner_tol" => s.inner_tol = num(&key, &v)?,
            "admm_abs_tol" => s.admm_abs_tol = num(&key, &v)?,
            "admm_rel_tol" => s.admm_rel_tol = num(&key, &v)?,
            other => {
                return Err(PyValueError::new_err(format!("unknown setting '{other}'")));
            }
        }
    }
    Ok(s)
}

/// Make a Python dict from a solve report.
fn report_dict<'py>(py: Python<'py>, r: &SolveReport) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("converged", r.converged)?;
    d.set_item("iterations", r.iterations)?;
    d.set_item("residual", r.residual)?;
    d.set_item("max_violation", r.max_violation)?;
    d.set_item("peak_bytes", r.peak_bytes)?;
    d.set_item("warnings", PyList::new(py, &r.warnings)?)?;
    Ok(d)
}

/// Read the weights of a dict. A key that is not in `allowed` raises `ValueError`.
/// Return the value for each allowed key, or `None` if the key is not set.
fn weights_from_dict(d: &Bound<'_, PyDict>, allowed: &[&str]) -> PyResult<Vec<Option<f64>>> {
    let mut out = vec![None; allowed.len()];
    for (k, v) in d.iter() {
        let key: String = k
            .extract()
            .map_err(|_| PyValueError::new_err("keys must be text"))?;
        match allowed.iter().position(|a| *a == key) {
            Some(i) => out[i] = Some(num(&key, &v)?),
            None => {
                return Err(PyValueError::new_err(format!(
                    "unknown key '{key}'; allowed keys: {}",
                    allowed.join(", ")
                )));
            }
        }
    }
    Ok(out)
}

/// Read a list of three numbers.
fn triple(name: &str, v: &Bound<'_, PyAny>) -> PyResult<[f64; 3]> {
    let x: Vec<f64> = v
        .extract()
        .map_err(|_| PyValueError::new_err(format!("'{name}' must have 3 numbers")))?;
    <[f64; 3]>::try_from(x)
        .map_err(|_| PyValueError::new_err(format!("'{name}' must have 3 numbers")))
}

/// A low-memory interpolation problem on a regular 3D grid.
#[pyclass(name = "Problem")]
struct PyProblem {
    /// `None` after `solve`. A problem can be solved one time.
    inner: Option<Problem>,
}

impl PyProblem {
    fn get(&mut self) -> PyResult<&mut Problem> {
        self.inner
            .as_mut()
            .ok_or_else(|| PyRuntimeError::new_err("problem already solved"))
    }
}

#[pymethods]
impl PyProblem {
    /// Make a problem. `origin`, `counts` and `step` each have 3 values.
    #[new]
    fn new(
        origin: &Bound<'_, PyAny>,
        counts: &Bound<'_, PyAny>,
        step: &Bound<'_, PyAny>,
    ) -> PyResult<Self> {
        let origin = triple("origin", origin)?;
        let step = triple("step", step)?;
        let counts: Vec<usize> = counts
            .extract()
            .map_err(|_| PyValueError::new_err("'counts' must have 3 positive integers"))?;
        let counts = <[usize; 3]>::try_from(counts)
            .map_err(|_| PyValueError::new_err("'counts' must have 3 positive integers"))?;
        let p = Problem::new(origin, counts, step).map_err(map_err)?;
        Ok(Self { inner: Some(p) })
    }

    /// Set the regularisation weights. Keys: dxx dyy dzz dxy dxz dyz.
    fn set_regularisation(&mut self, weights: &Bound<'_, PyDict>) -> PyResult<()> {
        let v = weights_from_dict(weights, &["dxx", "dyy", "dzz", "dxy", "dxz", "dyz"])?;
        let mut w = RegularisationWeights::default();
        if let Some(x) = v[0] {
            w.dxx = x;
        }
        if let Some(x) = v[1] {
            w.dyy = x;
        }
        if let Some(x) = v[2] {
            w.dzz = x;
        }
        if let Some(x) = v[3] {
            w.dxy = x;
        }
        if let Some(x) = v[4] {
            w.dxz = x;
        }
        if let Some(x) = v[5] {
            w.dyz = x;
        }
        self.get()?.set_regularisation(w);
        Ok(())
    }

    /// Set the border regularisation weights. Keys: dx dy dz.
    fn set_border_regularisation(&mut self, weights: &Bound<'_, PyDict>) -> PyResult<()> {
        let v = weights_from_dict(weights, &["dx", "dy", "dz"])?;
        let mut w = BorderRegularisationWeights::default();
        if let Some(x) = v[0] {
            w.dx = x;
        }
        if let Some(x) = v[1] {
            w.dy = x;
        }
        if let Some(x) = v[2] {
            w.dz = x;
        }
        self.get()?.set_border_regularisation(w);
        Ok(())
    }

    /// Set the ridge factor. 0 removes the ridge term. The default is 1e-8.
    fn set_ridge_factor(&mut self, ridge_factor: f64) -> PyResult<()> {
        self.get()?.set_ridge_factor(ridge_factor).map_err(map_err)
    }

    /// Set the start value of the solve. `x0` has one value for each node.
    fn set_start(&mut self, x0: PyReadonlyArray1<'_, f64>) -> PyResult<()> {
        let x = slice("x0", &x0)?;
        self.get()?.set_start(x).map_err(map_err)
    }

    /// Add value constraints.
    fn add_value(
        &mut self,
        points: PyReadonlyArray1<'_, f64>,
        values: PyReadonlyArray1<'_, f64>,
        weights: PyReadonlyArray1<'_, f64>,
    ) -> PyResult<()> {
        let (p, v, w) = (
            slice("points", &points)?,
            slice("values", &values)?,
            slice("weights", &weights)?,
        );
        self.get()?.add_value(p, v, w).map_err(map_err)
    }

    /// Add gradient constraints (direction only).
    fn add_gradient(
        &mut self,
        points: PyReadonlyArray1<'_, f64>,
        vectors: PyReadonlyArray1<'_, f64>,
        weights: PyReadonlyArray1<'_, f64>,
    ) -> PyResult<()> {
        let (p, v, w) = (
            slice("points", &points)?,
            slice("vectors", &vectors)?,
            slice("weights", &weights)?,
        );
        self.get()?.add_gradient(p, v, w).map_err(map_err)
    }

    /// Add normal constraints (the vector length counts).
    fn add_normal(
        &mut self,
        points: PyReadonlyArray1<'_, f64>,
        vectors: PyReadonlyArray1<'_, f64>,
        weights: PyReadonlyArray1<'_, f64>,
    ) -> PyResult<()> {
        let (p, v, w) = (
            slice("points", &points)?,
            slice("vectors", &vectors)?,
            slice("weights", &weights)?,
        );
        self.get()?.add_normal(p, v, w).map_err(map_err)
    }

    /// Add tangent constraints.
    fn add_tangent(
        &mut self,
        points: PyReadonlyArray1<'_, f64>,
        vectors: PyReadonlyArray1<'_, f64>,
        weights: PyReadonlyArray1<'_, f64>,
    ) -> PyResult<()> {
        let (p, v, w) = (
            slice("points", &points)?,
            slice("vectors", &vectors)?,
            slice("weights", &weights)?,
        );
        self.get()?.add_tangent(p, v, w).map_err(map_err)
    }

    /// Add interface constraints. Points with the same id have the same value.
    fn add_interface(
        &mut self,
        points: PyReadonlyArray1<'_, f64>,
        ids: PyReadonlyArray1<'_, f64>,
        weight: f64,
    ) -> PyResult<()> {
        let (p, i) = (slice("points", &points)?, slice("ids", &ids)?);
        self.get()?.add_interface(p, i, weight).map_err(map_err)
    }

    /// Add directional regularisation. `directions` has 3 values for each node.
    fn add_directional_regularisation(
        &mut self,
        directions: PyReadonlyArray1<'_, f64>,
        weight: f64,
    ) -> PyResult<()> {
        let d = slice("directions", &directions)?;
        self.get()?
            .add_directional_regularisation(d, weight)
            .map_err(map_err)
    }

    /// Add value bounds. A bound of -inf or +inf means no bound.
    fn add_value_bound(
        &mut self,
        points: PyReadonlyArray1<'_, f64>,
        lower: PyReadonlyArray1<'_, f64>,
        upper: PyReadonlyArray1<'_, f64>,
        weights: PyReadonlyArray1<'_, f64>,
    ) -> PyResult<()> {
        let (p, l, u, w) = (
            slice("points", &points)?,
            slice("lower", &lower)?,
            slice("upper", &upper)?,
            slice("weights", &weights)?,
        );
        self.get()?.add_value_bound(p, l, u, w).map_err(map_err)
    }

    /// Add pair bounds on the difference of the field at two points.
    fn add_pair_bound(
        &mut self,
        points_a: PyReadonlyArray1<'_, f64>,
        points_b: PyReadonlyArray1<'_, f64>,
        lower: PyReadonlyArray1<'_, f64>,
        upper: PyReadonlyArray1<'_, f64>,
        weights: PyReadonlyArray1<'_, f64>,
    ) -> PyResult<()> {
        let (a, b, l, u, w) = (
            slice("points_a", &points_a)?,
            slice("points_b", &points_b)?,
            slice("lower", &lower)?,
            slice("upper", &upper)?,
            slice("weights", &weights)?,
        );
        self.get()?.add_pair_bound(a, b, l, u, w).map_err(map_err)
    }

    /// Estimate the memory of a solve. Return a dict with the keys
    /// vectors, constraints, working and total (bytes).
    fn estimate<'py>(
        &mut self,
        py: Python<'py>,
        settings: &Bound<'py, PyDict>,
    ) -> PyResult<Bound<'py, PyDict>> {
        let s = parse_settings(settings)?;
        let e = self.get()?.estimate_memory(&s).map_err(map_err)?;
        let d = PyDict::new(py);
        d.set_item("vectors", e.vectors)?;
        d.set_item("constraints", e.constraints)?;
        d.set_item("working", e.working)?;
        d.set_item("total", e.total)?;
        Ok(d)
    }

    /// Solve the problem. Return the field (numpy array) and a report dict.
    /// The problem is used up. A second call raises `RuntimeError`.
    fn solve<'py>(
        &mut self,
        py: Python<'py>,
        settings: &Bound<'py, PyDict>,
    ) -> PyResult<(Bound<'py, PyAny>, Bound<'py, PyDict>)> {
        let s = parse_settings(settings)?;
        let problem = self
            .inner
            .take()
            .ok_or_else(|| PyRuntimeError::new_err("problem already solved"))?;
        // The solve does not use Python objects. Release the GIL.
        let (field, report) = py.detach(move || problem.solve(s)).map_err(map_err)?;
        let array: Bound<'py, PyAny> = match field.data {
            FieldData::F64(v) => PyArray1::from_vec(py, v).into_any(),
            FieldData::F32(v) => PyArray1::from_vec(py, v).into_any(),
        };
        Ok((array, report_dict(py, &report)?))
    }
}

// ---------------------------------------------------------------------------
// Fault kernel
// ---------------------------------------------------------------------------

/// Make a `ValueError`.
fn bad(ctx: &str, msg: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(format!("{ctx}: {msg}"))
}

/// Cast a value to a dict, or raise `ValueError`.
fn as_dict<'a, 'py>(ctx: &str, v: &'a Bound<'py, PyAny>) -> PyResult<&'a Bound<'py, PyDict>> {
    v.cast::<PyDict>().map_err(|_| bad(ctx, "must be a dict"))
}

/// Get a key of a dict, or raise `ValueError`. A key that is not in `allowed` also raises.
fn fields<'py>(
    ctx: &str,
    d: &Bound<'py, PyDict>,
    allowed: &[&str],
) -> PyResult<Vec<Option<Bound<'py, PyAny>>>> {
    for k in d.keys() {
        let key: String = k.extract().map_err(|_| bad(ctx, "keys must be text"))?;
        if !allowed.contains(&key.as_str()) {
            return Err(bad(ctx, format!("unknown key '{key}'")));
        }
    }
    allowed.iter().map(|k| d.get_item(*k)).collect()
}

/// Get a required value from the result of `fields`.
fn need<'py>(ctx: &str, key: &str, v: &Option<Bound<'py, PyAny>>) -> PyResult<Bound<'py, PyAny>> {
    v.clone()
        .ok_or_else(|| bad(ctx, format!("missing key '{key}'")))
}

/// Read a float, a list of floats or a numpy array of floats.
fn floats(ctx: &str, key: &str, v: &Bound<'_, PyAny>) -> PyResult<Vec<f64>> {
    if let Ok(a) = v.extract::<PyReadonlyArray1<'_, f64>>() {
        return Ok(a.as_array().to_vec());
    }
    v.extract::<Vec<f64>>()
        .map_err(|_| bad(ctx, format!("'{key}' must be an array of numbers")))
}

/// Read a fixed number of floats.
fn fixed<const N: usize>(ctx: &str, key: &str, v: &Bound<'_, PyAny>) -> PyResult<[f64; N]> {
    <[f64; N]>::try_from(floats(ctx, key, v)?)
        .map_err(|_| bad(ctx, format!("'{key}' must have {N} numbers")))
}

fn float(ctx: &str, key: &str, v: &Bound<'_, PyAny>) -> PyResult<f64> {
    v.extract::<f64>()
        .map_err(|_| bad(ctx, format!("'{key}' must be a number")))
}

fn index(ctx: &str, key: &str, v: &Bound<'_, PyAny>, nfaults: usize) -> PyResult<usize> {
    match v.extract::<usize>() {
        Ok(i) if i < nfaults => Ok(i),
        _ => Err(bad(
            ctx,
            format!("'{key}' must be an integer below {nfaults}"),
        )),
    }
}

fn indices(ctx: &str, key: &str, v: &Bound<'_, PyAny>, nfaults: usize) -> PyResult<Vec<usize>> {
    let items: Vec<Bound<'_, PyAny>> = v
        .extract()
        .map_err(|_| bad(ctx, format!("'{key}' must be a list of integers")))?;
    items.iter().map(|i| index(ctx, key, i, nfaults)).collect()
}

/// Read a grid dict.
fn parse_grid(ctx: &str, v: &Bound<'_, PyAny>) -> PyResult<GridField> {
    let d = fields(
        ctx,
        as_dict(ctx, v)?,
        &[
            "origin",
            "nsteps",
            "step_vector",
            "local_origin",
            "local_rotation",
            "values",
        ],
    )?;
    let k = |i: usize, key: &str| need(ctx, key, &d[i]);
    let origin = fixed::<3>(ctx, "origin", &k(0, "origin")?)?;
    let nsteps: Vec<usize> = k(1, "nsteps")?
        .extract()
        .map_err(|_| bad(ctx, "'nsteps' must have 3 positive integers"))?;
    let nsteps = <[usize; 3]>::try_from(nsteps)
        .ok()
        .filter(|n| n.iter().all(|&c| c >= 2))
        .ok_or_else(|| bad(ctx, "'nsteps' must have 3 integers of 2 or more"))?;
    let step_vector = fixed::<3>(ctx, "step_vector", &k(2, "step_vector")?)?;
    let local_origin = match &d[3] {
        Some(v) if !v.is_none() => Some(fixed::<3>(ctx, "local_origin", v)?),
        _ => None,
    };
    let local_rotation = match &d[4] {
        Some(v) if !v.is_none() => {
            let flat: Vec<f64> = if let Ok(a) = v.extract::<PyReadonlyArray2<'_, f64>>() {
                a.as_array().iter().copied().collect()
            } else {
                let rows: Vec<Vec<f64>> = v
                    .extract()
                    .map_err(|_| bad(ctx, "'local_rotation' must be a 3x3 matrix"))?;
                rows.concat()
            };
            if flat.len() != 9 {
                return Err(bad(ctx, "'local_rotation' must be a 3x3 matrix"));
            }
            Some(std::array::from_fn(|i| {
                std::array::from_fn(|j| flat[3 * i + j])
            }))
        }
        _ => None,
    };
    let values = floats(ctx, "values", &k(5, "values")?)?;
    if values.len() != nsteps.iter().product::<usize>() {
        return Err(bad(
            ctx,
            format!(
                "'values' has {} numbers; 'nsteps' needs {}",
                values.len(),
                nsteps.iter().product::<usize>()
            ),
        ));
    }
    Ok(GridField {
        origin,
        nsteps,
        step_vector,
        local_origin,
        local_rotation,
        values,
    })
}

/// Read a profile part.
fn parse_part(ctx: &str, v: &Bound<'_, PyAny>) -> PyResult<Profile1> {
    let dict = as_dict(ctx, v)?;
    let kind: String = need(ctx, "kind", &dict.get_item("kind")?)?
        .extract()
        .map_err(|_| bad(ctx, "'kind' must be text"))?;
    let allowed: &[&str] = match kind.as_str() {
        "ones" | "zeros" | "smooth_peak" => &["kind"],
        "cubic" => &["kind", "w", "lim", "min", "max", "value_points"],
        "composite" => &["kind", "positive", "negative"],
        other => return Err(bad(ctx, format!("unknown profile kind '{other}'"))),
    };
    let d = fields(ctx, dict, allowed)?;
    let k = |i: usize| need(ctx, allowed[i], &d[i]);
    Ok(match kind.as_str() {
        "ones" => Profile1::Ones,
        "zeros" => Profile1::Zeros,
        "smooth_peak" => Profile1::SmoothPeak,
        "cubic" => {
            let points: Vec<[f64; 2]> = k(5)?
                .extract()
                .map_err(|_| bad(ctx, "'value_points' must be a list of [x, y]"))?;
            Profile1::Cubic {
                w: fixed::<4>(ctx, "w", &k(1)?)?,
                lim: fixed::<2>(ctx, "lim", &k(2)?)?,
                min: float(ctx, "min", &k(3)?)?,
                max: float(ctx, "max", &k(4)?)?,
                value_points: points.iter().map(|p| (p[0], p[1])).collect(),
            }
        }
        _ => Profile1::Composite {
            positive: Box::new(parse_part(ctx, &k(1)?)?),
            negative: Box::new(parse_part(ctx, &k(2)?)?),
        },
    })
}

/// Read a region dict.
fn parse_region(ctx: &str, v: &Bound<'_, PyAny>, nfaults: usize) -> PyResult<Region> {
    let dict = as_dict(ctx, v)?;
    let kind: String = need(ctx, "kind", &dict.get_item("kind")?)?
        .extract()
        .map_err(|_| bad(ctx, "'kind' must be text"))?;
    match kind.as_str() {
        "positive" | "negative" => {
            let d = fields(ctx, dict, &["kind", "fault", "vector", "point"])?;
            let opt = |i: usize, key: &str| match &d[i] {
                Some(v) if !v.is_none() => fixed::<3>(ctx, key, v).map(Some),
                _ => Ok(None),
            };
            Ok(Region::Sign {
                positive: kind == "positive",
                fault: index(ctx, "fault", &need(ctx, "fault", &d[1])?, nfaults)?,
                vector: opt(2, "vector")?,
                point: opt(3, "point")?,
            })
        }
        "unconformity" => {
            let d = fields(ctx, dict, &["kind", "field", "faults", "value", "sign"])?;
            Ok(Region::Unconformity {
                field: parse_grid(ctx, &need(ctx, "field", &d[1])?)?,
                faults: indices(ctx, "faults", &need(ctx, "faults", &d[2])?, nfaults)?,
                value: float(ctx, "value", &need(ctx, "value", &d[3])?)?,
                sign: need(ctx, "sign", &d[4])?
                    .extract()
                    .map_err(|_| bad(ctx, "'sign' must be a bool"))?,
            })
        }
        other => Err(bad(ctx, format!("unknown region kind '{other}'"))),
    }
}

/// Read a fault dict.
fn parse_fault(i: usize, v: &Bound<'_, PyAny>, nfaults: usize) -> PyResult<Fault> {
    let probe = format!("faults[{i}]");
    let d = fields(
        &probe,
        as_dict(&probe, v)?,
        &[
            "name",
            "displacement",
            "steps",
            "frame",
            "cutters",
            "regions",
            "profile",
        ],
    )?;
    // The name is only for error messages.
    let name: String = match &d[0] {
        Some(n) => n
            .extract()
            .map_err(|_| bad(&probe, "'name' must be text"))?,
        None => String::new(),
    };
    let ctx = format!("fault {i} '{name}'");
    let k = |j: usize, key: &str| need(&ctx, key, &d[j]);
    let frame: Vec<Bound<'_, PyAny>> = k(3, "frame")?
        .extract()
        .map_err(|_| bad(&ctx, "'frame' must be a list of 3 grids"))?;
    if frame.len() != 3 {
        return Err(bad(&ctx, "'frame' must be a list of 3 grids"));
    }
    let grid = |j: usize| parse_grid(&format!("{ctx} key 'frame'[{j}]"), &frame[j]);
    let frame = [grid(0)?, grid(1)?, grid(2)?];
    let rctx = format!("{ctx} key 'regions'");
    let lists: Vec<Vec<Bound<'_, PyAny>>> = k(5, "regions")?
        .extract()
        .map_err(|_| bad(&rctx, "must be a list of 3 lists"))?;
    if lists.len() != 3 {
        return Err(bad(&rctx, "must be a list of 3 lists"));
    }
    let region = |j: usize| -> PyResult<Vec<Region>> {
        lists[j]
            .iter()
            .map(|r| parse_region(&rctx, r, nfaults))
            .collect()
    };
    let pctx = format!("{ctx} key 'profile'");
    let pd = fields(
        &pctx,
        as_dict(&pctx, &k(6, "profile")?)?,
        &["scale", "gx", "gy", "gz"],
    )?;
    let part = |j: usize, key: &str| parse_part(&pctx, &need(&pctx, key, &pd[j])?);
    let steps = k(2, "steps")?
        .extract::<usize>()
        .ok()
        .filter(|&s| s > 0)
        .ok_or_else(|| bad(&ctx, "'steps' must be a positive integer"))?;
    Ok(Fault {
        frame,
        cutters: indices(&ctx, "cutters", &k(4, "cutters")?, nfaults)?,
        regions: [region(0)?, region(1)?, region(2)?],
        profile: Profile {
            scale: float(&pctx, "scale", &need(&pctx, "scale", &pd[0])?)?,
            gx: part(1, "gx")?,
            gy: part(2, "gy")?,
            gz: part(3, "gz")?,
        },
        displacement: float(&ctx, "displacement", &k(1, "displacement")?)?,
        steps,
    })
}

/// Read the spec of a fault chain.
fn parse_chain(spec: &Bound<'_, PyAny>) -> PyResult<FaultChain> {
    let d = fields(
        "spec",
        as_dict("spec", spec)?,
        &["version", "block", "faults"],
    )?;
    let version = need("spec", "version", &d[0])?
        .extract::<i64>()
        .map_err(|_| bad("spec", "'version' must be an integer"))?;
    if version != 1 {
        return Err(bad("spec", format!("unsupported version {version}; use 1")));
    }
    let bctx = "spec key 'block'";
    let b = fields(
        bctx,
        as_dict(bctx, &need("spec", "block", &d[1])?)?,
        &["side_epsilon", "extension", "exact", "limit", "step"],
    )?;
    let num = |i: usize, key: &str| float(bctx, key, &need(bctx, key, &b[i])?);
    let extension = match need(bctx, "extension", &b[1])?
        .extract::<String>()
        .as_deref()
    {
        Ok("smooth") => Extension::Smooth,
        Ok("mirror") => Extension::Mirror,
        _ => return Err(bad(bctx, "'extension' must be \"smooth\" or \"mirror\"")),
    };
    let block = BlockSettings {
        side_epsilon: num(0, "side_epsilon")?,
        extension,
        exact: num(2, "exact")?,
        limit: num(3, "limit")?,
        step: num(4, "step")?,
    };
    let items: Vec<Bound<'_, PyAny>> = need("spec", "faults", &d[2])?
        .extract()
        .map_err(|_| bad("spec", "'faults' must be a list"))?;
    let faults = items
        .iter()
        .enumerate()
        .map(|(i, f)| parse_fault(i, f, items.len()))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(FaultChain { faults, block })
}

/// Copy an (N, 3) float64 array.
fn read_points(points: &PyReadonlyArray2<'_, f64>) -> PyResult<Vec<[f64; 3]>> {
    if points.shape()[1] != 3 {
        return Err(PyValueError::new_err("'points' must have shape (N, 3)"));
    }
    Ok(points
        .as_array()
        .rows()
        .into_iter()
        .map(|r| [r[0], r[1], r[2]])
        .collect())
}

/// Run `work` without the GIL. A panic becomes `RuntimeError`, with no partial result.
fn guarded<T: Send>(py: Python<'_>, work: impl FnOnce() -> T + Send) -> PyResult<T> {
    py.detach(|| std::panic::catch_unwind(std::panic::AssertUnwindSafe(work)))
        .map_err(|e| {
            let msg = e
                .downcast_ref::<&str>()
                .map(|s| s.to_string())
                .or_else(|| e.downcast_ref::<String>().cloned())
                .unwrap_or_else(|| "unknown error".into());
            PyRuntimeError::new_err(format!("fault kernel error: {msg}"))
        })
}

fn rows_array<'py>(py: Python<'py>, rows: Vec<[f64; 3]>) -> PyResult<Bound<'py, PyAny>> {
    let n = rows.len();
    let flat: Vec<f64> = rows.into_iter().flatten().collect();
    Ok(PyArray1::from_vec(py, flat).reshape([n, 3])?.into_any())
}

/// One value for each fault.
fn per_fault<T: Copy>(name: &str, a: &[T], n: usize) -> PyResult<()> {
    if a.len() == n {
        Ok(())
    } else {
        Err(PyValueError::new_err(format!(
            "'{name}' needs {n} values (one for each fault), got {}",
            a.len()
        )))
    }
}

fn read_sides(sides: &PyReadonlyArray1<'_, i8>, n: usize) -> PyResult<Vec<Side>> {
    let s = sides.as_array().to_vec();
    per_fault("sides", &s, n)?;
    s.iter()
        .map(|&v| match v {
            2 => Ok(Side::Free),
            0 => Ok(Side::Absent),
            1 => Ok(Side::Hanging),
            -1 => Ok(Side::Foot),
            other => Err(PyValueError::new_err(format!(
                "'sides' has value {other}; use 2, 0, 1 or -1"
            ))),
        })
        .collect()
}

fn read_restored(restored: &PyReadonlyArray1<'_, bool>, n: usize) -> PyResult<Vec<bool>> {
    let r = restored.as_array().to_vec();
    per_fault("restored", &r, n)?;
    Ok(r)
}

/// A set of faults. It is read only after it is made.
#[pyclass(name = "FaultChain", frozen)]
struct PyFaultChain {
    inner: FaultChain,
}

#[pymethods]
impl PyFaultChain {
    /// Make a chain from a spec dict. A bad spec raises `ValueError`.
    #[new]
    fn new(spec: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self {
            inner: parse_chain(spec)?,
        })
    }

    /// Move the points through the faults in `order`. Return an (N, 3) array.
    #[allow(clippy::too_many_arguments)] // the Python signature is fixed by the plan
    fn restore<'py>(
        &self,
        py: Python<'py>,
        points: PyReadonlyArray2<'py, f64>,
        order: PyReadonlyArray1<'py, i64>,
        reverse: bool,
        sides: PyReadonlyArray1<'py, i8>,
        restored: PyReadonlyArray1<'py, bool>,
        threads: usize,
    ) -> PyResult<Bound<'py, PyAny>> {
        let n = self.inner.faults.len();
        let pts = read_points(&points)?;
        let sides = read_sides(&sides, n)?;
        let restored = read_restored(&restored, n)?;
        let order = order
            .as_array()
            .iter()
            .map(|&i| {
                usize::try_from(i).ok().filter(|&i| i < n).ok_or_else(|| {
                    PyValueError::new_err(format!("'order' has index {i}; the spec has {n} faults"))
                })
            })
            .collect::<PyResult<Vec<usize>>>()?;
        let out = guarded(py, || {
            self.inner
                .restore_batch(&order, &pts, reverse, &sides, &restored, threads)
        })?;
        rows_array(py, out)
    }

    /// Return one frame coordinate of a fault at the points. Return an (N,) array.
    #[allow(clippy::too_many_arguments)] // the Python signature is fixed by the plan
    fn frame_value<'py>(
        &self,
        py: Python<'py>,
        fault: usize,
        coordinate: usize,
        points: PyReadonlyArray2<'py, f64>,
        sides: PyReadonlyArray1<'py, i8>,
        restored: PyReadonlyArray1<'py, bool>,
        threads: usize,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let n = self.inner.faults.len();
        if fault >= n {
            return Err(PyValueError::new_err(format!(
                "'fault' is {fault}; the spec has {n} faults"
            )));
        }
        if coordinate >= 3 {
            return Err(PyValueError::new_err("'coordinate' must be 0, 1 or 2"));
        }
        let pts = read_points(&points)?;
        let sides = read_sides(&sides, n)?;
        let restored = read_restored(&restored, n)?;
        let out = guarded(py, || {
            self.inner
                .frame_value_batch(fault, coordinate, &pts, &sides, &restored, threads)
        })?;
        Ok(PyArray1::from_vec(py, out))
    }
}

/// A scalar field on a regular grid. It is read only after it is made.
#[pyclass(name = "GridField", frozen)]
struct PyGridField {
    inner: GridField,
}

impl PyGridField {
    /// Apply `f` to each point. Use the pool with the feature `parallel`.
    fn map<T: Send>(
        &self,
        points: &[[f64; 3]],
        threads: usize,
        f: impl Fn(&GridField, [f64; 3]) -> T + Sync,
    ) -> Vec<T> {
        #[cfg(feature = "parallel")]
        {
            use rayon::prelude::*;
            crate::fault::pool::pool(threads)
                .install(|| points.par_iter().map(|&p| f(&self.inner, p)).collect())
        }
        #[cfg(not(feature = "parallel"))]
        {
            let _ = threads;
            points.iter().map(|&p| f(&self.inner, p)).collect()
        }
    }
}

#[pymethods]
impl PyGridField {
    /// Make a field from a grid dict. A bad spec raises `ValueError`.
    #[new]
    fn new(spec: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self {
            inner: parse_grid("spec", spec)?,
        })
    }

    /// Return the field values at the points. Return an (N,) array.
    fn values<'py>(
        &self,
        py: Python<'py>,
        points: PyReadonlyArray2<'py, f64>,
        threads: usize,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let pts = read_points(&points)?;
        let out = guarded(py, || self.map(&pts, threads, |g, p| g.value(p)))?;
        Ok(PyArray1::from_vec(py, out))
    }

    /// Return the field gradients at the points. Return an (N, 3) array.
    fn gradients<'py>(
        &self,
        py: Python<'py>,
        points: PyReadonlyArray2<'py, f64>,
        threads: usize,
    ) -> PyResult<Bound<'py, PyAny>> {
        let pts = read_points(&points)?;
        let out = guarded(py, || self.map(&pts, threads, |g, p| g.gradient(p)))?;
        rows_array(py, out)
    }
}

#[pymodule]
fn loop_interpolation_core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyProblem>()?;
    m.add_class::<PyFaultChain>()?;
    m.add_class::<PyGridField>()?;
    m.add("FAULT_KERNEL_VERSION", 1)?;
    Ok(())
}
