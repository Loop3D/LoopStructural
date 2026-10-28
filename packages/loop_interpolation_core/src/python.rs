//! Python bindings (PyO3). Compiled only with the feature `python`.
//!
//! The module `loop_interpolation_core` has one class, `Problem`.
//! Arrays come in as contiguous `float64` numpy arrays.
//! `solve` releases the GIL and consumes the problem.

use numpy::{PyArray1, PyReadonlyArray1};
use pyo3::exceptions::{PyMemoryError, PyNotImplementedError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use crate::error::{BoundKind, CoreError};
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

#[pymodule]
fn loop_interpolation_core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyProblem>()?;
    Ok(())
}
