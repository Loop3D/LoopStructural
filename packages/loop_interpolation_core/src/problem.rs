//! The problem: grid, constraints and weights.

use std::collections::HashMap;
use std::fmt;

use crate::error::CoreError;
use crate::field::Field;
use crate::grid::StructuredGrid3D;
use crate::interpolator::{
    BorderRegularisationWeights, FiniteDifferenceInterpolator, RegularisationWeights,
};
use crate::memory::{ProblemCounts, estimate, max_grid_for_budget};
use crate::report::{MemoryEstimate, SolveReport, Step};
use crate::session::SolveSession;
use crate::settings::{Precision, Settings};

/// The most nodes that a grid can have: `isize::MAX / 8` (FR-018).
const MAX_NODES: u64 = isize::MAX as u64 / 8;

/// From this many value rows the conflict check is skipped.
const CONFLICT_CHECK_LIMIT: usize = 1_000_000;

/// Make the error for a bad row. `family` is a name such as "value".
fn bad_row(family: &str, row: usize, reason: &str) -> CoreError {
    CoreError::InvalidInput {
        constraint: format!("{family}[{row}]"),
        reason: reason.into(),
    }
}

/// Make the error for a bad array or a bad scalar. No row is named.
fn bad_input(family: &str, reason: &str) -> CoreError {
    CoreError::InvalidInput {
        constraint: family.into(),
        reason: reason.into(),
    }
}

/// Check that a flat array holds 3 values for each point. Return the count.
fn count_points(family: &str, name: &str, flat: &[f64]) -> Result<usize, CoreError> {
    if flat.len() % 3 != 0 {
        return Err(bad_input(
            family,
            &format!("{name} length {} is not a multiple of 3", flat.len()),
        ));
    }
    Ok(flat.len() / 3)
}

/// Check that an array has `n` values.
fn check_len(family: &str, name: &str, len: usize, n: usize) -> Result<(), CoreError> {
    if len != n {
        return Err(bad_input(
            family,
            &format!("{name} has {len} values; {n} are needed"),
        ));
    }
    Ok(())
}

/// Get row `i` of a flat array of 3-vectors.
fn row3(flat: &[f64], i: usize) -> [f64; 3] {
    [flat[3 * i], flat[3 * i + 1], flat[3 * i + 2]]
}

/// Check that all 3 values of row `i` are finite.
fn check_finite3(family: &str, flat: &[f64], i: usize, what: &str) -> Result<(), CoreError> {
    if row3(flat, i).iter().all(|v| v.is_finite()) {
        Ok(())
    } else {
        Err(bad_row(family, i, &format!("{what} is not finite")))
    }
}

/// Check that a weight is finite and more than zero.
fn check_weight(family: &str, row: usize, w: f64) -> Result<(), CoreError> {
    if w.is_finite() && w > 0.0 {
        Ok(())
    } else {
        Err(bad_row(
            family,
            row,
            "weight must be finite and more than 0",
        ))
    }
}

/// Check that a vector row is finite and not zero.
fn check_vector(family: &str, flat: &[f64], i: usize) -> Result<(), CoreError> {
    check_finite3(family, flat, i, "vector")?;
    let v = row3(flat, i);
    if v[0] * v[0] + v[1] * v[1] + v[2] * v[2] > 0.0 {
        Ok(())
    } else {
        Err(bad_row(family, i, "vector has zero length"))
    }
}

/// Check a pair of bounds. `lower` and `upper` can be infinite.
fn check_bounds(family: &str, row: usize, lower: f64, upper: f64) -> Result<(), CoreError> {
    if lower.is_nan() || upper.is_nan() {
        return Err(bad_row(family, row, "bound is NaN"));
    }
    if lower > upper {
        return Err(bad_row(family, row, "lower bound is more than upper bound"));
    }
    Ok(())
}

/// A flat-input interpolation problem. Arrays are flat: points have 3n values.
///
/// The problem holds the ported interpolator. The interpolator is not
/// `Clone`, so the problem is not `Clone`.
pub struct Problem {
    /// Always `Some`, except for a short time inside `with_interp`.
    interp: Option<FiniteDifferenceInterpolator>,
    pub(crate) warnings: Vec<String>,
    /// First value seen at each point (bits of the coordinates).
    /// `None` after the conflict check is skipped.
    seen_values: Option<HashMap<[u64; 3], f64>>,
    /// Number of value rows seen, also the dropped ones.
    value_rows: usize,
    /// True after `add_directional_regularisation`.
    uses_directional: bool,
}

impl fmt::Debug for Problem {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Problem")
            .field("grid", &self.interp().grid)
            .field("warnings", &self.warnings)
            .finish_non_exhaustive()
    }
}

impl Problem {
    fn interp(&self) -> &FiniteDifferenceInterpolator {
        self.interp.as_ref().expect("interpolator is present")
    }

    fn interp_mut(&mut self) -> &mut FiniteDifferenceInterpolator {
        self.interp.as_mut().expect("interpolator is present")
    }

    /// Change the interpolator with a method that takes it by value.
    fn with_interp(
        &mut self,
        f: impl FnOnce(FiniteDifferenceInterpolator) -> FiniteDifferenceInterpolator,
    ) {
        let old = self.interp.take().expect("interpolator is present");
        self.interp = Some(f(old));
    }

    /// Number of constraint rows of all families.
    fn n_constraints(&self) -> usize {
        let i = self.interp();
        i.n_value_constraints()
            + i.n_gradient_constraints()
            + i.n_norm_constraints()
            + i.n_tangent_constraints()
            + i.n_interface_constraints()
            + i.n_inequality_constraints()
    }

    /// Count the points that are outside the grid.
    fn count_outside(&self, points: &[f64]) -> usize {
        let grid = &self.interp().grid;
        (0..points.len() / 3)
            .filter(|&i| grid.trilinear_coeffs(row3(points, i)).is_none())
            .count()
    }

    /// Add one warning if `outside` is more than zero.
    fn warn_outside(&mut self, family: &str, outside: usize) {
        if outside > 0 {
            self.warnings.push(format!(
                "{family}: {outside} point(s) outside the grid were dropped"
            ));
        }
    }

    pub fn new(origin: [f64; 3], counts: [usize; 3], step: [f64; 3]) -> Result<Problem, CoreError> {
        // Check the sizes first. Do this before any large allocation.
        for (axis, &nodes) in counts.iter().enumerate() {
            if nodes < 3 {
                return Err(CoreError::GridTooSmall { axis, nodes });
            }
        }
        let nodes = (counts[0] as u128)
            .saturating_mul(counts[1] as u128)
            .saturating_mul(counts[2] as u128);
        if nodes > MAX_NODES as u128 {
            return Err(CoreError::GridTooLarge {
                nodes: u64::try_from(nodes).unwrap_or(u64::MAX),
                max_nodes: MAX_NODES,
            });
        }
        if !origin.iter().all(|v| v.is_finite()) {
            return Err(bad_input("grid", "origin is not finite"));
        }
        if !step.iter().all(|v| v.is_finite() && *v > 0.0) {
            return Err(bad_input("grid", "step must be finite and more than 0"));
        }
        let mut warnings = Vec::new();
        let max_step = step.iter().cloned().fold(f64::MIN, f64::max);
        let min_step = step.iter().cloned().fold(f64::MAX, f64::min);
        if max_step / min_step > 100.0 {
            warnings.push(format!(
                "grid step ratio {:.1} is more than 100; the result can be poor",
                max_step / min_step
            ));
        }
        // The ported constructor builds tables (about 7 bytes per node).
        // Task T018 makes this lazy.
        let interp = FiniteDifferenceInterpolator::new(StructuredGrid3D::new(origin, counts, step));
        Ok(Problem {
            interp: Some(interp),
            warnings,
            seen_values: Some(HashMap::new()),
            value_rows: 0,
            uses_directional: false,
        })
    }

    pub fn set_regularisation(&mut self, w: RegularisationWeights) {
        self.with_interp(|i| i.with_regularisation(w));
    }

    pub fn set_border_regularisation(&mut self, w: BorderRegularisationWeights) {
        self.with_interp(|i| i.with_border_regularisation(w));
    }

    pub fn add_value(
        &mut self,
        points: &[f64],
        values: &[f64],
        weights: &[f64],
    ) -> Result<(), CoreError> {
        const F: &str = "value";
        let n = count_points(F, "points", points)?;
        check_len(F, "values", values.len(), n)?;
        check_len(F, "weights", weights.len(), n)?;
        for i in 0..n {
            check_finite3(F, points, i, "point")?;
            if !values[i].is_finite() {
                return Err(bad_row(F, i, "value is not finite"));
            }
            check_weight(F, i, weights[i])?;
        }
        // Decide if the conflict check runs.
        if self.seen_values.is_some() && self.value_rows.saturating_add(n) >= CONFLICT_CHECK_LIMIT {
            self.seen_values = None;
            self.warnings.push(
                "value: 1000000 rows or more; the check for conflicting values was skipped".into(),
            );
        }
        self.value_rows = self.value_rows.saturating_add(n);
        let mut outside = 0usize;
        let mut conflicts = 0usize;
        // Upper bound: rows outside the grid are dropped later.
        self.interp_mut().reserve_value(n);
        for i in 0..n {
            let p = row3(points, i);
            if self.interp().grid.trilinear_coeffs(p).is_none() {
                outside += 1;
                continue;
            }
            if let Some(seen) = self.seen_values.as_mut() {
                // Adding 0.0 turns -0.0 into 0.0, so both give the same key.
                let key = [
                    (p[0] + 0.0).to_bits(),
                    (p[1] + 0.0).to_bits(),
                    (p[2] + 0.0).to_bits(),
                ];
                match seen.get(&key) {
                    Some(&v) if v != values[i] => conflicts += 1,
                    Some(_) => {}
                    None => {
                        seen.insert(key, values[i]);
                    }
                }
            }
            self.interp_mut()
                .add_value_constraint(p, values[i], weights[i]);
        }
        self.warn_outside(F, outside);
        if conflicts > 0 {
            self.warnings.push(format!(
                "value: {conflicts} row(s) have a different value from an earlier row at the same point"
            ));
        }
        Ok(())
    }

    /// Shared code for the gradient, normal and tangent families.
    fn add_vectors(
        &mut self,
        family: &str,
        points: &[f64],
        vectors: &[f64],
        weights: &[f64],
        add: impl Fn(&mut FiniteDifferenceInterpolator, [f64; 3], [f64; 3], f64),
    ) -> Result<(), CoreError> {
        let n = count_points(family, "points", points)?;
        check_len(family, "vectors", vectors.len(), 3 * n)?;
        check_len(family, "weights", weights.len(), n)?;
        for i in 0..n {
            check_finite3(family, points, i, "point")?;
            check_vector(family, vectors, i)?;
            check_weight(family, i, weights[i])?;
        }
        let outside = self.count_outside(points);
        // Reserve exactly for the rows that are in the grid.
        let valid = n - outside;
        match family {
            "gradient" => self.interp_mut().reserve_gradient(valid),
            "normal" => self.interp_mut().reserve_norm(valid),
            "tangent" => self.interp_mut().reserve_tangent(valid),
            _ => {}
        }
        for i in 0..n {
            add(
                self.interp_mut(),
                row3(points, i),
                row3(vectors, i),
                weights[i],
            );
        }
        self.warn_outside(family, outside);
        Ok(())
    }

    /// Gradient constraint, direction only.
    pub fn add_gradient(
        &mut self,
        points: &[f64],
        vectors: &[f64],
        weights: &[f64],
    ) -> Result<(), CoreError> {
        self.add_vectors("gradient", points, vectors, weights, |i, p, v, w| {
            i.add_gradient_constraint(p, v, w)
        })
    }

    /// Normal constraint, with length.
    pub fn add_normal(
        &mut self,
        points: &[f64],
        vectors: &[f64],
        weights: &[f64],
    ) -> Result<(), CoreError> {
        self.add_vectors("normal", points, vectors, weights, |i, p, v, w| {
            i.add_norm_constraint(p, v, w)
        })
    }

    pub fn add_tangent(
        &mut self,
        points: &[f64],
        vectors: &[f64],
        weights: &[f64],
    ) -> Result<(), CoreError> {
        self.add_vectors("tangent", points, vectors, weights, |i, p, v, w| {
            i.add_tangent_constraint(p, v, w)
        })
    }

    pub fn add_interface(
        &mut self,
        points: &[f64],
        ids: &[f64],
        weight: f64,
    ) -> Result<(), CoreError> {
        const F: &str = "interface";
        let n = count_points(F, "points", points)?;
        check_len(F, "ids", ids.len(), n)?;
        if !(weight.is_finite() && weight > 0.0) {
            return Err(bad_input(F, "weight must be finite and more than 0"));
        }
        for i in 0..n {
            check_finite3(F, points, i, "point")?;
        }
        let outside = self.count_outside(points);
        let pts: Vec<[f64; 3]> = (0..n).map(|i| row3(points, i)).collect();
        // Upper bound: at most one row for each point.
        self.interp_mut().reserve_interface(n);
        self.interp_mut()
            .add_interface_constraints(&pts, ids, weight);
        self.warn_outside(F, outside);
        Ok(())
    }

    pub fn add_directional_regularisation(
        &mut self,
        directions: &[f64],
        weight: f64,
    ) -> Result<(), CoreError> {
        const F: &str = "directional_regularisation";
        let dof = self.interp().dof();
        check_len(F, "directions", directions.len(), 3 * dof)?;
        if !(weight.is_finite() && weight > 0.0) {
            return Err(bad_input(F, "weight must be finite and more than 0"));
        }
        for i in 0..dof {
            check_finite3(F, directions, i, "direction")?;
        }
        let dirs: Vec<[f64; 3]> = (0..dof).map(|i| row3(directions, i)).collect();
        self.interp_mut()
            .add_directional_regularisation(&dirs, weight);
        self.uses_directional = true;
        Ok(())
    }

    pub fn add_value_bound(
        &mut self,
        points: &[f64],
        lower: &[f64],
        upper: &[f64],
        weights: &[f64],
    ) -> Result<(), CoreError> {
        const F: &str = "value_bound";
        let n = count_points(F, "points", points)?;
        check_len(F, "lower", lower.len(), n)?;
        check_len(F, "upper", upper.len(), n)?;
        check_len(F, "weights", weights.len(), n)?;
        for i in 0..n {
            check_finite3(F, points, i, "point")?;
            check_bounds(F, i, lower[i], upper[i])?;
            check_weight(F, i, weights[i])?;
        }
        let outside = self.count_outside(points);
        self.interp_mut().reserve_value_inequality(n - outside);
        for i in 0..n {
            self.interp_mut().add_value_inequality_constraint(
                row3(points, i),
                lower[i],
                upper[i],
                weights[i],
            );
        }
        self.warn_outside(F, outside);
        Ok(())
    }

    pub fn add_pair_bound(
        &mut self,
        points_a: &[f64],
        points_b: &[f64],
        lower: &[f64],
        upper: &[f64],
        weights: &[f64],
    ) -> Result<(), CoreError> {
        const F: &str = "pair_bound";
        let n = count_points(F, "points_a", points_a)?;
        check_len(F, "points_b", points_b.len(), 3 * n)?;
        check_len(F, "lower", lower.len(), n)?;
        check_len(F, "upper", upper.len(), n)?;
        check_len(F, "weights", weights.len(), n)?;
        for i in 0..n {
            check_finite3(F, points_a, i, "point a")?;
            check_finite3(F, points_b, i, "point b")?;
            check_bounds(F, i, lower[i], upper[i])?;
            check_weight(F, i, weights[i])?;
        }
        // A row is dropped when point a or point b is outside the grid.
        let outside = {
            let grid = &self.interp().grid;
            (0..n)
                .filter(|&i| {
                    grid.trilinear_coeffs(row3(points_a, i)).is_none()
                        || grid.trilinear_coeffs(row3(points_b, i)).is_none()
                })
                .count()
        };
        self.interp_mut().reserve_pair_inequality(n - outside);
        for i in 0..n {
            self.interp_mut().add_inequality_pair_constraint(
                row3(points_a, i),
                row3(points_b, i),
                lower[i],
                upper[i],
                weights[i],
            );
        }
        self.warn_outside(F, outside);
        Ok(())
    }

    pub fn warnings(&self) -> &[String] {
        &self.warnings
    }

    /// Fill the counts of the problem. This allocates nothing.
    fn counts(&self) -> ProblemCounts {
        let i = self.interp();
        ProblemCounts {
            dof: i.dof() as u64,
            n_value: i.n_value_constraints() as u64,
            n_gradient: i.n_gradient_constraints() as u64,
            n_normal: i.n_norm_constraints() as u64,
            n_tangent: i.n_tangent_constraints() as u64,
            n_interface: i.n_interface_constraints() as u64,
            n_pair_bound: i.n_inequality_pair_constraints() as u64,
            n_value_bound: i.n_value_inequality_constraints() as u64,
            uses_directional: self.uses_directional,
        }
    }

    pub fn estimate_memory(&self, s: &Settings) -> Result<MemoryEstimate, CoreError> {
        Ok(estimate(&self.counts(), s))
    }

    pub fn residual_f64(&self, _field: &Field) -> f64 {
        unimplemented!()
    }

    /// Start a solve session. This consumes the problem.
    pub fn start(self, s: Settings) -> Result<SolveSession, CoreError> {
        // Check the budget first. Nothing is allocated here.
        if let Some(budget) = s.budget_bytes {
            let counts = self.counts();
            let e = estimate(&counts, &s);
            if e.total > budget {
                return Err(CoreError::OverBudget {
                    estimate_bytes: e.total,
                    budget_bytes: budget,
                    max_nodes_per_axis: max_grid_for_budget(&counts, &s, budget),
                });
            }
        }
        if self.n_constraints() == 0 {
            return Err(CoreError::NoConstraints);
        }
        if s.precision != Precision::F64 {
            return Err(CoreError::Unsupported {
                setting: "f32 precision".into(),
            });
        }
        SolveSession::new(self, &s)
    }

    /// The interpolator, for the solve session.
    pub(crate) fn interpolator(&self) -> &FiniteDifferenceInterpolator {
        self.interp()
    }

    /// The interpolator, for the solve session (setup).
    pub(crate) fn interpolator_mut(&mut self) -> &mut FiniteDifferenceInterpolator {
        self.interp_mut()
    }

    /// Solve in one call: a loop over `step`.
    pub fn solve(self, s: Settings) -> Result<(Field, SolveReport), CoreError> {
        let mut session = self.start(s)?;
        // One step is large, so the loop has few turns. The CG state stops
        // at `max_iter` by itself.
        while let Step::Running(_) = session.step(1 << 20)? {}
        session.finish()
    }
}
