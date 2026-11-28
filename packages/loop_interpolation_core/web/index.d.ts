// Type declarations for the browser package of loop_interpolation_core.
//
// Two parts:
//   1. The messages of web/worker.js (task T042, plan D6).
//   2. The classes of the wasm package pkg/ (made by wasm-bindgen, see
//      src/wasm.rs for the JavaScript signatures).
//
// Coordinates: every point and vector is in grid coordinates. The solver
// does not rotate or move the data (spec FR-020). The caller does that first.

// ---------------------------------------------------------------------------
// Shared types
// ---------------------------------------------------------------------------

/** Identifies one request. The worker copies it into each reply. */
export type MessageId = number | string;

/** Solver precision. 'f32' stores the field in single precision. */
export type Precision = 'f64' | 'f32';

/** Grid definition. Counts are node counts per axis. */
export interface Grid {
  origin: Float64Array;
  counts: Uint32Array;
  step: Float64Array;
}

/**
 * Settings of estimate and solve. precision, tol and maxIter are always
 * given in a message. The other keys are optional.
 */
export interface Settings {
  precision: Precision;
  tol: number;
  maxIter: number;
  budgetBytes?: number;
  rho?: number;
  innerMaxIter?: number;
  innerTol?: number;
  admmAbsTol?: number;
  admmRelTol?: number;
}

/** Regularisation weights. A missing key means the default weight 1.0. */
export interface Regularisation {
  dxx?: number;
  dyy?: number;
  dzz?: number;
  dxy?: number;
  dxz?: number;
  dyz?: number;
  /** Border weights. */
  dx?: number;
  dy?: number;
  dz?: number;
}

/** Value constraints: the field equals values at points. */
export interface ValueConstraints {
  points: Float64Array;
  values: Float64Array;
  weights: Float64Array;
}

/** Gradient constraints: the gradient of the field equals vectors at points. */
export interface GradientConstraints {
  points: Float64Array;
  vectors: Float64Array;
  weights: Float64Array;
}

/** Normal constraints: the gradient is parallel to vectors at points. */
export interface NormalConstraints {
  points: Float64Array;
  vectors: Float64Array;
  weights: Float64Array;
}

/** Tangent constraints: the gradient is orthogonal to vectors at points. */
export interface TangentConstraints {
  points: Float64Array;
  vectors: Float64Array;
  weights: Float64Array;
}

/** Interface constraints: points with the same id lie on one surface. */
export interface InterfaceConstraints {
  points: Float64Array;
  ids: Float64Array;
  weight: number;
}

/** Inequality bounds on the value at points. */
export interface ValueBoundConstraints {
  points: Float64Array;
  lower: Float64Array;
  upper: Float64Array;
  weights: Float64Array;
}

/** Inequality bounds on the difference of values at two points. */
export interface PairBoundConstraints {
  pointsA: Float64Array;
  pointsB: Float64Array;
  lower: Float64Array;
  upper: Float64Array;
  weights: Float64Array;
}

/** All constraint families. A missing family adds no constraints. */
export interface Constraints {
  value?: ValueConstraints;
  gradient?: GradientConstraints;
  normal?: NormalConstraints;
  tangent?: TangentConstraints;
  interface?: InterfaceConstraints;
  valueBound?: ValueBoundConstraints;
  pairBound?: PairBoundConstraints;
}

/** Memory estimate in bytes, as the worker posts it (a plain object). */
export interface MemoryEstimate {
  total: number;
  vectors: number;
  constraints: number;
  working: number;
}

// ---------------------------------------------------------------------------
// Messages in (page -> worker)
// ---------------------------------------------------------------------------

/** Starts the module load. The worker sends no reply. */
export interface InitMessage {
  type: 'init';
}

/** Asks for the memory estimate. The reply is a ResultMessage with estimate. */
export interface EstimateMessage {
  type: 'estimate';
  id: MessageId;
  grid: Grid;
  constraints?: Constraints;
  settings: Settings;
}

/** Starts a solve. The reply is progress messages, then a result or an error. */
export interface SolveMessage {
  type: 'solve';
  id: MessageId;
  grid: Grid;
  regularisation?: Regularisation;
  constraints?: Constraints;
  settings: Settings;
}

/**
 * Cancels a solve. Without id, the running solve is cancelled. The worker
 * answers with 'cancelled'. If no solve runs, it answers with 'cancelled' too.
 */
export interface CancelMessage {
  type: 'cancel';
  id?: MessageId;
}

export type WorkerInMessage = InitMessage | EstimateMessage | SolveMessage | CancelMessage;

// ---------------------------------------------------------------------------
// Messages out (worker -> page)
// ---------------------------------------------------------------------------

/** Solve progress. Sent after the first step and then at most every 100 ms. */
export interface ProgressMessage {
  type: 'progress';
  id: MessageId;
  stage: 'setup' | 'solve' | 'check';
  iteration: number;
  residual: number;
  maxViolation: number;
}

/** Result of a finished solve. The field buffer is transferred, not copied. */
export interface SolveResult {
  type: 'result';
  id: MessageId;
  field: Float64Array | Float32Array;
  memoryBytes: number;
  done: true;
  converged: boolean;
  peakBytes: number;
  iterations: number;
  residual: number;
  maxViolation: number;
  warnings: string[];
}

/** Result of an estimate message. */
export interface EstimateResult {
  type: 'result';
  id: MessageId;
  estimate: MemoryEstimate;
}

export type ResultMessage = SolveResult | EstimateResult;

/** Confirms a cancel. */
export interface CancelledMessage {
  type: 'cancelled';
  id?: MessageId;
}

/** Error body. JSON numbers that are not finite arrive as null. */
export type ErrorBody =
  | { kind: 'OverBudget'; estimate: MemoryEstimate; budget: number; maxNodesPerAxis: number }
  | { kind: 'SolveRunning' }
  | { kind: 'InvalidInput'; constraint: string; reason: string }
  | { kind: 'GridTooSmall'; axis: number; nodes: number }
  | { kind: 'GridTooLarge'; nodes: number; maxNodes: number }
  | { kind: 'NoConstraints' }
  | { kind: 'OutOfMemory'; requestedBytes: number }
  | {
      kind: 'Infeasible';
      worst: Array<{ index: number; kind: 'value' | 'pair'; violation: number | null }>;
    }
  | { kind: 'PointOutsideGrid'; index: number }
  | { kind: 'Unsupported'; setting: string }
  | { kind: 'Panic'; message: string }
  | { kind: 'WorkerLoadError'; message: string }
  | { kind: 'InternalError'; message: string };

/** Error kinds. */
export type ErrorKind = ErrorBody['kind'];

/** Error reply. id is missing when the error is not tied to a request. */
export type ErrorMessage = { type: 'error'; id?: MessageId } & ErrorBody;

export type WorkerOutMessage =
  | ProgressMessage
  | ResultMessage
  | CancelledMessage
  | ErrorMessage;

// ---------------------------------------------------------------------------
// Classes of the wasm package (pkg/loop_interpolation_core.js)
// ---------------------------------------------------------------------------

/** Memory estimate with a free() method (the wasm object). */
export declare class Estimate {
  readonly total: number;
  readonly vectors: number;
  readonly constraints: number;
  readonly working: number;
  free(): void;
}

/** Result of one Session.step call. converged, peakBytes and warnings are valid when done is true. */
export declare class StepResult {
  readonly done: boolean;
  readonly stage: 'setup' | 'solve' | 'check' | 'done';
  readonly iteration: number;
  readonly residual: number;
  readonly violation: number;
  readonly converged: boolean;
  readonly peakBytes: number;
  readonly warnings: string[];
  free(): void;
}

/** An interpolation problem. */
export declare class Problem {
  constructor(origin: Float64Array, counts: Uint32Array, step: Float64Array);
  setRegularisation(weights: Regularisation): void;
  addValue(points: Float64Array, values: Float64Array, weights: Float64Array): void;
  addGradient(points: Float64Array, vectors: Float64Array, weights: Float64Array): void;
  addNormal(points: Float64Array, vectors: Float64Array, weights: Float64Array): void;
  addTangent(points: Float64Array, vectors: Float64Array, weights: Float64Array): void;
  addInterface(points: Float64Array, ids: Float64Array, weight: number): void;
  addValueBound(
    points: Float64Array,
    lower: Float64Array,
    upper: Float64Array,
    weights: Float64Array,
  ): void;
  addPairBound(
    pointsA: Float64Array,
    pointsB: Float64Array,
    lower: Float64Array,
    upper: Float64Array,
    weights: Float64Array,
  ): void;
  /** Memory estimate. Throws a JsError whose message is JSON with a kind key. */
  estimate(settings?: Partial<Settings>): Estimate;
  /** Starts a solve. This consumes the problem. Do not use it again. */
  start(settings?: Partial<Settings>): Session;
  free(): void;
}

/** A running solve. */
export declare class Session {
  /** Runs at most maxIters iterations. After done, it keeps returning done. */
  step(maxIters: number): StepResult;
  /** Copies the field out. Before the end, it is the best result so far. */
  field(): Float64Array | Float32Array;
  evaluate(points: Float64Array): Float64Array;
  evaluateGradient(points: Float64Array): Float64Array;
  /** Drops the session. */
  free(): void;
}

/** Size of the wasm linear memory in bytes. The memory never shrinks. */
export declare function memoryBytes(): number;
