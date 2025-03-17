#include "scurve.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>

#define EPSILON (1e-6f)
#define MAX_DIST_ERROR (6.0)
#define NEWTON_ITERATIONS (10)
#define NEWTON_TOLERANCE (1e-3f)
#define NEWTON_RELAX_FACTOR (1.5f)
#define VP_GUESS_ITERATIONS (5)

#define ON_END_RESET_VAJ

typedef struct
{
    float sign;   ///< Sign of the motion (1 or -1).
    float dist;   ///< Total distance to travel.
    float s1_end; // (seconds) cumulative end time of segment 1
    float s2_end; // segment 2 end, etc.
    float s3_end;
    float s4_end;
    float s5_end;
    float s6_end;
    float totalTime; // end of segment 7
} scurve_param_t;

typedef struct scurve_segment_state
{
    float x0;                          // position at segment start
    float v0;                          // velocity at segment start
    float a0;                          // acceleration at segment start
    float j0;                          // jerk at segment start
    float t0;                          // time offset (when the segment started)
    float tEnd;                        // time when the segment ends
    int segmentIndex;                  // 1..7 or idle
    struct scurve_segment_state *next; // pointer to next segment state
} scurve_segment_state_t;

struct scurve_t_
{
    scurve_state_e state; ///< Current motion profile state (IDLE, INIT, CONF, BUSY, QUICKSTOP).
    float Ts;             ///< Sample time (seconds).
    float pos_target;     ///< Target setpoint.
    float vel_start;      ///< Starting velocity.
    float vel_stop;       ///< Stopping velocity.
    float vel_max;        ///< Maximum allowed velocity in the profile.
    float acc_max;        ///< Maximum allowed acceleration in the profile.
    float dec_max;        ///< Maximum allowed deceleration in the profile.
    float jerk_max;       ///< Maximum allowed jerk in the profile.
    float vel_peak;       ///< Peak velocity in the profile.
    float vel_start_clip; ///< Start velocity after clipping.
    float pos_out;        ///< Profiled position output.
    float vel_out;        ///< Profiled velocity output.
    float acc_out;        ///< Profiled acceleration output.
    float jerk_out;       ///< Profiled jerk output.

    uint32_t profile_ticks; ///< Current tick counter for profile progression.

    // keep track of current segment
    int segStateIndex;
    scurve_segment_state_t segState[7];

    scurve_param_t param; ///< Precomputed parameters and thresholds for the motion profile.
};

/* ------------------------------------------------------------------
   Helper: Computes distance and time for a velocity transition from V1 to V2
   with jerk-limited acceleration/deceleration.

   Two cases are handled:
   1. Full trapezoidal profile: When Δv > A²/J
      - Three phases: jerk-up, constant acceleration, jerk-down
      - Jerk phases: Tj = A/J each
      - Constant acceleration phase: Tm = (Δv - A·Tj)/A
      - Total time: Ta = 2·Tj + Tm
      
   2. Triangular profile: When Δv ≤ A²/J
      - Two phases only: jerk-up and jerk-down
      - No constant acceleration phase
      - Total time: Tshort = 2·sqrt(Δv/J)
      
   Returns:
   - distance: Total distance covered during the transition
   - time: Total duration of the transition
   - Tj: Duration of jerk phase
   - Ta: Total acceleration phase duration
   ------------------------------------------------------------------ */
typedef struct
{
    float distance; // total distance covered
    float time;     // total time (T1 + T2 + T3)
    float Tj;       // jerk-up or jerk-down time
    float Ta;       // entire accelerate phase (Tj + Tm + Tj)
} seg_scurve_t;

// Integrates motion equations for jerk-limited segments (1,3,5,7)
// Uses standard jerk-limited polynomials to compute position, velocity and acceleration
// at time dt from initial conditions x0, v0, a0 with constant jerk j
static inline void __attribute__((optimize("O3"))) integrate_segment1357(float x0, float v0, float a0, float j, float dt, float *xEnd, float *vEnd, float *aEnd)
{
    // Standard jerk-limited motion equations:
    float a_local = a0 + j * dt;                                              // a(t) = a0 + jt
    float v_local = v0 + a0 * dt + 0.5f * j * (dt * dt);                     // v(t) = v0 + a0·t + ½jt²
    float x_local = x0 + v0 * dt + 0.5f * a0 * (dt * dt) + (j * (dt * dt * dt)) / 6.0f; // x(t) = x0 + v0·t + ½a0·t² + ⅙jt³

    *xEnd = x_local;
    *vEnd = v_local;
    *aEnd = a_local;
}

// Integrates motion equations for constant acceleration segments (2,6)
// Uses standard constant acceleration equations to compute position, velocity
// and acceleration at time dt from initial conditions x0, v0, a0
static inline void __attribute__((optimize("O3"))) integrate_segment26(float x0, float v0, float a0, float dt, float *xEnd, float *vEnd, float *aEnd)
{
    // Standard constant acceleration equations:
    float a_local = a0;                                    // a(t) = a0
    float v_local = v0 + a0 * dt;                         // v(t) = v0 + a0·t
    float x_local = x0 + v0 * dt + 0.5f * a0 * (dt * dt); // x(t) = x0 + v0·t + ½a0·t²

    *xEnd = x_local;
    *vEnd = v_local;
    *aEnd = a_local;
}

// Integrates motion equations for constant velocity segment (4)
// Uses simple constant velocity equations to compute position and velocity
// at time dt from initial conditions x0, v0
static inline void __attribute__((optimize("O3"))) integrate_segment4(float x0, float v0, float dt, float *xEnd, float *vEnd)
{
    // Constant velocity equations:
    float v_local = v0;           // v(t) = v0
    float x_local = x0 + v0 * dt; // x(t) = x0 + v0·t

    *xEnd = x_local;
    *vEnd = v_local;
}

// Compute s-curve distance/time from Vs -> Ve
// with max accel = A, max jerk = J.
// This formula is well-known; it breaks down into two cases:
static void compute_scurve_distance(float Vs, float Ve, float A, float J, seg_scurve_t *out)
{
    float dv = (Ve - Vs); // We assume Ve > Vs in typical usage

    if (dv < EPSILON)
    {
        // No real velocity change
        out->distance = 0.0f;
        out->time = 0.0f;
        out->Tj = 0.0f;
        out->Ta = 0.0f;
        return;
    }

    float dv_threshold = (A * A) / J; // threshold for full s-curve (trapezoidal)
    if (dv > dv_threshold + EPSILON)
    {
        // Full (trapezoidal) profile: maximum acceleration A is reached
        float Tj = A / J;                // time for jerk phase
        float Tflat = (dv - A * Tj) / A; // duration of constant acceleration phase
        float Ta = 2.0f * Tj + Tflat;    // total time for the acceleration phase
        // Distance computed via trapezoidal integration of velocity
        float dist_acc = (Vs * Ta) + 0.5f * dv * Ta;

        out->distance = dist_acc;
        out->time = Ta;
        out->Tj = Tj;
        out->Ta = Ta;
    }
    else
    {
        // Triangular profile: not enough dv to reach full acceleration A
        float Tshort = 2.0f * sqrtf(dv / J);
        float dist_acc = (Vs * Tshort) + (dv * Tshort / 2.0f);

        out->distance = dist_acc;
        out->time = Tshort;
        out->Tj = 0.5f * Tshort; // jerk phase duration
        out->Ta = Tshort;
    }
}

/*
 * Computes the optimal peak velocity for an asymmetric 7-segment S-curve profile.
 * 
 * This function finds the peak velocity Vp that satisfies the total distance constraint D
 * while respecting acceleration and deceleration limits. The solution handles both
 * trapezoidal and triangular acceleration/deceleration phases.
 * 
 * Algorithm:
 * 1. First attempts a closed-form solution assuming full trapezoidal phases
 * 2. If the trapezoidal solution is invalid, uses Newton's method to refine
 *    the solution considering the actual acceleration/deceleration regimes
 * 3. Handles special cases where Vp must be clipped to Vs or Ve
 * 
 * Parameters:
 * - Vs: Start velocity
 * - Vm: Maximum allowed velocity
 * - Ve: End velocity
 * - A: Maximum acceleration
 * - Dec: Maximum deceleration
 * - J: Maximum jerk
 * - D: Total distance to travel
 * 
 * Returns:
 * The optimal peak velocity that satisfies the distance constraint
 */
static float scurve_findVpeakClosedForm(float Vs, float Vm, float Ve, float A, float Dec, float J, float D)
{
    // Thresholds for full trapezoidal phases
    float Vp_trap_acc = Vs + (A * A) / J;
    float Vp_trap_dec = Ve + (Dec * Dec) / J;
    // Vp must be at least the max of Vs and Ve
    float Vp_min = fmaxf(Vs, Ve);

    // ---- 1) Compute the full trapezoidal candidate via quadratic equation ----
    // Using: C2 * Vp^2 + C1 * Vp + C0 = 0, where
    //   C2 = 0.5*(1/A + 1/Dec)
    //   C1 = (A + Dec)/(2*J)
    //   C0 = -0.5*(Vs^2/A + Ve^2/Dec) + (A*Vs+Dec*Ve)/(2*J) - D
    float C2 = 0.5f * (1.0f / A + 1.0f / Dec);
    float C1 = (A + Dec) / (2.0f * J);
    float C0 = -0.5f * (Vs * Vs / A + Ve * Ve / Dec) + (A * Vs + Dec * Ve) / (2.0f * J) - D;
    float disc = C1 * C1 - 4.0f * C2 * C0;
    if (disc < 0.0f)
    {
        disc = 0.0f;
    }
    float Vp_candidate = (-C1 + sqrtf(disc)) / (2.0f * C2);

    // If both phases can be full trapezoidal, return candidate.
    if ((Vp_candidate >= Vp_trap_acc) && (Vp_candidate >= Vp_trap_dec))
    {
        return Vp_candidate;
    }

    // Determine regimes
    bool accTri = (Vp_candidate < Vp_trap_acc);
    bool decTri = (Vp_candidate < Vp_trap_dec);

    // Otherwise, refine Vp using Newton's method with the appropriate distance functions.
    float Vp = fmaxf(Vp_min, Vp_candidate);
    if ((Vp - Vp_min) < EPSILON)
    {
        Vp = Vp_min + 1e-3f;
    }
    for (int i = 0; i < VP_GUESS_ITERATIONS; i++)
    {
        // Compute regime flags based on current Vp.
        bool localAccTri = (Vp < Vp_trap_acc);
        bool localDecTri = (Vp < Vp_trap_dec);

        // Compute f(Vp) = f_acc + f_dec - D and its derivative df.
        float f_val = 0.0f;
        float df = 0.0f;

        // Acceleration phase function
        if (localAccTri)
        {
            // Triangular acceleration phase:
            // f_acc = (Vp + Vs) * sqrt((Vp - Vs)/J)
            float diff = Vp - Vs;
            float sqrt_term = sqrtf(diff / J);
            f_val += (Vp + Vs) * sqrt_term;
            df += sqrt_term + (Vp + Vs) / (2.0f * sqrtf(J * diff));
        }
        else
        {
            // Trapezoidal acceleration phase:
            // f_acc = 0.5 * (Vp + Vs) * [ (Vp - Vs)/A + (A/J) ]
            float diff = Vp - Vs;
            float term = diff / A + A / J;
            f_val += 0.5f * (Vp + Vs) * term;
            df += Vp / A + 0.5f * (A / J);
        }

        // Deceleration phase function
        if (localDecTri)
        {
            // Triangular deceleration phase:
            // f_dec = (Vp + Ve) * sqrt((Vp - Ve)/J)
            float diff = Vp - Ve;
            float sqrt_term = sqrtf(diff / J);
            f_val += (Vp + Ve) * sqrt_term;
            df += sqrt_term + (Vp + Ve) / (2.0f * sqrtf(J * diff));
        }
        else
        {
            // Trapezoidal deceleration phase:
            // f_dec = 0.5 * (Vp + Ve) * [ (Vp - Ve)/Dec + (Dec/J) ]
            float diff = Vp - Ve;
            float term = diff / Dec + Dec / J;
            f_val += 0.5f * (Vp + Ve) * term;
            df += Vp / Dec + 0.5f * (Dec / J);
        }

        // Total function value F(Vp)= f_acc + f_dec - D
        f_val -= D;

        // Newton update:
        float newVp = Vp - (f_val / df);
        // Instead of always clamping to Vp_min+epsilon,
        // if newVp < Vp_min, take the average so the descent is more gradual.
        if (newVp < Vp_min)
        {
            newVp = 0.5f * (Vp + Vp_min);
        }
        Vp = newVp;
    }
    return Vp;
}

/* 
 * Calculates parameters for a 7-segment S-curve motion profile
 * 
 * The profile consists of:
 * 1. Acceleration phase (segments 1-3):
 *    - S1: Jerk-up to max acceleration
 *    - S2: Constant acceleration
 *    - S3: Jerk-down to zero acceleration
 * 2. Constant velocity phase (segment 4)
 * 3. Deceleration phase (segments 5-7):
 *    - S5: Jerk-up to max deceleration
 *    - S6: Constant deceleration
 *    - S7: Jerk-down to zero deceleration
 * 
 * The function handles special cases:
 * - When start velocity > end velocity: Only deceleration phase
 * - When end velocity > start velocity: Only acceleration phase
 * - When distance is short: No constant velocity phase
 * - When velocity changes are small: Triangular acceleration/deceleration profiles
 */
static void calc_param(scurve_t *obj)
{
    float Ts = obj->Ts;
    float dist = obj->pos_target - obj->pos_out;
    obj->param.dist = dist;
    obj->param.sign = (dist >= 0.0f) ? +1.0f : -1.0f;
    float D = fabsf(dist);     // absolute distance to be covered
    float Vs = obj->vel_start; // start velocity
    float Ve = obj->vel_stop;  // end velocity
    float Vm = obj->vel_max;   // allowable max velocity
    float A = obj->acc_max;    // allowable max acceleration
    float Dec = obj->dec_max;  // allowable max deceleration
    float J = obj->jerk_max;   // allowable max jerk

    /* ------------------------------------------------------------------
       Compute minimal distance/time to:
          - accelerate from Vs -> Vm with (A, J)
          - decelerate from Vm -> Ve with (Dec, J)
       ------------------------------------------------------------------ */
    seg_scurve_t accPart = {0};
    seg_scurve_t decPart = {0};

    // accelerate from Vs to Vm
    compute_scurve_distance(Vs, Vm, A, J, &accPart);
    // decelerate from Vm to Ve
    compute_scurve_distance(Ve, Vm, Dec, J, &decPart);

    float s_acc = accPart.distance;
    float s_dec = decPart.distance;
    float t_acc = accPart.time;
    float t_dec = decPart.time;
    float s_sum = s_acc + s_dec;

    /* ------------------------------------------------------------------
       4) Check if we have room for a cruise segment
       ------------------------------------------------------------------ */
    float s_cruise = 0.0f;
    float t_cruise = 0.0f;
    if (s_sum < D)
    {
        s_cruise = D - s_sum;
        // time cruising (constant velocity = Vm):
        t_cruise = s_cruise / Vm;
        obj->vel_peak = Vm;
        obj->vel_start_clip = Vs;
    }
    else
    {
        bool vstart_clipped = false;
        if (Vs > Ve)
        {
            // We assume s_acc = 0, so we only have a deceleration phase from Vs down to Ve.
            // If that distance alone exceeds the target D, we must reduce Vs.
            compute_scurve_distance(Ve, Vs, Dec, J, &decPart);
            if (decPart.distance > D)
            {
                float dv = Vs - Ve;
                float dv_threshold = (Dec * Dec) / J;
                float VsCandidate;
                if (dv <= dv_threshold + EPSILON)
                {
                    // Triangular deceleration profile (no constant deceleration phase)
                    // -------------------------------------------------------
                    // Original formulation was numerically sensitive.
                    // Here we solve:
                    //    (x + 2*Ve) * sqrt(x/J) = D, with x = VsCandidate - Ve.
                    // Let w = sqrt(x). Then:
                    //    w^3 + 2*Ve*w = sqrt(J)*D.
                    //
                    // Solve for w:
                    float half_term = (sqrtf(J) * D) / 2.0f;
                    float delta_val = half_term * half_term + (8.0f * Ve * Ve * Ve) / 27.0f;
                    float w = cbrtf(half_term + sqrtf(delta_val)) +
                              cbrtf(half_term - sqrtf(delta_val));
                    // -- NEWTON REFINEMENT --
                    for (int i = 0; i < NEWTON_ITERATIONS; i++)
                    {
                        float f_val = w * w * w + 2.0f * Ve * w - sqrtf(J) * D;
                        if (fabsf(f_val) < NEWTON_TOLERANCE)
                            break; // Converged sufficiently
                        float df = 3.0f * w * w + 2.0f * Ve;
                        float dw = f_val / df;
                        w = w - NEWTON_RELAX_FACTOR * dw;
                    }
                    // Recover x and then VsCandidate:
                    float x = w * w;
                    VsCandidate = Ve + x;
                }
                else
                {
                    // Trapezoidal deceleration profile is possible (with a constant deceleration phase)
                    // -------------------------------------------------------
                    // We want newVs such that:
                    //     ((Vs + Ve - δ) / 2) * ( ((Vs - Ve - δ) / Dec) + (Dec / J) ) = D
                    // with δ = Vs - newVs.
                    // Define:
                    float C = Vs + Ve;
                    float Vdelta = Vs - Ve;
                    float E = Vdelta + (Dec * Dec) / J;
                    float disc_val = (C + E) * (C + E) - 4.0f * (C * E - 2.0f * Dec * D);
                    if (disc_val < 0.0f)
                        disc_val = 0.0f;
                    // Choose the smaller positive root:
                    float delta_trap = ((C + E) - sqrtf(disc_val)) / 2.0f;
                    VsCandidate = Vs - delta_trap;

                    if ((VsCandidate - Ve) < dv_threshold)
                    {
                        // If the computed reduction is too small, recalculate using the triangular approach.
                        float half_term = (sqrtf(J) * D) / 2.0f;
                        float delta_val = half_term * half_term + (8.0f * Ve * Ve * Ve) / 27.0f;
                        float w = cbrtf(half_term + sqrtf(delta_val)) +
                                  cbrtf(half_term - sqrtf(delta_val));
                        // Refine w using Newton iteration for better accuracy:
                        for (int i = 0; i < NEWTON_ITERATIONS; i++)
                        {
                            float f_val = w * w * w + 2.0f * Ve * w - sqrtf(J) * D;
                            if (fabsf(f_val) < NEWTON_TOLERANCE)
                                break; // Converged sufficiently
                            float df = 3.0f * w * w + 2.0f * Ve;
                            float dw = f_val / df;
                            w = w - NEWTON_RELAX_FACTOR * dw;
                        }
                        float x = w * w;
                        VsCandidate = Ve + x;
                    }
                }

                // Clamp if necessary.
                if (VsCandidate < Ve)
                    VsCandidate = Ve;
                else if (VsCandidate > Vs)
                    VsCandidate = Vs;

                // -------------------------------------------------------
                // 2) Check if this candidate actually produces a deceleration distance ~ D.
                //    If so, mark that we clipped Vs.
                // -------------------------------------------------------
                compute_scurve_distance(Ve, VsCandidate, Dec, J, &decPart);
                float distError = decPart.distance - D;
                if (fabsf(distError) <= MAX_DIST_ERROR)
                {
                    // Mark that we clipped Vs
                    accPart.distance = 0;
                    accPart.time = 0;
                    accPart.Tj = 0;
                    vstart_clipped = true;

                    // Update Vs to the clipped value.
                    Vs = VsCandidate;
                }
            }
        }
        else if (Ve > Vs)
        {
            // We assume s_dec = 0, so we only have an acceleration phase from Vs up to Ve.
            // If that distance alone exceeds the target D, we must reduce Ve.
            compute_scurve_distance(Vs, Ve, A, J, &accPart);
            if (accPart.distance > D)
            {
                float dv = Ve - Vs;
                float dv_threshold = A * A / J;
                float VeCandidate, delta;

                if (dv <= dv_threshold + EPSILON)
                {
                    // Triangular acceleration profile (no constant acceleration phase)
                    // -------------------------------------------------------
                    // 1) Cubic formula approach.
                    // Solve: delta^3 + 4*Vs*delta^2 + 4*Vs^2*delta - J*D^2 = 0,
                    // with delta = newVe - Vs.
                    // -------------------------------------------------------
                    float a_coef = 4.0f * Vs; // helps form "4*Vs*delta^2"
                    float p_coef = -(4.0f * Vs * Vs) / 3.0f;
                    float q_coef = -(16.0f * Vs * Vs * Vs) / 27.0f - (J * D * D);
                    float half_q = q_coef / 2.0f;
                    float disc = half_q * half_q + powf(p_coef / 3.0f, 3);
                    float tol = 1e-6f * fabsf(half_q * half_q);
                    if (disc < 0.0f && fabsf(disc) < tol)
                    {
                        // tiny negative discriminant -> treat as zero
                        disc = 0.0f;
                    }
                    float sqrt_val = sqrtf(disc);

                    // Principal real root
                    float y = cbrtf(-half_q + sqrt_val) + cbrtf(-half_q - sqrt_val);
                    delta = y - (a_coef / 3.0f);
                    VeCandidate = Vs + delta;

                    // If that candidate is out of [Vs, Ve], try the alternative sign.
                    if (VeCandidate < Vs || VeCandidate > Ve)
                    {
                        float y_alt = cbrtf(-half_q - sqrt_val) + cbrtf(-half_q + sqrt_val);
                        float delta_alt = y_alt - (a_coef / 3.0f);
                        float VeCandidateAlt = Vs + delta_alt;
                        // Use whichever candidate is in [Vs, Ve].
                        if (VeCandidateAlt >= Vs && VeCandidateAlt <= Ve)
                            VeCandidate = VeCandidateAlt;
                    }
                }
                else
                {
                    // Trapezoidal acceleration profile is possible (with a constant acceleration phase)
                    // -------------------------------------------------------
                    // Here we solve for delta = newVe - Vs using the relation:
                    //    (Vs + 0.5*delta) * (delta/A + A/J) = D.
                    //
                    // Expanding and rearranging leads to a quadratic:
                    //    delta^2 + (2*Vs + A^2/J)*delta + (2*Vs*A^2/J - 2*A*D) = 0.
                    // -------------------------------------------------------
                    float b = 2.0f * Vs + (A * A) / J;
                    float c_term = (2.0f * Vs * A * A) / J - 2.0f * A * D;
                    float disc = b * b - 4.0f * c_term;
                    if (disc < 0.0f)
                        disc = 0.0f;
                    float sqrt_disc = sqrtf(disc);
                    float delta1 = (-b + sqrt_disc) / 2.0f;
                    float delta2 = (-b - sqrt_disc) / 2.0f;
                    // Select a valid positive root that is no smaller than the trapezoidal limit.
                    if (delta1 >= 0.0f && delta1 >= dv_threshold && (Vs + delta1) <= Ve)
                        delta = delta1;
                    else if (delta2 >= 0.0f && delta2 >= dv_threshold && (Vs + delta2) <= Ve)
                        delta = delta2;
                    else
                        delta = fmaxf(delta1, delta2); // fallback: choose the larger root

                    VeCandidate = Vs + delta;

                    if (delta < dv_threshold)
                    {
                        // The computed δ is too small; switch to the triangular (cubic) solution
                        float a_coef = 4.0f * Vs;
                        float p_coef = -(4.0f * Vs * Vs) / 3.0f;
                        float q_coef = -(16.0f * Vs * Vs * Vs) / 27.0f - (J * D * D);
                        float half_q = q_coef / 2.0f;
                        float disc = half_q * half_q + powf(p_coef / 3.0f, 3);
                        float tol = 1e-6f * fabsf(half_q * half_q);
                        if (disc < 0.0f && fabsf(disc) < tol)
                        {
                            disc = 0.0f;
                        }
                        float sqrt_val = sqrtf(disc);
                        float y = cbrtf(-half_q + sqrt_val) + cbrtf(-half_q - sqrt_val);
                        delta = y - (a_coef / 3.0f);
                        VeCandidate = Vs + delta;
                    }
                }

                // Clamp if necessary.
                if (VeCandidate < Vs)
                    VeCandidate = Vs;
                else if (VeCandidate > Ve)
                    VeCandidate = Ve;

                // -------------------------------------------------------
                // 2) Check if this candidate actually produces an acceleration distance ~ D.
                //    If still overshooting, a bracket/bisection fallback may be used.
                // -------------------------------------------------------
                compute_scurve_distance(Vs, VeCandidate, A, J, &accPart);
                float distError = accPart.distance - D;
                if (fabsf(distError) <= MAX_DIST_ERROR)
                {
                    // Mark that we clipped Ve
                    decPart.distance = 0;
                    decPart.time = 0;
                    decPart.Tj = 0;
                    vstart_clipped = true;
                }
            }
        }
        if (!vstart_clipped)
        {
            float Vp = scurve_findVpeakClosedForm(Vs, Vm, Ve, A, Dec, J, D);
            if (Vp <= Ve)
            {
                Ve = Vp;
                Vm = Vp;
                compute_scurve_distance(Vs, Vm, A, J, &accPart);
                decPart.distance = 0.0f;
                decPart.time = 0.0f;
                decPart.Tj = 0.0f;
            }
            else if (Vp < Vs)
            {
                Vs = Vp;
                Vm = Vp;
                compute_scurve_distance(Vs, Vm, Dec, J, &decPart);
                accPart.distance = 0.0f;
                accPart.time = 0.0f;
                accPart.Tj = 0.0f;
            }
            else
            {
                Vm = Vp;
                compute_scurve_distance(Vs, Vm, A, J, &accPart);
                compute_scurve_distance(Ve, Vm, Dec, J, &decPart);
            }
        }
        obj->vel_peak = Vm;
        obj->vel_start_clip = Vs;
        s_acc = accPart.distance;
        s_dec = decPart.distance;
        t_acc = accPart.time;
        t_dec = decPart.time;
        if (s_acc + s_dec < D)
        {
            s_cruise = D - (s_acc + s_dec);
            t_cruise = s_cruise / Vm;
        }
        else
        {
            s_cruise = 0.0f;
            t_cruise = 0.0f;
        }
    }

    /* ------------------------------------------------------------------
       5) We now have total times for 7 "segments":
          Seg1..3: accelerate from Vs..Vm   => total time t_acc
          Seg4:    constant velocity        => time t_cruise
          Seg5..7: decelerate from Vm..Ve   => total time t_dec
       ------------------------------------------------------------------ */
    float T1 = accPart.Tj;                // jerk-up time
    float T2 = t_acc - 2.0f * accPart.Tj; // "flat accel" time
    float T3 = accPart.Tj;                // jerk-down time
    float T4 = t_cruise;                  // constant velocity
    float T5 = decPart.Tj;                // dec jerk up
    float T6 = t_dec - 2.0f * decPart.Tj; // dec flat
    float T7 = decPart.Tj;                // dec jerk down

    // Store in your param struct:
    obj->param.s1_end = T1;
    obj->param.s2_end = obj->param.s1_end + T2;
    obj->param.s3_end = obj->param.s2_end + T3;
    obj->param.s4_end = obj->param.s3_end + T4;
    obj->param.s5_end = obj->param.s4_end + T5;
    obj->param.s6_end = obj->param.s5_end + T6;
    obj->param.totalTime = obj->param.s6_end + T7;
}

// Evaluates the motion state (position, velocity, acceleration, jerk)
// at the given time t_sec within a segment
// The segment type determines which integration function is used:
// - Segments 1,3,5,7: Jerk-limited motion
// - Segments 2,6: Constant acceleration
// - Segment 4: Constant velocity
static void segment_eval(scurve_segment_state_t *seg, float t_sec, float *pos, float *vel, float *acc, float *jrk)
{
    float dt = t_sec - seg->t0; // local time within the segment

    // Select appropriate integration based on segment type
    switch (seg->segmentIndex)
    {
    case 1:
    case 3:
    case 5:
    case 7:
        // Jerk-limited segments
        integrate_segment1357(seg->x0, seg->v0, seg->a0,
                              seg->j0, dt, pos, vel, acc);
        *jrk = seg->j0;
        break;
    case 2:
    case 6:
        // Constant acceleration segments
        integrate_segment26(seg->x0, seg->v0, seg->a0,
                            dt, pos, vel, acc);
        *jrk = 0.0f;
        break;
    case 4:
        // Constant velocity segment
        integrate_segment4(seg->x0, seg->v0,
                           dt, pos, vel);
        *jrk = 0.0f;
        *acc = 0.0f;
        break;
    }
    return;
}

/**
 * Initializes an S-curve motion profile generator
 * 
 * @param buffer Pointer to memory buffer for the scurve object
 * @param buffer_size Size of the provided buffer (must be >= sizeof(scurve_t))
 * @param sampleTime_s Sample time in seconds for profile generation
 * @return Pointer to initialized scurve object, or NULL if buffer too small
 */
scurve_t *scurve_init(void *buffer, size_t buffer_size, float sampleTime_s)
{
    if (buffer_size < sizeof(scurve_t))
    {
        return NULL; // Buffer too small
    }

    scurve_t *obj = (scurve_t *)buffer;

    obj->state = CURVE_IDLE;
    obj->Ts = sampleTime_s;
    obj->pos_target = 0.0f;
    obj->pos_out = 0.0f;
    obj->vel_out = 0.0f;
    obj->acc_out = 0.0f;
    obj->jerk_out = 0.0f;
    obj->profile_ticks = 0;
    obj->vel_start = 0.0f;
    obj->vel_stop = 0.0f;
    obj->vel_max = 0.0f;
    obj->acc_max = 0.0f;
    obj->jerk_max = 0.0f;
    return obj;
}

// State machine implementation for S-curve profile generation
// States:
// - IDLE: Waiting for new profile
// - CONF: Computing profile parameters
// - BUSY: Executing profile segments
// - ONEND: Profile complete, resetting outputs
static scurve_state_e scurve_run(scurve_t *obj)
{
    switch (obj->state)
    {
    case CURVE_IDLE:
        // do nothing
        obj->profile_ticks = 0;
        break;

    case CURVE_CONF:
    {
        calc_param(obj);

        float segEnd[7] = {
            obj->param.s1_end,
            obj->param.s2_end,
            obj->param.s3_end,
            obj->param.s4_end,
            obj->param.s5_end,
            obj->param.s6_end,
            obj->param.totalTime};

        // Create the first segment always:
        int segCount = 0;      // count of created segments
        float tCurrent = 0.0f; // current start time marker
        scurve_segment_state_t *firstSeg = &obj->segState[segCount];
        segCount++; // first segment (segment 1)

        // Initialize the segment 1
        firstSeg->segmentIndex = 1;
        firstSeg->t0 = 0.0f;
        firstSeg->tEnd = segEnd[0];
        firstSeg->x0 = obj->pos_out;
        firstSeg->v0 = obj->param.sign * obj->vel_start_clip; // vel_start <= vel_stop <= vel_peak <= vel_max
        firstSeg->a0 = 0.0f;
        firstSeg->j0 = obj->param.sign * obj->jerk_max;

        // Keep track of the last created segment:
        scurve_segment_state_t *prevCreatedSeg = firstSeg;
        tCurrent = segEnd[0];

        // For each of the remaining 6 segments, use segEnd[] values.
        for (int i = 1; i < 7; i++)
        {
            float tEnd = segEnd[i];
            float duration = tEnd - tCurrent;

            if (duration > EPSILON)
            {
                scurve_segment_state_t *seg = &obj->segState[segCount];
                segCount++;

                seg->segmentIndex = i + 1; // segments numbered 2..7
                seg->t0 = tCurrent;
                seg->tEnd = tEnd;
                seg->next = NULL;

                // Evaluate position/velocity/acceleration from end of the previous segment.
                float prev_xEnd, prev_vEnd, prev_aEnd, prev_jEnd;
                segment_eval(prevCreatedSeg, prevCreatedSeg->tEnd,
                             &prev_xEnd, &prev_vEnd, &prev_aEnd, &prev_jEnd);
                seg->x0 = prev_xEnd;
                seg->v0 = prev_vEnd;
                seg->a0 = prev_aEnd;

                // Set initial jerk based on segment type:
                switch (seg->segmentIndex)
                {
                case 2:
                    seg->j0 = 0.0f;
                    break;
                case 3:
                    seg->j0 = -obj->param.sign * obj->jerk_max;
                    break;
                case 4:
                    seg->j0 = 0.0f;
                    break;
                case 5:
                    seg->j0 = -obj->param.sign * obj->jerk_max;
                    break;
                case 6:
                    seg->j0 = 0.0f;
                    break;
                case 7:
                    seg->j0 = obj->param.sign * obj->jerk_max;
                    break;
                default:
                    seg->j0 = 0.0f;
                    break;
                }

                // Link the previous segment with this new segment.
                prevCreatedSeg->next = seg;
                prevCreatedSeg = seg;
            }

            tCurrent = tEnd;
        }

        obj->segStateIndex = 0;
        obj->profile_ticks = 0;
        obj->state = CURVE_BUSY;
        // break;
    }

    case CURVE_BUSY:
    {
        float t_sec = (float)obj->profile_ticks++ * obj->Ts;
        scurve_segment_state_t *seg = &obj->segState[obj->segStateIndex];

        if (t_sec >= seg->tEnd)
        {
            if (seg->next == NULL)
            {
                // finalize
                float p, v, a, j;
                segment_eval(seg, obj->param.totalTime, &p, &v, &a, &j);

                // Snap outputs to that final position/velocity/accel/jerk:
                if (obj->param.sign * (p - obj->pos_target) > 0.0f)
                {
                    p = obj->pos_target;
                }
                obj->pos_out = p;
                // obj->pos_out = obj->pos_target;
                obj->vel_out = v;
                obj->acc_out = a;
                obj->jerk_out = j;

                // End of profile
                obj->state = CURVE_ONEND;
                break;
            }
            // move to next segment
            obj->segStateIndex++;
            seg = seg->next;
        }

        // Compute the profile at the current time
        float p, v, a, j;
        segment_eval(seg, t_sec, &p, &v, &a, &j);

        // if overshoot, snap to target
        // if target is positive, overshoot is positive
        // if target is negative, overshoot is negative
        if (obj->param.sign * (p - obj->pos_target) >= 0.0f)
        {
            p = obj->pos_target;
            obj->state = CURVE_ONEND;
        }

        obj->pos_out = p;
        obj->vel_out = v;
        obj->acc_out = a;
        obj->jerk_out = j;
        break;
    }

    case CURVE_ONEND:
    {
        obj->pos_out = obj->pos_target;
#ifdef ON_END_RESET_VAJ
        obj->vel_out = 0.0f;
        obj->acc_out = 0.0f;
        obj->jerk_out = 0.0f;
#endif
        obj->state = CURVE_IDLE;
        break;
    }

    default:
        break;
    }
    return obj->state;
}

/**
 * Starts execution of a new motion profile
 * 
 * Validates parameters and transitions from IDLE to CONF state.
 * The profile will only start if:
 * - Current state is IDLE
 * - Maximum velocity > start and stop velocities
 * - All limits (velocity, acceleration, deceleration, jerk) are positive
 * 
 * @param obj Pointer to scurve object
 * @return CURVE_SUCCESS if profile started successfully,
 *         CURVE_INVALID_STATE if not in IDLE state,
 *         CURVE_INVALID_PARAMS if parameters are invalid
 */
scurve_return_status_e scurve_startProfile(scurve_t *obj)
{
    if (obj->state != CURVE_IDLE)
    {
        return CURVE_INVALID_STATE;
    }

    if (obj->vel_max <= obj->vel_start || obj->vel_max <= obj->vel_stop)
    {
        return CURVE_INVALID_PARAMS;
    }

    if (obj->vel_max <= 0.0f || obj->acc_max <= 0.0f || obj->dec_max <= 0.0f || obj->jerk_max <= 0.0f)
    {
        return CURVE_INVALID_PARAMS;
    }

    obj->state = CURVE_CONF;
    return CURVE_SUCCESS;
}

scurve_state_e scurve_getState(scurve_t *obj)
{
    return obj->state;
}

void scurve_getOutputs(scurve_t *obj, float *pos, float *vel, float *acc, float *jrk)
{
    *pos = obj->pos_out;
    *vel = obj->vel_out;
    *acc = obj->acc_out;
    *jrk = obj->jerk_out;
}

void scurve_setSampleTime(scurve_t *obj, float sampleTime_s)
{
    if (obj->state == CURVE_IDLE)
        obj->Ts = sampleTime_s;
}

/**
 * Sets target position for the motion profile
 * Only allowed when in IDLE state
 */
void scurve_setPositionTarget(scurve_t *obj, float pos_target)
{
    if (obj->state == CURVE_IDLE)
        obj->pos_target = pos_target;
}

/**
 * Sets current position output
 * Only allowed when in IDLE state
 */
void scurve_setPositionOutput(scurve_t *obj, float pos_out)
{
    if (obj->state == CURVE_IDLE)
        obj->pos_out = pos_out;
}

/**
 * Sets initial velocity for the motion profile
 * Only allowed when in IDLE state
 */
void scurve_setVelocityStart(scurve_t *obj, float vel_start)
{
    if (obj->state == CURVE_IDLE)
        obj->vel_start = vel_start;
}

/**
 * Sets final velocity for the motion profile
 * Only allowed when in IDLE state
 */
void scurve_setVelocityStop(scurve_t *obj, float vel_stop)
{
    if (obj->state == CURVE_IDLE)
        obj->vel_stop = vel_stop;
}

/**
 * Sets maximum allowed velocity
 * Only allowed when in IDLE state
 * Must be greater than both start and stop velocities
 */
void scurve_setVelocityMax(scurve_t *obj, float vel)
{
    if (obj->state == CURVE_IDLE)
        obj->vel_max = vel;
}

/**
 * Sets maximum allowed acceleration
 * Only allowed when in IDLE state
 * Must be positive
 */
void scurve_setAccelerationMax(scurve_t *obj, float acc)
{
    if (obj->state == CURVE_IDLE)
        obj->acc_max = acc;
}

/**
 * Sets maximum allowed deceleration
 * Only allowed when in IDLE state
 * Must be positive
 */
void scurve_setDecelerationMax(scurve_t *obj, float dec)
{
    if (obj->state == CURVE_IDLE)
        obj->dec_max = dec;
}

/**
 * Sets maximum allowed jerk
 * Only allowed when in IDLE state
 * Must be positive
 */
void scurve_setJerkMax(scurve_t *obj, float jerk)
{
    if (obj->state == CURVE_IDLE)
        obj->jerk_max = jerk;
}

float scurve_getPositionOutput(scurve_t *obj)
{
    return obj->pos_out;
}

float scurve_getVelocityOutput(scurve_t *obj)
{
    return obj->vel_out;
}

float scurve_getAccelerationOutput(scurve_t *obj)
{
    return obj->acc_out;
}

float scurve_getJerkOutput(scurve_t *obj)
{
    return obj->jerk_out;
}

uint32_t scurve_getProfileTick(scurve_t *obj)
{
    return obj->profile_ticks;
}