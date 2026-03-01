/*
 * integrator_perf_test.cpp
 *
 * Compares integration step counts and wall-clock runtime between the RK4
 * and RK45/DOPRI5 integrators for an identical set of rays (lamppost point
 * source → Kerr disc).
 *
 * Metrics reported
 * ----------------
 *   - Wall-clock time for the propagation phase only (excludes PointSource
 *     setup and redshift calculation, which are common to both methods).
 *   - Per-ray step count statistics: min, mean, median, 75th/95th percentile,
 *     and max.
 *   - Estimated total function evaluations: RK4 costs 4 evaluations per step;
 *     DOPRI5 costs 6 per step (the implementation is FSAL so the 7th stage
 *     of step n reuses k1 of step n+1 — but that saving is not exploited
 *     here, so 6 is the correct count).
 *   - ASCII histogram of the step-count distribution.
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>
#include <chrono>
using namespace std;

#include "raytracer/pointsource.h"
#include "include/kerr.h"

// ---------------------------------------------------------------------------
// Source / grid parameters
// ---------------------------------------------------------------------------
static const double SPIN        = 0.998;
static const double SOURCE[]    = {0.0, 5.0, 1e-3, 0.0};
static const double V           = 0.0;
static const double COSALPHA0   = -0.995;
static const double COSALPHAMAX =  0.995;
static const double DCOSALPHA   =  0.05;
static const double BETA0       = -M_PI;
static const double BETAMAX     =  M_PI;
static const double DBETA       =  0.05;
static const double R_MAX       = 1000.0;
static const int    SHOW_PROGRESS = 0;

// DOPRI5 stages per step (6 new evaluations; 7th stage is FSAL but not used
// as such in this implementation)
static const int FEVALS_PER_STEP_RK4  = 4;
static const int FEVALS_PER_STEP_RK45 = 6;

// ---------------------------------------------------------------------------
struct StepStats
{
    long   n_valid;
    long   n_invalid;
    long   total_steps;
    double mean;
    double pct25, pct50, pct75, pct95;
    int    min_steps, max_steps;
    double wall_seconds;
    vector<int> sorted;   // sorted step counts of valid rays
};

// ---------------------------------------------------------------------------
static StepStats run_and_collect(bool use_rk45)
{
    double pos[4] = {SOURCE[0], SOURCE[1], SOURCE[2], SOURCE[3]};
    PointSource<double> src(pos, V, SPIN, TOL,
                            DCOSALPHA, DBETA, COSALPHA0, COSALPHAMAX, BETA0, BETAMAX);

    src.redshift_start();   // initialise momenta — not timed

    auto t0 = chrono::high_resolution_clock::now();
    if (use_rk45)
        src.run_raytrace_rk45(R_MAX, M_PI_2, SHOW_PROGRESS);
    else
        src.run_raytrace_rk4 (R_MAX, M_PI_2, SHOW_PROGRESS);
    auto t1 = chrono::high_resolution_clock::now();

    StepStats s{};
    s.wall_seconds = chrono::duration<double>(t1 - t0).count();

    int n = src.get_count();
    s.sorted.reserve(n);
    for (int i = 0; i < n; i++)
    {
        int st = src.rays[i].steps;
        if (st <= 0) { ++s.n_invalid; continue; }
        ++s.n_valid;
        s.total_steps += st;
        s.sorted.push_back(st);
    }
    sort(s.sorted.begin(), s.sorted.end());

    if (!s.sorted.empty())
    {
        s.min_steps = s.sorted.front();
        s.max_steps = s.sorted.back();
        s.mean      = (double)s.total_steps / s.n_valid;
        auto pct = [&](double p) -> double {
            double idx = p * (s.sorted.size() - 1);
            int lo = (int)idx;
            int hi = min(lo + 1, (int)s.sorted.size() - 1);
            double frac = idx - lo;
            return s.sorted[lo] * (1 - frac) + s.sorted[hi] * frac;
        };
        s.pct25 = pct(0.25);
        s.pct50 = pct(0.50);
        s.pct75 = pct(0.75);
        s.pct95 = pct(0.95);
    }
    return s;
}

// ---------------------------------------------------------------------------
// Print an ASCII histogram aligned in two columns (rk4 | rk45)
// ---------------------------------------------------------------------------
static void print_histogram(const StepStats& rk4, const StepStats& rk45,
                             int nbins = 12, int bar_width = 30)
{
    int lo = min(rk4.min_steps, rk45.min_steps);
    int hi = max(rk4.max_steps, rk45.max_steps);
    if (lo >= hi) return;

    // log-spaced bins so we cover the wide range nicely
    double log_lo = log10(max(lo, 1));
    double log_hi = log10(hi) + 1e-9;
    double dlog   = (log_hi - log_lo) / nbins;

    auto bin_of = [&](int v) -> int {
        if (v <= 0) return -1;
        int b = (int)((log10(v) - log_lo) / dlog);
        return max(0, min(nbins - 1, b));
    };

    vector<int> cnt4(nbins, 0), cnt45(nbins, 0);
    for (int v : rk4.sorted)   cnt4 [bin_of(v)]++;
    for (int v : rk45.sorted)  cnt45[bin_of(v)]++;

    int peak = 0;
    for (int b = 0; b < nbins; b++)
        peak = max(peak, max(cnt4[b], cnt45[b]));

    cout << "\nStep-count distribution (log10-spaced bins):\n";
    cout << setw(18) << "range"
         << setw(8) << "N_rk4"
         << "  " << left << setw(bar_width) << "RK4 bar"
         << setw(8) << "N_rk45"
         << "  RK45 bar" << right << "\n";
    cout << string(18 + 8 + 2 + bar_width + 8 + 2 + bar_width, '-') << "\n";

    for (int b = 0; b < nbins; b++)
    {
        double lo_b = pow(10, log_lo + b * dlog);
        double hi_b = pow(10, log_lo + (b + 1) * dlog);
        int len4  = (peak > 0) ? cnt4 [b] * bar_width / peak : 0;
        int len45 = (peak > 0) ? cnt45[b] * bar_width / peak : 0;

        cout << fixed << setprecision(0)
             << setw(8) << lo_b << " -"
             << setw(8) << hi_b
             << setw(8) << cnt4[b]
             << "  " << left
             << string(len4, '#') << string(bar_width - len4, ' ')
             << setw(8) << cnt45[b]
             << "  " << string(len45, '#') << right << "\n";
    }
}

// ---------------------------------------------------------------------------
int main()
{
    cout << fixed << setprecision(3);
    cout << "Integrator performance comparison\n";
    cout << "spin = " << SPIN
         << ",  source r = " << SOURCE[1]
         << ",  source theta = " << SOURCE[2] << "\n";
    cout << "angular grid: dcosalpha = " << DCOSALPHA
         << ",  dbeta = " << DBETA << "\n\n";

    cout << "Running RK4...\n";
    StepStats rk4 = run_and_collect(false);

    cout << "\nRunning RK45...\n";
    StepStats rk45 = run_and_collect(true);

    // -----------------------------------------------------------------------
    // Summary table
    // -----------------------------------------------------------------------
    cout << "\n";
    cout << string(52, '=') << "\n";
    cout << left << setw(30) << "Metric"
         << right << setw(10) << "RK4"
         << setw(12) << "RK45" << "\n";
    cout << string(52, '-') << "\n";

    long n_total = rk4.n_valid + rk4.n_invalid;
    cout << left << setw(30) << "Total rays"
         << right << setw(10) << n_total
         << setw(12) << (rk45.n_valid + rk45.n_invalid) << "\n";

    cout << left << setw(30) << "Valid rays"
         << right << setw(10) << rk4.n_valid
         << setw(12) << rk45.n_valid << "\n";

    cout << left << setw(30) << "Invalid rays (skipped)"
         << right << setw(10) << rk4.n_invalid
         << setw(12) << rk45.n_invalid << "\n";

    cout << string(52, '-') << "\n";

    cout << left << setw(30) << "Wall time (s)"
         << right << setw(10) << setprecision(3) << rk4.wall_seconds
         << setw(12) << rk45.wall_seconds << "\n";

    if (rk4.wall_seconds > 0)
    {
        double ratio = rk45.wall_seconds / rk4.wall_seconds;
        cout << left << setw(30) << "RK45 / RK4 time ratio"
             << right << setw(10) << ""
             << setw(12) << setprecision(2) << ratio << "x\n";
    }

    cout << string(52, '-') << "\n";
    cout << setprecision(1);

    cout << left << setw(30) << "Steps per ray: min"
         << right << setw(10) << rk4.min_steps
         << setw(12) << rk45.min_steps << "\n";
    cout << left << setw(30) << "Steps per ray: mean"
         << right << setw(10) << rk4.mean
         << setw(12) << rk45.mean << "\n";
    cout << left << setw(30) << "Steps per ray: median (50th)"
         << right << setw(10) << rk4.pct50
         << setw(12) << rk45.pct50 << "\n";
    cout << left << setw(30) << "Steps per ray: 75th pct"
         << right << setw(10) << rk4.pct75
         << setw(12) << rk45.pct75 << "\n";
    cout << left << setw(30) << "Steps per ray: 95th pct"
         << right << setw(10) << rk4.pct95
         << setw(12) << rk45.pct95 << "\n";
    cout << left << setw(30) << "Steps per ray: max"
         << right << setw(10) << rk4.max_steps
         << setw(12) << rk45.max_steps << "\n";

    cout << string(52, '-') << "\n";

    cout << left << setw(30) << "Total steps"
         << right << setw(10) << rk4.total_steps
         << setw(12) << rk45.total_steps << "\n";

    long fevals4  = rk4.total_steps  * FEVALS_PER_STEP_RK4;
    long fevals45 = rk45.total_steps * FEVALS_PER_STEP_RK45;
    cout << left << setw(30) << "Total func evals (est.)"
         << right << setw(10) << fevals4
         << setw(12) << fevals45 << "\n";

    if (fevals4 > 0)
    {
        cout << left << setw(30) << "Func evals ratio RK45/RK4"
             << right << setw(10) << ""
             << setw(12) << setprecision(2) << (double)fevals45 / fevals4 << "x\n";
    }

    if (rk4.wall_seconds > 0 && rk45.wall_seconds > 0)
    {
        double thr4  = fevals4  / rk4.wall_seconds  / 1e6;
        double thr45 = fevals45 / rk45.wall_seconds / 1e6;
        cout << left << setw(30) << "Throughput (Mfevals/s): RK4"
             << right << setw(10) << setprecision(1) << thr4
             << setw(12) << "" << "\n";
        cout << left << setw(30) << "Throughput (Mfevals/s): RK45"
             << right << setw(10) << ""
             << setw(12) << thr45 << "\n";
    }

    cout << string(52, '=') << "\n";

    print_histogram(rk4, rk45);

    return 0;
}
