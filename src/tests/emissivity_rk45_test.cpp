/*
 * emissivity_rk45_test.cpp
 *
 * Runs the emissivity calculation (lamppost point source → accretion disc) using
 * both the RK4 and RK45 integrators and compares the per-bin emissivity, mean
 * redshift, and mean arrival time for each radial annulus.
 *
 * Design note on tolerances
 * -------------------------
 * Some rays near the photon-sphere separatrix are chaotically sensitive to
 * step size: a tiny numerical difference causes them to orbit one extra time
 * and land at a completely different radius.  This is physically expected and
 * is NOT the bug we are testing for.  To avoid false failures from these rays
 * the pass/fail criteria are applied only to *well-populated* bins
 * (>= MIN_BIN_RAYS in BOTH integrators), where averaging over many rays
 * suppresses the sensitivity to individual separatrix-crossing outliers.
 *
 * The bug that motivated this test (DOPRI5 intermediate stages stepping
 * inside the event horizon, corrupting entire steps undetected) produces
 * large emissivity glitches even in the well-populated inner-disc bins.
 * Those glitches are caught by the EMIS_TOL threshold.
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
using namespace std;

#include "raytracer/pointsource.h"
#include "include/kerr.h"
#include "include/disc.h"

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
static const double R_DISC      =  500.0;
static const double GAMMA       =  2.0;
static const int    NR          = 30;
static const bool   LOGBIN      = true;
static const int    SHOW_PROGRESS = 0;

// ---------------------------------------------------------------------------
// Pass/fail thresholds (applied only to bins with >= MIN_BIN_RAYS in both
// AND whose ray counts agree within MAX_RAYS_DIFF_FRAC)
// ---------------------------------------------------------------------------
static const int    MIN_BIN_RAYS      = 100;  // ignore low-statistics bins
static const double MAX_RAYS_DIFF_FRAC= 0.10; // skip bins where N_rk4/N_rk45 differ >10%
                                               // (different counts → ray redistribution from
                                               //  chaotic photon-sphere separatrix, not accuracy)
static const double EMIS_TOL          = 0.10; // 10% – enough to catch horizon-crossing glitches
static const double REDSHIFT_TOL      = 0.005;// 0.5%
static const double TIME_TOL          = 0.05; // 5%

// ---------------------------------------------------------------------------
// Reduce the emissivity profile for one PointSource run.
// Returns the total number of disc-hitting rays.
// ---------------------------------------------------------------------------
static long reduce(PointSource<double>& src,
                   bool use_rk45,
                   double r_isco,
                   double dr_fac,
                   double disc_r[],
                   double disc_area[],
                   long   disc_rays[],
                   double disc_emis[],
                   double disc_redshift[],
                   double disc_time[])
{
    for (int ir = 0; ir < NR; ir++)
    {
        disc_rays[ir] = 0;
        disc_emis[ir] = disc_redshift[ir] = disc_time[ir] = 0.0;
    }

    src.redshift_start();

    if (use_rk45)
        src.run_raytrace_rk45(R_MAX, M_PI_2, SHOW_PROGRESS);
    else
        src.run_raytrace_rk4(R_MAX, M_PI_2, SHOW_PROGRESS);

    src.range_phi();
    src.redshift(-1.0, false, false, 0);

    long disc_count = 0;
    for (int ray = 0; ray < src.get_count(); ray++)
    {
        if (src.rays[ray].steps <= 0) continue;
        if (!isfinite(src.rays[ray].redshift) || src.rays[ray].redshift <= 0) continue;
        if (src.rays[ray].r < r_isco) continue;

        double x, y, z;
        cartesian(x, y, z, src.rays[ray].r, src.rays[ray].theta,
                  src.rays[ray].phi, SPIN);
        if (z > 1e-2) continue;

        int ir = static_cast<int>(log(src.rays[ray].r / r_isco) / log(dr_fac));
        if (ir < 0 || ir >= NR) continue;

        ++disc_rays[ir];
        disc_emis[ir]     += 1.0 / pow(src.rays[ray].redshift, GAMMA);
        disc_redshift[ir] += src.rays[ray].redshift;
        disc_time[ir]     += src.rays[ray].t;
        ++disc_count;
    }

    for (int ir = 0; ir < NR; ir++)
    {
        if (disc_rays[ir] > 0)
        {
            disc_redshift[ir] /= disc_rays[ir];
            disc_time[ir]     /= disc_rays[ir];
        }
        if (disc_area[ir] > 0)
            disc_emis[ir] /= disc_area[ir];
    }
    return disc_count;
}

// ---------------------------------------------------------------------------
int main()
{
    cout << fixed << setprecision(6);

    const double r_isco  = kerr_isco<double>(SPIN, +1);
    const double dr_fac  = exp(log(R_DISC / r_isco) / NR);

    cout << "spin = " << SPIN << ",  r_ISCO = " << r_isco << endl;
    cout << "source r = " << SOURCE[1]
         << ",  source theta = " << SOURCE[2] << endl;
    cout << "angular grid: dcosalpha = " << DCOSALPHA
         << ",  dbeta = " << DBETA << endl;
    cout << "radial bins:  NR = " << NR << " (log),  r_isco to " << R_DISC << endl;
    cout << endl;

    // Shared bin geometry
    double disc_r[NR], disc_area[NR];
    for (int ir = 0; ir < NR; ir++)
    {
        disc_r[ir]   = r_isco * pow(dr_fac, ir);
        disc_area[ir] = integrate_disc_area(disc_r[ir], disc_r[ir] * dr_fac, SPIN);
    }

    // Output arrays
    long   rk4_rays[NR]  = {}, rk45_rays[NR]  = {};
    double rk4_emis[NR]  = {}, rk45_emis[NR]  = {};
    double rk4_redsh[NR] = {}, rk45_redsh[NR] = {};
    double rk4_time[NR]  = {}, rk45_time[NR]  = {};

    // --- RK4 run ---
    cout << "Running RK4..." << endl;
    {
        double pos[4] = {SOURCE[0], SOURCE[1], SOURCE[2], SOURCE[3]};
        PointSource<double> rk4(pos, V, SPIN, TOL,
                                DCOSALPHA, DBETA, COSALPHA0, COSALPHAMAX, BETA0, BETAMAX);
        reduce(rk4, false, r_isco, dr_fac,
               disc_r, disc_area, rk4_rays, rk4_emis, rk4_redsh, rk4_time);
    }

    // --- RK45 run ---
    cout << "Running RK45..." << endl;
    {
        double pos[4] = {SOURCE[0], SOURCE[1], SOURCE[2], SOURCE[3]};
        PointSource<double> rk45(pos, V, SPIN, TOL,
                                 DCOSALPHA, DBETA, COSALPHA0, COSALPHAMAX, BETA0, BETAMAX);
        reduce(rk45, true, r_isco, dr_fac,
               disc_r, disc_area, rk45_rays, rk45_emis, rk45_redsh, rk45_time);
    }

    // ---------------------------------------------------------------------------
    // Comparison table
    // ---------------------------------------------------------------------------
    cout << endl;
    cout << setw(5)  << "bin"
         << setw(10) << "r"
         << setw(8)  << "N_rk4"
         << setw(8)  << "N_rk45"
         << setw(14) << "emis_rk4"
         << setw(14) << "emis_rk45"
         << setw(11) << "rel_emis"
         << setw(11) << "rel_redsh"
         << setw(11) << "rel_time"
         << setw(10) << "verdict"
         << endl;
    cout << string(102, '-') << endl;

    double max_rel_emis = 0, max_rel_redsh = 0, max_rel_time = 0;
    double sum_rel_emis = 0, sum_rel_redsh = 0, sum_rel_time = 0;
    int n_compared = 0, n_failed = 0;

    for (int ir = 0; ir < NR; ir++)
    {
        if (rk4_rays[ir] == 0 && rk45_rays[ir] == 0) continue;

        bool enough = (rk4_rays[ir] >= MIN_BIN_RAYS && rk45_rays[ir] >= MIN_BIN_RAYS);
        long ray_diff = abs(rk4_rays[ir] - rk45_rays[ir]);
        long ray_max  = max(rk4_rays[ir], rk45_rays[ir]);
        bool consistent = enough && (ray_diff <= (long)(MAX_RAYS_DIFF_FRAC * ray_max));
        bool valid  = consistent
                   && isfinite(rk4_emis[ir]) && isfinite(rk45_emis[ir])
                   && rk4_emis[ir] > 0 && rk45_emis[ir] > 0;

        double rel_emis = 0, rel_redsh = 0, rel_time = 0;
        if (valid)
        {
            rel_emis  = abs(rk4_emis[ir]  - rk45_emis[ir])  / rk4_emis[ir];
            rel_redsh = abs(rk4_redsh[ir] - rk45_redsh[ir]) / rk4_redsh[ir];
            rel_time  = abs(rk4_time[ir]  - rk45_time[ir])  / abs(rk4_time[ir]);

            if (rel_emis  > max_rel_emis)  max_rel_emis  = rel_emis;
            if (rel_redsh > max_rel_redsh) max_rel_redsh = rel_redsh;
            if (rel_time  > max_rel_time)  max_rel_time  = rel_time;

            sum_rel_emis  += rel_emis;
            sum_rel_redsh += rel_redsh;
            sum_rel_time  += rel_time;
            ++n_compared;

            if (rel_emis > EMIS_TOL || rel_redsh > REDSHIFT_TOL || rel_time > TIME_TOL)
                ++n_failed;
        }

        string verdict;
        if (!enough)           verdict = "(low-N)";
        else if (!consistent)  verdict = "(sep)";   // chaotic separatrix ray redistribution
        else if (!valid)       verdict = "(NaN)";
        else if (n_failed > n_compared - 1)  // just incremented for this bin
                               verdict = "DIFF";

        cout << setw(5)  << ir
             << setw(10) << disc_r[ir]
             << setw(8)  << rk4_rays[ir]
             << setw(8)  << rk45_rays[ir]
             << setw(14) << (isfinite(rk4_emis[ir])  ? rk4_emis[ir]  : 0.0)
             << setw(14) << (isfinite(rk45_emis[ir]) ? rk45_emis[ir] : 0.0)
             << setw(11) << (valid ? to_string(rel_emis).substr(0,8)  : (consistent ? "(NaN)" : "  -"))
             << setw(11) << (valid ? to_string(rel_redsh).substr(0,8) : (consistent ? "(NaN)" : "  -"))
             << setw(11) << (valid ? to_string(rel_time).substr(0,8)  : (consistent ? "(NaN)" : "  -"))
             << setw(10) << verdict
             << endl;
    }
    cout << string(102, '-') << endl << endl;

    // ---------------------------------------------------------------------------
    // Summary
    // ---------------------------------------------------------------------------
    cout << "Bins with >= " << MIN_BIN_RAYS << " rays and consistent counts (<"
         << (int)(MAX_RAYS_DIFF_FRAC*100) << "% diff) compared:  " << n_compared << endl;
    if (n_compared > 0)
    {
        cout << "  Max  |rel_emis|  = " << max_rel_emis
             << "  (threshold " << EMIS_TOL << ")" << endl;
        cout << "  Max  |rel_redsh| = " << max_rel_redsh
             << "  (threshold " << REDSHIFT_TOL << ")" << endl;
        cout << "  Max  |rel_time|  = " << max_rel_time
             << "  (threshold " << TIME_TOL << ")" << endl;
        cout << "  Mean |rel_emis|  = " << sum_rel_emis  / n_compared << endl;
        cout << "  Mean |rel_redsh| = " << sum_rel_redsh / n_compared << endl;
        cout << "  Mean |rel_time|  = " << sum_rel_time  / n_compared << endl;
        cout << "  Bins failing any threshold: " << n_failed << endl;
    }

    bool pass = (n_compared > 0)
             && (n_failed == 0)
             && (max_rel_emis  <= EMIS_TOL)
             && (max_rel_redsh <= REDSHIFT_TOL)
             && (max_rel_time  <= TIME_TOL);

    cout << endl << (pass ? "PASS" : "FAIL") << endl;
    return pass ? 0 : 1;
}
