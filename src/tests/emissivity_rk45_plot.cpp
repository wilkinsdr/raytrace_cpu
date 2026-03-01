/*
 * emissivity_rk45_plot.cpp
 *
 * Runs RK4 and RK45 emissivity calculations and writes per-bin data to a CSV
 * file for plotting.  Intended to be invoked by emissivity_rk45_plot.py.
 *
 * Usage:  emissivity_rk45_plot [output_csv_path]
 *
 * If output_csv_path is omitted the data is written to emissivity_rk45_data.csv
 * in the current working directory.  All library progress messages go to stdout
 * as usual; only the CSV goes to the specified file.
 *
 * CSV columns:
 *   r  N_rk4  N_rk45  emis_rk4  emis_rk45
 *   redsh_rk4  redsh_rk45  time_rk4  time_rk45
 */

#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#include "raytracer/pointsource.h"
#include "include/kerr.h"
#include "include/disc.h"

// ---------------------------------------------------------------------------
// Source / grid parameters  (match emissivity_rk45_test.cpp)
// ---------------------------------------------------------------------------
static const double SPIN        = 0.998;
static const double SOURCE[]    = {0.0, 5.0, 1e-3, 0.0};
static const double V           = 0.0;
static const double COSALPHA0   = -0.995;
static const double COSALPHAMAX =  0.995;
static const double DCOSALPHA   =  0.01;
static const double BETA0       = -M_PI;
static const double BETAMAX     =  M_PI;
static const double DBETA       =  0.01;
static const double R_MAX       = 1000.0;
static const double R_DISC      =  500.0;
static const double GAMMA       =  2.0;
static const int    NR          = 30;
static const int    SHOW_PROGRESS = 0;

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
int main(int argc, char* argv[])
{
    const char* outpath      = (argc > 1) ? argv[1]       : "emissivity_rk45_data.csv";
    const double rk45_tol_arg = (argc > 2) ? atof(argv[2]) : 1e-8;

    const double r_isco  = kerr_isco<double>(SPIN, +1);
    const double dr_fac  = exp(log(R_DISC / r_isco) / NR);

    double disc_r[NR], disc_area[NR];
    for (int ir = 0; ir < NR; ir++)
    {
        disc_r[ir]    = r_isco * pow(dr_fac, ir);
        disc_area[ir] = integrate_disc_area(disc_r[ir], disc_r[ir] * dr_fac, SPIN);
    }

    long   rk4_rays[NR]  = {}, rk45_rays[NR]  = {};
    double rk4_emis[NR]  = {}, rk45_emis[NR]  = {};
    double rk4_redsh[NR] = {}, rk45_redsh[NR] = {};
    double rk4_time[NR]  = {}, rk45_time[NR]  = {};

    cout << "Running RK4..." << endl;
    {
        double pos[4] = {SOURCE[0], SOURCE[1], SOURCE[2], SOURCE[3]};
        PointSource<double> rk4(pos, V, SPIN, TOL,
                                DCOSALPHA, DBETA, COSALPHA0, COSALPHAMAX, BETA0, BETAMAX);
        reduce(rk4, false, r_isco, dr_fac,
               disc_r, disc_area, rk4_rays, rk4_emis, rk4_redsh, rk4_time);
    }

    cout << "Running RK45 (tol=" << rk45_tol_arg << ")..." << endl;
    {
        double pos[4] = {SOURCE[0], SOURCE[1], SOURCE[2], SOURCE[3]};
        PointSource<double> rk45(pos, V, SPIN, TOL,
                                 DCOSALPHA, DBETA, COSALPHA0, COSALPHAMAX, BETA0, BETAMAX);
        rk45.set_rk45_tol(rk45_tol_arg);
        reduce(rk45, true, r_isco, dr_fac,
               disc_r, disc_area, rk45_rays, rk45_emis, rk45_redsh, rk45_time);
    }

    ofstream fout(outpath);
    if (!fout)
    {
        cerr << "ERROR: cannot open output file: " << outpath << endl;
        return 1;
    }

    fout << "# rk45_tol=" << rk45_tol_arg << "\n";
    fout << "r N_rk4 N_rk45 emis_rk4 emis_rk45 redsh_rk4 redsh_rk45 time_rk4 time_rk45\n";
    fout.precision(8);
    fout << scientific;
    for (int ir = 0; ir < NR; ir++)
    {
        if (rk4_rays[ir] == 0 && rk45_rays[ir] == 0) continue;
        fout << disc_r[ir]
             << " " << rk4_rays[ir]
             << " " << rk45_rays[ir]
             << " " << (isfinite(rk4_emis[ir])   ? rk4_emis[ir]   : 0.0)
             << " " << (isfinite(rk45_emis[ir])  ? rk45_emis[ir]  : 0.0)
             << " " << (isfinite(rk4_redsh[ir])  ? rk4_redsh[ir]  : 0.0)
             << " " << (isfinite(rk45_redsh[ir]) ? rk45_redsh[ir] : 0.0)
             << " " << (isfinite(rk4_time[ir])   ? rk4_time[ir]   : 0.0)
             << " " << (isfinite(rk45_time[ir])  ? rk45_time[ir]  : 0.0)
             << "\n";
    }

    cout << "Data written to " << outpath << endl;
    return 0;
}
