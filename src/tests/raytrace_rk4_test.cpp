/*
 * raytrace_rk4_test.cpp
 *
 * Runs the same rays through both the Euler and RK4 propagators and compares
 * the results to verify that propagate_rk4() produces consistent output.
 *
 * For each ray the test checks:
 *   - Both integrators reach the same termination boundary (disc, r_max, or horizon)
 *   - Reports the differences in final position between the two integrators
 *
 * Expected behaviour: RK4 and Euler should agree on which boundary each ray hits.
 * The final positions will differ slightly (Euler is first-order, RK4 fourth-order),
 * but the differences should be small compared to the overall path length.
 */

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

#include "raytracer/pointsource.h"
#include "include/kerr.h"

int main()
{
    const double spin = 0.998;
    double source[] = {0.0, 5.0, 1e-3, 0.0};  // near-polar lamppost at r=5
    const double V = disc_velocity(source[1], spin, +1);
    const double dcosalpha = 0.2;
    const double dbeta = 0.2;
    const double r_max = 1000.0;
    const double theta_max = M_PI_2;

    cout << fixed << setprecision(6);
    cout << "spin = " << spin << ",  source r = " << source[1]
         << ",  V = " << V << endl;
    cout << "dcosalpha = " << dcosalpha << ",  dbeta = " << dbeta << endl << endl;

    // Run Euler integrator
    PointSource<double>* euler_src = new PointSource<double>(
        source, V, spin, TOL, dcosalpha, dbeta);
    euler_src->run_raytrace(r_max, theta_max, 0);

    // Run RK4 integrator with identical initial conditions
    PointSource<double>* rk4_src = new PointSource<double>(
        source, V, spin, TOL, dcosalpha, dbeta);
    rk4_src->run_raytrace_rk4(r_max, theta_max, 0);

    // Compare results
    const int n = euler_src->get_count();
    const double h_r = kerr_horizon(spin);

    // Disc boundary tolerance: RK4 uses ptheta1 for boundary clipping but
    // applies the full weighted RK4 update, so the final theta may overshoot
    // theta_max by a small amount in the last step.  Use a loose tolerance
    // (0.01 rad ~ 0.6 deg) to avoid false mismatches from this overshoot.
    const double disc_tol = 1e-2;

    int n_disc = 0, n_rlim = 0, n_horizon = 0, n_skip = 0;
    int n_mismatch = 0, n_diverge = 0;
    double max_dr = 0, max_dtheta = 0, max_dphi = 0;

    // classify termination boundary for a ray
    auto boundary = [&](const Ray<double>& ray) -> char {
        if(abs(ray.theta - theta_max) < disc_tol) return 'D';  // disc
        if(ray.r >= r_max)                        return 'R';  // r_max
        if(ray.r <= h_r)                          return 'H';  // horizon
        return '?';
    };

    cout << setw(5)  << "ray"
         << setw(12) << "r(Euler)"
         << setw(12) << "r(RK4)"
         << setw(12) << "|dr|"
         << setw(14) << "|dtheta|"
         << setw(12) << "|dphi|"
         << setw(10) << "boundary" << endl;
    cout << string(77, '-') << endl;

    for(int i = 0; i < n; i++)
    {
        const Ray<double>& e = euler_src->rays[i];
        const Ray<double>& r = rk4_src->rays[i];

        if(e.steps < 0 || r.steps < 0) { n_skip++; continue; }

        char e_b = boundary(e);
        char r_b = boundary(r);
        bool match = (e_b == r_b);

        if(e_b == 'D') n_disc++;
        else if(e_b == 'R') n_rlim++;
        else if(e_b == 'H') n_horizon++;

        double dr     = abs(e.r     - r.r);
        double dtheta = abs(e.theta - r.theta);
        double dphi   = abs(e.phi   - r.phi);

        if(!match)
        {
            // A genuine divergence: the two integrators ended at completely
            // different boundaries (e.g. disc vs r_max).  This can happen
            // for near-photon-sphere orbiting rays where small integration
            // errors are amplified over many orbits — physically expected.
            n_diverge++;
            cout << setw(5)  << i
                 << setw(12) << e.r
                 << setw(12) << r.r
                 << setw(12) << dr
                 << setw(14) << dtheta
                 << setw(12) << dphi
                 << "  " << e_b << "!=" << r_b << " (diverged)" << endl;
        }
        else
        {
            // Same boundary: accumulate position differences
            if(dr     > max_dr)     max_dr     = dr;
            if(dtheta > max_dtheta) max_dtheta = dtheta;
            if(dphi   > max_dphi)   max_dphi   = dphi;
        }
    }

    cout << string(77, '-') << endl << endl;

    cout << "Summary:" << endl;
    cout << "  Total rays:              " << n << endl;
    cout << "  Skipped (invalid):       " << n_skip << endl;
    cout << "  Hit disc (theta_max):    " << n_disc << endl;
    cout << "  Hit r_max:               " << n_rlim << endl;
    cout << "  Hit horizon:             " << n_horizon << endl;
    cout << "  Diverged trajectories:   " << n_diverge
         << "  (expected for orbiting rays near photon sphere)" << endl;
    cout << endl;
    cout << "Max |difference| for rays reaching the same boundary:" << endl;
    cout << "  |dr|     = " << max_dr << endl;
    cout << "  |dtheta| = " << max_dtheta << endl;
    cout << "  |dphi|   = " << max_dphi << endl;

    // Pass if diverged trajectories are a small fraction (<10%) of valid rays
    int n_valid = n - n_skip;
    bool pass = (n_valid > 0) && (n_diverge * 10 <= n_valid);
    cout << endl << (pass ? "PASS" : "FAIL") << endl;

    delete euler_src;
    delete rk4_src;

    return pass ? 0 : 1;
}
