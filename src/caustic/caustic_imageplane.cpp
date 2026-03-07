//
// caustic_imageplane.cpp
//
// Sky-plane critical curve mapping for Kerr spacetime.
//
// Traces rays backward from a 2D image plane to the accretion disc and computes
// the Jacobian J = d(r_disc, phi_disc)/d(x_img, y_img) at each pixel via central
// finite differences on the neighboring pixels.
//
// Critical curves on the sky are the zero-contours of det(J).  These map to
// caustics on the disc under the lens map.
//
// The rdot_flips count (number of radial turning points during propagation) is
// stored per pixel to label image order: 0 = direct image, 1 = first photon
// ring, 2 = second ring, etc.
//
// FITS output extensions:
//   DET_J      -- det(J) in image-plane coordinates; NaN where undefined
//   SIGN_J     -- sign of det(J): +1/-1; 0 where undefined
//   ORDER      -- rdot_flips (image order); -1 where no disc hit
//   HIT        -- 1 if ray hit disc, 0 otherwise
//   RADIUS     -- disc radius (gravitational radii) of the disc hit
//   PHI        -- disc azimuthal angle (radians, range [-pi,pi])
//   REDSHIFT   -- photon energy ratio E_obs/E_emit (includes gravitational + Doppler)
//
// NOTE on array dimensions:
//   The user specifies Nx, Ny grid steps.  ImagePlane creates (Nx+1)*(Ny+1) rays
//   (fencepost convention).  All output arrays therefore have dimensions
//   img_Nx = Nx+1, img_Ny = Ny+1.
//

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <limits>

using namespace std;

#include "../raytracer/imageplane.h"
#include "../raytracer/ray_destination.h"

#include "../include/fits_output.h"
#include "../include/array.h"
#include "../include/par_file.h"
#include "../include/par_args.h"
#include "../include/kerr.h"

static const double NaN = std::numeric_limits<double>::quiet_NaN();

// Wraps a phi difference into [-pi, pi] to handle the phi = +-pi branch cut.
static inline double wrap_dphi(double dphi)
{
    while (dphi >  M_PI) dphi -= 2*M_PI;
    while (dphi < -M_PI) dphi += 2*M_PI;
    return dphi;
}


int main(int argc, char **argv)
{
    ParameterArgs par_args(argc, argv);

    string par_filename = par_args.key_exists("--parfile")
                          ? par_args.get_string_parameter("--parfile")
                          : "../par/caustic_imageplane.par";

    ParameterFile par_file(par_filename.c_str());
    string out_filename = par_args.key_exists("--outfile")
                          ? par_args.get_parameter<string>("--outfile")
                          : par_file.get_parameter<string>("outfile");

    double dist  = par_file.get_parameter<double>("dist");
    double incl  = par_args.key_exists("--incl")
                   ? par_args.get_parameter<double>("--incl")
                   : par_file.get_parameter<double>("incl");
    double plane_phi0 = par_file.get_parameter<double>("plane_phi0", 0);
    double spin  = par_args.key_exists("--spin")
                   ? par_args.get_parameter<double>("--spin")
                   : par_file.get_parameter<double>("spin");
    double r_disc = par_file.get_parameter<double>("r_disc");

    double x0   = par_file.get_parameter<double>("x0",   -1 * r_disc);
    double xmax = par_file.get_parameter<double>("xmax",  r_disc);
    int    Nx   = par_file.get_parameter<int>("Nx");
    double y0   = par_file.get_parameter<double>("y0",   x0);
    double ymax = par_file.get_parameter<double>("ymax", xmax);
    int    Ny   = par_file.get_parameter<int>("Ny", Nx);

    string integrator_str = par_file.get_parameter<string>("integrator", "rk45");
    double rk45_tol = par_file.get_parameter<double>("rk45_tol", 1e-8);
    double precision = par_file.get_parameter<double>("precision", PRECISION);
    int show_progress = par_args.key_exists("--show_progress")
                        ? par_args.get_parameter<int>("--show_progress")
                        : par_file.get_parameter<int>("show_progress", 1);

    Integrator integrator;
    if (integrator_str == "euler")
    {
        cerr << "Warning: Euler integrator does not support RayDestination; using RK45" << endl;
        integrator = Integrator::RK45;
    }
    else if (integrator_str == "rk4")
        integrator = Integrator::RK4;
    else
        integrator = Integrator::RK45;

    // Step sizes for the ray grid.
    // ImagePlane internally computes Nx_internal = (xmax-x0)/dx + 1 = Nx + 1 (fencepost).
    double dx = (xmax - x0) / Nx;
    double dy = (ymax - y0) / Ny;

    // Actual output array dimensions (one extra point per axis due to fencepost)
    int img_Nx = Nx + 1;
    int img_Ny = Ny + 1;

    double r_isco = kerr_isco<double>(spin, +1);
    cout << "ISCO at r = " << r_isco << endl;
    cout << "Image plane: " << img_Nx << " x " << img_Ny << " = " << img_Nx*img_Ny << " rays" << endl;

    // --- Ray tracing ---

    DiscWithISCODestination<double> disc(r_isco, r_disc);

    ImagePlane<double> raytrace_source(dist, incl, x0, xmax, dx, y0, ymax, dy, spin,
                                       plane_phi0, precision);
    if (integrator == Integrator::RK45)
        raytrace_source.set_rk45_tol(rk45_tol);

    raytrace_source.redshift_start();
    raytrace_source.run_raytrace(&disc, integrator, 1.1 * dist, show_progress);
    raytrace_source.redshift(&disc, true);   // reverse=true for ImagePlane backward trace
    // NOTE: range_phi() is called after data collection so we can read the accumulated
    // phi (not yet normalised) to compute image order = floor(|Δφ| / π).

    // --- Gather per-pixel data into 2D arrays ---
    //
    // get_x_index(ray) returns values in [0, img_Nx-1]
    // get_y_index(ray) returns values in [0, img_Ny-1]

    Array2D<int>    hit_map(img_Nx, img_Ny);
    Array2D<double> r_disc_map(img_Nx, img_Ny);
    Array2D<double> phi_map(img_Nx, img_Ny);
    Array2D<int>    order_map(img_Nx, img_Ny);
    Array2D<double> redshift_map(img_Nx, img_Ny);

    hit_map.zero();
    r_disc_map.zero();
    phi_map.zero();
    order_map.zero();
    redshift_map.zero();

    long disc_count = 0;

    for (int ray = 0; ray < raytrace_source.get_count(); ray++)
    {
        int ix = raytrace_source.get_x_index(ray);
        int iy = raytrace_source.get_y_index(ray);

        // Ray must have terminated successfully (steps > 0) with a valid disc hit:
        // within the disc annulus [r_isco, r_disc] and a positive redshift.
        if (raytrace_source.rays[ray].steps > 0 &&
            raytrace_source.rays[ray].r >= r_isco &&
            raytrace_source.rays[ray].r <  r_disc &&
            raytrace_source.rays[ray].redshift > 0)
        {
            double phi_acc = raytrace_source.rays[ray].phi;  // accumulated (not yet normalised)
            hit_map[ix][iy]      = 1;
            r_disc_map[ix][iy]   = raytrace_source.rays[ray].r;
            phi_map[ix][iy]      = atan2(sin(phi_acc), cos(phi_acc));  // normalise to (-π, π]
            order_map[ix][iy]    = (int)floor(fabs(phi_acc) / M_PI);   // image order from accumulated phi
            redshift_map[ix][iy] = raytrace_source.rays[ray].redshift;
            ++disc_count;
        }
        else
        {
            hit_map[ix][iy]      = 0;
            r_disc_map[ix][iy]   = 0;
            phi_map[ix][iy]      = 0;
            order_map[ix][iy]    = -1;
            redshift_map[ix][iy] = 0;
        }
    }

    cout << disc_count << " rays hit the disc" << endl;

    // Diagnostic: count ray failure modes for non-disc rays
    long cnt_horizon = 0, cnt_rlim = 0, cnt_steplim = 0, cnt_other = 0, cnt_outofrange = 0;
    for (int ray = 0; ray < raytrace_source.get_count(); ray++)
    {
        auto& r = raytrace_source.rays[ray];
        if (r.steps > 0 && (r.r < r_isco || r.r >= r_disc || r.redshift <= 0))
            ++cnt_outofrange;
        else if (r.steps <= 0 || !(r.status & RAY_STATUS_DEST))
        {
            if (r.status & RAY_STATUS_HORIZON)  ++cnt_horizon;
            else if (r.status & RAY_STATUS_RLIM)    ++cnt_rlim;
            else if (r.status & RAY_STATUS_STEPLIM) ++cnt_steplim;
            else ++cnt_other;
        }
    }
    cout << "  -> horizon=" << cnt_horizon << " rlim=" << cnt_rlim
         << " steplim=" << cnt_steplim << " out-of-range=" << cnt_outofrange
         << " other=" << cnt_other << endl;

    // --- Compute Jacobian det(J) via central finite differences ---
    //
    // For pixel (ix,iy), we compute:
    //   dr/dx  = (r[ix+1][iy] - r[ix-1][iy]) / (2*dx)
    //   dr/dy  = (r[ix][iy+1] - r[ix][iy-1]) / (2*dy)
    //   dphi/dx = wrap(phi[ix+1][iy] - phi[ix-1][iy]) / (2*dx)
    //   dphi/dy = wrap(phi[ix][iy+1] - phi[ix][iy-1]) / (2*dy)
    //   det(J) = (dr/dx)(dphi/dy) - (dr/dy)(dphi/dx)
    //
    // The Jacobian is defined only when all four cardinal neighbours hit the disc
    // AND share the same image order as the centre pixel.  Otherwise NaN.
    //
    // At image-order boundaries (neighbours have different image orders), the
    // Jacobian is formally undefined — but these boundaries ARE critical curves
    // (the photon-ring boundaries).  We mark them with a large sentinel value
    // so they are visible in the output.

    const double SENTINEL = 1e30;   // marks order-boundary critical curves

    Array2D<double> det_J_map(img_Nx, img_Ny);
    Array2D<double> sign_J_map(img_Nx, img_Ny);

    for (int ix = 0; ix < img_Nx; ix++)
    {
        for (int iy = 0; iy < img_Ny; iy++)
        {
            det_J_map[ix][iy]  = NaN;
            sign_J_map[ix][iy] = 0;

            // Skip pixels without a disc hit or at the image-plane boundary
            if (!hit_map[ix][iy]) continue;
            if (ix == 0 || ix == img_Nx-1 || iy == 0 || iy == img_Ny-1) continue;

            int  ord     = order_map[ix][iy];
            bool all_hit = hit_map[ix+1][iy] && hit_map[ix-1][iy] &&
                           hit_map[ix][iy+1] && hit_map[ix][iy-1];

            if (!all_hit) continue;

            // Check if any neighbour is in a different image order.
            bool order_match = (order_map[ix+1][iy] == ord) &&
                               (order_map[ix-1][iy] == ord) &&
                               (order_map[ix][iy+1] == ord) &&
                               (order_map[ix][iy-1] == ord);

            if (!order_match)
            {
                // Pixel sits on an image-order boundary = critical curve
                det_J_map[ix][iy] = SENTINEL;
                sign_J_map[ix][iy] = 0;
                continue;
            }

            // Central finite differences
            double dr_dx   = (r_disc_map[ix+1][iy]   - r_disc_map[ix-1][iy])   / (2*dx);
            double dr_dy   = (r_disc_map[ix][iy+1]   - r_disc_map[ix][iy-1])   / (2*dy);
            double dphi_dx = wrap_dphi(phi_map[ix+1][iy] - phi_map[ix-1][iy])  / (2*dx);
            double dphi_dy = wrap_dphi(phi_map[ix][iy+1] - phi_map[ix][iy-1])  / (2*dy);

            double det = dr_dx * dphi_dy - dr_dy * dphi_dx;
            det_J_map[ix][iy]  = det;
            sign_J_map[ix][iy] = (det > 0) ? 1.0 : (det < 0) ? -1.0 : 0.0;
        }
    }

    // --- FITS output ---

    FITSOutput<double> fits(out_filename);
    fits.create_primary();
    fits.write_comment("Kerr spacetime caustic / critical curve mapping (image plane)");
    fits.write_keyword("GENERATOR", "Simulation results were generated by this software", "caustic_imageplane");
    fits.write_keyword("DIST",   "Distance to image plane (rg)", dist);
    fits.write_keyword("INCL",   "Inclination (degrees)", incl);
    fits.write_keyword("SPIN",   "Black hole spin parameter a/M", spin);
    fits.write_keyword("ISCO",   "Innermost stable circular orbit (rg)", r_isco);
    fits.write_keyword("RDISC",  "Outer disc radius (rg)", r_disc);
    fits.write_keyword("NRAYS",  "Total number of rays", img_Nx * img_Ny);
    fits.write_keyword("DISC_N", "Rays that hit the disc", disc_count);

    // Helper lambda to write common axis keywords
    auto write_axis_keywords = [&](const char* axis1_label)
    {
        fits.write_keyword("AXIS1", "Quantity along X axis", axis1_label);
        fits.write_keyword("AXIS2", "Quantity along Y axis", "Image plane Y (rg)");
        fits.write_keyword("X0",   "Start of X axis (rg)", x0);
        fits.write_keyword("XMAX", "End of X axis (rg)", xmax);
        fits.write_keyword("DX",   "X step (rg)", dx);
        fits.write_keyword("NX",   "Number of pixels in X", img_Nx);
        fits.write_keyword("Y0",   "Start of Y axis (rg)", y0);
        fits.write_keyword("YMAX", "End of Y axis (rg)", ymax);
        fits.write_keyword("DY",   "Y step (rg)", dy);
        fits.write_keyword("NY",   "Number of pixels in Y", img_Ny);
    };

    // DET_J: Jacobian determinant.
    fits.write_image(det_J_map, img_Nx, img_Ny, false);
    fits.set_ext_name("DET_J");
    fits.write_comment("Jacobian determinant det(d(r,phi)/d(x,y)); zero-crossings = critical curves");
    write_axis_keywords("Image plane X (rg)");
    fits.write_keyword("SENTINL", "Value used at image-order boundaries (also critical curves)", SENTINEL);

    // SIGN_J: sign of det(J).
    fits.write_image(sign_J_map, img_Nx, img_Ny, false);
    fits.set_ext_name("SIGN_J");
    fits.write_comment("Sign of Jacobian determinant (+1/-1/0); parity flips at critical curves");
    write_axis_keywords("Image plane X (rg)");

    // ORDER: image order (rdot_flips count).
    Array2D<double> order_out(img_Nx, img_Ny);
    for (int ix = 0; ix < img_Nx; ix++)
        for (int iy = 0; iy < img_Ny; iy++)
            order_out[ix][iy] = static_cast<double>(order_map[ix][iy]);
    fits.write_image(order_out, img_Nx, img_Ny, false);
    fits.set_ext_name("ORDER");
    fits.write_comment("Image order (rdot_flips): 0=direct, 1=first photon ring, -1=no hit");
    write_axis_keywords("Image plane X (rg)");

    // HIT: disc-hit mask.
    Array2D<double> hit_out(img_Nx, img_Ny);
    for (int ix = 0; ix < img_Nx; ix++)
        for (int iy = 0; iy < img_Ny; iy++)
            hit_out[ix][iy] = static_cast<double>(hit_map[ix][iy]);
    fits.write_image(hit_out, img_Nx, img_Ny, false);
    fits.set_ext_name("HIT");
    fits.write_comment("1 if ray hit disc, 0 otherwise");
    write_axis_keywords("Image plane X (rg)");

    // RADIUS: disc emission radius.
    fits.write_image(r_disc_map, img_Nx, img_Ny, false);
    fits.set_ext_name("RADIUS");
    fits.write_comment("Disc emission radius (rg); 0 if no disc hit");
    write_axis_keywords("Image plane X (rg)");
    fits.write_keyword("PIXVAL",  "Pixel value quantity", "RADIUS");
    fits.write_keyword("PIXUNIT", "Pixel value unit",     "RG");

    // PHI: disc azimuthal coordinate.
    fits.write_image(phi_map, img_Nx, img_Ny, false);
    fits.set_ext_name("PHI");
    fits.write_comment("Disc azimuthal angle (radians, range [-pi,pi]); 0 if no disc hit");
    write_axis_keywords("Image plane X (rg)");
    fits.write_keyword("PIXVAL",  "Pixel value quantity", "PHI");
    fits.write_keyword("PIXUNIT", "Pixel value unit",     "RAD");

    // REDSHIFT: photon energy ratio E_obs/E_emit.
    fits.write_image(redshift_map, img_Nx, img_Ny, false);
    fits.set_ext_name("REDSHIFT");
    fits.write_comment("Photon energy ratio E_obs/E_emit (gravitational + Doppler); 0 if no disc hit");
    write_axis_keywords("Image plane X (rg)");
    fits.write_keyword("PIXVAL",  "Pixel value quantity", "REDSHIFT");
    fits.write_keyword("PIXUNIT", "Pixel value unit",     "E_OBS/E_EMIT");

    fits.close();

    cout << "Done. Output: " << out_filename << endl;
    return 0;
}
