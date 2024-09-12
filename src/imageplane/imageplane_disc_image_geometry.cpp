//
// Created by drw on 04/04/23.
//

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "raytrace_destination.h"
#include <fitsio.h>
#include <vector>

using namespace std;

#include "../raytracer/imageplane.h"
#include "../raytracer/log_imageplane.h"

#include "../include/fits_output.h"
#include "../include/array.h"
#include "../include/par_file.h"
#include "../include/par_args.h"

double powerlaw3(double r, double q1, double rb1, double q2, double rb2, double q3)
{
    if(r < rb1)
        return pow(r, -1 * q1);
    else if(r < rb2)
        return pow(rb1, q2 - q1) * pow(r, -1 * q2);
    else
        return pow(rb1, q2 - q1) * pow(rb2, q3 - q2) * pow(r, -1 * q3);
}

//int main(int argc, char **argv)
//{
//    ParameterArgs par_args(argc, argv);
//
//    // parameter configuration file
//    char default_par_filename[] = "../par/imageplane_disc_image.par";
//    char *par_filename;
//    if(par_args.key_exists("--parfile"))
//    {
//        string par_filename_str = par_args.get_string_parameter("--parfile");
//        par_filename = (char *) par_filename_str.c_str();
//    }
//    else
//        par_filename = default_par_filename;
//
//    ParameterFile par_file(par_filename);
//    string out_filename = (par_args.key_exists("--outfile")) ? par_args.get_parameter<string>("--outfile")
//                                                             : par_file.get_parameter<string>("outfile");
//    double dist = par_file.get_parameter<double>("dist");
//    double incl = (par_args.key_exists("--incl")) ? par_args.get_parameter<double>("--incl")
//                                                  : par_file.get_parameter<double>("incl");
//
//    double plane_phi0 = (par_args.key_exists("--plane_phi0")) ? par_args.get_parameter<double>("--plane_phi0")
//                                                  : par_file.get_parameter<double>("plane_phi0", 0);
//    double spin = (par_args.key_exists("--spin")) ? par_args.get_parameter<double>("--spin")
//                                                  : par_file.get_parameter<double>("spin");
//    //int num_threads = par_file.get_parameter<int>("num_threads");
//
//    double r_disc = par_file.get_parameter<double>("r_disc");   // radial distance that the disc spans
//    double r_angle_disc_dis = par_file.get_parameter<double>("r_angle_disc_dis");
//    double r_torus = par_file.get_parameter<double>("r_torus");   // radius of the torus
//
//    double thetalim = (par_args.key_exists("--thetalim")) ? par_args.get_parameter<double>("--thetalim")
//            :  par_file.get_parameter<double>("thetalim", 0.001);
//
//    // double theta_lim_1 = par_file.get_parameter<double>("theta_lim_1");   // limiting theta value of the part of the warped disc to the left of the x-z plane at origin
//    // double theta_lim_2 = par_file.get_parameter<double>("theta_lim_2");   // limiting theta value of the part of the warped disc to the right of the x-z plane at origin
//
//    double major_axis = par_file.get_parameter<double>("major_axis");   // major axis size for ellipse geometry
//    double minor_axis = par_file.get_parameter<double>("minor_axis");   // minor axis size for ellipse geometry
//
//    double efficiency = par_file.get_parameter<double>("efficiency");
//    double edd = (par_args.key_exists("--edd")) ? par_args.get_parameter<double>("--edd")
//                                                :  par_file.get_parameter<double>("edd", 0.3); // accretion rate
//    double x0 = par_file.get_parameter<double>("x0", -1 * r_disc);
//    double xmax = par_file.get_parameter<double>("xmax", r_disc);
//    int Nx = par_file.get_parameter<int>("Nx");
//    double y0 = par_file.get_parameter<double>("y0", x0);
//    double ymax = par_file.get_parameter<double>("ymax", xmax);
//    int Ny = par_file.get_parameter<int>("Ny", Nx);
//    int img_Nx = par_file.get_parameter<int>("img_Nx", Nx);
//    int img_Ny = par_file.get_parameter<int>("img_Ny", img_Nx);
//    double q1 = par_file.get_parameter<double>("q1", 3);
//    double rb1 = par_file.get_parameter<double>("rb1", 4);
//    double q2 = par_file.get_parameter<double>("q2", 3);
//    double rb2 = par_file.get_parameter<double>("rb2", 10);
//    double q3 = par_file.get_parameter<double>("q3", 3);
//    double precision = par_file.get_parameter<double>("precision", PRECISION);
//    double max_tstep = par_file.get_parameter<double>("max_tstep", MAXDT);
//    bool flip_image = par_file.get_parameter<bool>("flip_image", true);
//
//    double dx = (xmax - x0) / Nx;
//    double dy = (ymax - y0) / Ny;
//
//    double img_dx = (xmax - x0) / img_Nx;
//    double img_dy = (ymax - y0) / img_Ny;
//
//    long disc_count = 0;
//
//    Array2D<double> disc_flux(img_Nx, img_Ny);
//    Array2D<double> disc_r(img_Nx, img_Ny);
//    Array2D<double> disc_phi(img_Nx, img_Ny);
//    Array2D<double> disc_enshift(img_Nx, img_Ny);
//    Array2D<double> disc_time(img_Nx, img_Ny);
//    Array2D<double> disc_emis(img_Nx, img_Ny);
//    Array2D<int> disc_Nrays(img_Nx, img_Ny);
//    Array2D<double> disc_x(img_Nx, img_Ny);
//    Array2D<double> disc_y(img_Nx, img_Ny);
////    Array2D<double> disc_angle1(img_Nx, img_Ny);
////    Array2D<double> disc_angle2(img_Nx, img_Ny);
//
//    disc_flux.zero();
//    disc_r.zero();
//    disc_phi.zero();
//    disc_enshift.zero();
//    disc_time.zero();
//    disc_emis.zero();
//    disc_Nrays.zero();
//    disc_x.zero();
//    disc_y.zero();
////    disc_angle1.zero();
////    disc_angle2.zero();
//
//    double r_isco = kerr_isco<double>(spin, +1);
//    cout << "ISCO at " << r_isco << endl;
//
////    int threads = 1;
////    omp_set_num_threads(8);
////    #pragma omp parallel
////    threads = omp_get_num_threads();
////    cout << "Threads " << threads << endl;
//    ImagePlane<double> raytrace_source(dist, incl, x0, xmax, dx, y0, ymax, dy, spin, plane_phi0,
//                                                precision);
//    //raytrace_source.set_max_tstep(max_tstep);
//
//    ZDestination<double>* my_destination = new ZDestination<double>(thetalim, r_disc);
//    //AngledDiscsDestination<double>* my_destination = new AngledDiscsDestination<double>(thetalim, thetalim, r_angle_disc_dis);
//    //TorusDiscDestination<double>* my_destination = new TorusDiscDestination<double>(r_torus, r_disc, r_isco);
//    //EllipseDiscDestination<double>* my_destination = new EllipseDiscDestination<double>(r_disc, r_isco, major_axis, minor_axis);
//    //SinDiscDestination<double>* my_destination = new SinDiscDestination<double>(r_disc);
//    //ShakuraDiscDestination<double>* my_destination = new ShakuraDiscDestination<double>(efficiency, edd, r_isco);
//
//
//    raytrace_source.redshift_start();
//    raytrace_source.run_raytrace(my_destination, 1.1 * dist, r_disc);
//    raytrace_source.redshift(my_destination, -1);
//    raytrace_source.range_phi();
//    //raytrace_source.calculate_ray_angles(-1, true);
//
//
//    int num = raytrace_source.get_count();
//    std::vector<double> xValues(num);
//    std::vector<double> yValues(num);
//    std::vector<double> zValues(num);
//
////    for (int ray = 0; ray < raytrace_source.get_count(); ray++) {
////        double x = raytrace_source.rays[ray].alpha;
////        double y = raytrace_source.rays[ray].beta;
////        int ix = static_cast<int>((x - x0) / img_dx);
////        int iy = static_cast<int>((y - y0) / img_dy);
////
////        if (flip_image) {
////            iy = img_Ny - iy - 1;
////        }
////
////        if (ix >= 0 && ix < img_Nx && iy >= 0 && iy < img_Ny)
////            ++disc_Nrays[ix][iy];
////
////        disc_x[ix][iy] += raytrace_source.rays[ray].alpha;
////        disc_y[ix][iy] += raytrace_source.rays[ray].beta;
////    }
//
//    for(int ray = 0; ray < raytrace_source.get_count(); ray++) {
//        if (raytrace_source.rays[ray].steps > 0) {
//            double xx, yy, zz;
//            cartesian(xx, yy, zz, raytrace_source.rays[ray].r, raytrace_source.rays[ray].theta,
//                      raytrace_source.rays[ray].phi, spin);
////            if(z < 1E-2 && raytrace_source.rays[ray].r >= r_isco && raytrace_source.rays[ray].r < r_disc &&
////               raytrace_source.rays[ray].redshift > 0)
//            //  if(raytrace_source.rays[ray].r >= r_isco && raytrace_source.rays[ray].r < r_disc &&
//            //   raytrace_source.rays[ray].redshift > 0)
//            if (raytrace_source.rays[ray].r >= r_isco && raytrace_source.rays[ray].r < r_disc &&
//                raytrace_source.rays[ray].redshift > 0 && raytrace_source.rays[ray].status == RAY_STOP_DEST) {
//                xValues[ray] = xx;  //  saving the positions of the rays into a file
//                yValues[ray] = yy;
//                zValues[ray] = zz;
//
//                double x = raytrace_source.rays[ray].alpha;
//                double y = raytrace_source.rays[ray].beta;
//                int ix = static_cast<int>((x - x0) / img_dx);
//                int iy = static_cast<int>((y - y0) / img_dy);
//
//                if (flip_image) {
//                    iy = img_Ny - iy - 1;
//                }
//
//                if (ix >= 0 && ix < img_Nx && iy >= 0 && iy < img_Ny)
//                    ++disc_Nrays[ix][iy];
//
//                double emis = powerlaw3(raytrace_source.rays[ray].r, q1, rb1, q2, rb2, q3);
//
//                disc_flux[ix][iy] += emis / pow(raytrace_source.rays[ray].redshift, 3);
//                disc_r[ix][iy] += raytrace_source.rays[ray].r;
//                disc_phi[ix][iy] += raytrace_source.rays[ray].phi;
//                disc_enshift[ix][iy] += 1. / raytrace_source.rays[ray].redshift;
//                disc_time[ix][iy] += raytrace_source.rays[ray].t;
//                disc_emis[ix][iy] += emis;
//                //disc_x[ix][iy] += raytrace_source.rays[ray].alpha;
//                //disc_y[ix][iy] += raytrace_source.rays[ray].beta;
//                //disc_angle1[ix][iy] += (180 / M_PI) * raytrace_source.rays[ray].angle1;
//                //disc_angle2[ix][iy] += (180 / M_PI) * raytrace_source.rays[ray].angle2;
//                ++disc_count;
//            }
//        }
//    }
//
//    const char* filename = "output.csv";
//
//    std::ofstream outputFile(filename);
//    if (!outputFile) {
//        std::cout << "Error opening the file." << std::endl;
//        return 1;
//    }
//
//    outputFile << "x, y, z" << std::endl;
//    for (int i = 0; i < num; i++) {
//        outputFile << xValues[i] << ", " << yValues[i] << ", " << zValues[i] << std::endl;
//    }
//
//    outputFile.close();
//
//    cout << disc_count << " rays hit the disc" << endl;
//
//    for(int ix = 0; ix < img_Nx; ix++)
//        for(int iy = 0; iy < img_Ny; iy++)
//            if(disc_Nrays[ix][iy] > 0)
//                disc_flux[ix][iy] /= disc_Nrays[ix][iy];
//
//    disc_r /= disc_Nrays;
//    disc_phi /= disc_Nrays;
//    disc_enshift /= disc_Nrays;
//    disc_time /= disc_Nrays;
//    disc_emis /= disc_Nrays;
////    disc_angle1 /= disc_Nrays;
////    disc_angle2 /= disc_Nrays;
//
//    FITSOutput<double> fits(out_filename);
//    fits.create_primary();
//    fits.write_comment("Raytraced images of accretion disc");
//    fits.write_keyword("GENERATOR", "Simulation results were generated by this software", "imageplane_disc_image");
//    fits.write_keyword("raytrace_source", "Raytrace source configuration", "imageplane");
//    fits.write_keyword("DIST", "Distance to image plane", dist);
//    fits.write_keyword("INCL", "Inclination of line of sight", incl);
//    fits.write_keyword("SPIN", "Black hole spin", spin);
//    fits.write_keyword("ISCO", "Innermost stable circular orbit", r_isco);
//
//    fits.write_keyword("RDISC", "Maximum radius of disc", r_disc);
//    fits.write_keyword("RANGLEDISC", "Distance from center to angled part of disc", r_angle_disc_dis);
//    fits.write_keyword("RTORUS", "Radius of torus", r_torus);
//
//    fits.write_keyword("Q1", "Emissivity profile inner power law index", q1);
//    fits.write_keyword("RB1", "Emissivity profile inner break radius", rb1);
//    fits.write_keyword("Q2", "Emissivity profile middle power law index", q2);
//    fits.write_keyword("RB2", "Emissivity profile outer break radius", rb2);
//    fits.write_keyword("Q3", "Emissivity profile outer power law index", q3);
//    fits.write_keyword("NRAYS", "Number of rays", Nx * Ny);
//    fits.write_keyword("DISCRAYS", "Rays hitting disc", disc_count);
//
//    fits.write_image(disc_flux, img_Nx, img_Ny, false);
//    fits.set_ext_name("FLUX");
//    fits.write_comment("Observed flux");
//    fits.write_keyword("AXIS1", "Quantity along X axis of image", "Image plane X");
//    fits.write_keyword("AXIS2", "Quantity along Y axis of image", "Image plane Y");
//    fits.write_keyword("PIXVAL", "Pixel value quantity", "FLUX");
//    fits.write_keyword("PIXUNIT", "Pixel value unit", "PH*EN^2/DT");
//    fits.write_keyword("X0", "Start of X axis", x0);
//    fits.write_keyword("XMAX", "End of X axis", xmax);
//    fits.write_keyword("DX", "X step", dx);
//    fits.write_keyword("NX", "Number of pixels in X", img_Nx);
//    fits.write_keyword("XUNIT", "Units of X axis", "RG");
//    fits.write_keyword("Y0", "Start of Y axis", y0);
//    fits.write_keyword("YMAX", "End of Y axis", ymax);
//    fits.write_keyword("DY", "Y step", dy);
//    fits.write_keyword("NY", "Number of pixels in Y", img_Ny);
//    fits.write_keyword("YUNIT", "Units of Y axis", "RG");
//
//    fits.write_image(disc_r, img_Nx, img_Ny, false);
//    fits.set_ext_name("RADIUS");
//    fits.write_comment("Emission radius");
//    fits.write_keyword("AXIS1", "Quantity along X axis of image", "Image plane X");
//    fits.write_keyword("AXIS2", "Quantity along Y axis of image", "Image plane Y");
//    fits.write_keyword("PIXVAL", "Pixel value quantity", "RADIUS");
//    fits.write_keyword("PIXUNIT", "Pixel value unit", "RG");
//    fits.write_keyword("X0", "Start of X axis", x0);
//    fits.write_keyword("XMAX", "End of X axis", xmax);
//    fits.write_keyword("DX", "X step", dx);
//    fits.write_keyword("NX", "Number of pixels in X", img_Nx);
//    fits.write_keyword("XUNIT", "Units of X axis", "RG");
//    fits.write_keyword("Y0", "Start of Y axis", y0);
//    fits.write_keyword("YMAX", "End of Y axis", ymax);
//    fits.write_keyword("DY", "Y step", dy);
//    fits.write_keyword("NY", "Number of pixels in Y", img_Ny);
//    fits.write_keyword("YUNIT", "Units of Y axis", "RG");
//
//    fits.write_image(disc_phi, img_Nx, img_Ny, false);
//    fits.set_ext_name("PHI");
//    fits.write_comment("Phi co-ordinate of emission");
//    fits.write_keyword("AXIS1", "Quantity along X axis of image", "Image plane X");
//    fits.write_keyword("AXIS2", "Quantity along Y axis of image", "Image plane Y");
//    fits.write_keyword("PIXVAL", "Pixel value quantity", "PHI");
//    fits.write_keyword("PIXUNIT", "Pixel value unit", "RAD");
//    fits.write_keyword("X0", "Start of X axis", x0);
//    fits.write_keyword("XMAX", "End of X axis", xmax);
//    fits.write_keyword("DX", "X step", dx);
//    fits.write_keyword("NX", "Number of pixels in X", img_Nx);
//    fits.write_keyword("XUNIT", "Units of X axis", "RG");
//    fits.write_keyword("Y0", "Start of Y axis", y0);
//    fits.write_keyword("YMAX", "End of Y axis", ymax);
//    fits.write_keyword("DY", "Y step", dy);
//    fits.write_keyword("NY", "Number of pixels in Y", img_Ny);
//    fits.write_keyword("YUNIT", "Units of Y axis", "RG");
//
//    fits.write_image(disc_enshift, img_Nx, img_Ny, false);
//    fits.set_ext_name("ENSHIFT");
//    fits.write_comment("Observed photon energy");
//    fits.write_keyword("AXIS1", "Quantity along X axis of image", "Image plane X");
//    fits.write_keyword("AXIS2", "Quantity along Y axis of image", "Image plane Y");
//    fits.write_keyword("PIXVAL", "Pixel value quantity", "ENERGY");
//    fits.write_keyword("PIXUNIT", "Pixel value unit", "OBS/EMIT");
//    fits.write_keyword("X0", "Start of X axis", x0);
//    fits.write_keyword("XMAX", "End of X axis", xmax);
//    fits.write_keyword("DX", "X step", dx);
//    fits.write_keyword("NX", "Number of pixels in X", img_Nx);
//    fits.write_keyword("XUNIT", "Units of X axis", "RG");
//    fits.write_keyword("Y0", "Start of Y axis", y0);
//    fits.write_keyword("YMAX", "End of Y axis", ymax);
//    fits.write_keyword("DY", "Y step", dy);
//    fits.write_keyword("NY", "Number of pixels in Y", img_Ny);
//    fits.write_keyword("YUNIT", "Units of Y axis", "RG");
//
//    fits.write_image(disc_time, img_Nx, img_Ny, false);
//    fits.set_ext_name("TIME");
//    fits.write_comment("Ray travel time to image plane");
//    fits.write_keyword("AXIS1", "Quantity along X axis of image", "Image plane X");
//    fits.write_keyword("AXIS2", "Quantity along Y axis of image", "Image plane Y");
//    fits.write_keyword("PIXVAL", "Pixel value quantity", "TIME");
//    fits.write_keyword("PIXUNIT", "Pixel value unit", "GM/c^3");
//    fits.write_keyword("X0", "Start of X axis", x0);
//    fits.write_keyword("XMAX", "End of X axis", xmax);
//    fits.write_keyword("DX", "X step", dx);
//    fits.write_keyword("NX", "Number of pixels in X", img_Nx);
//    fits.write_keyword("XUNIT", "Units of X axis", "RG");
//    fits.write_keyword("Y0", "Start of Y axis", y0);
//    fits.write_keyword("YMAX", "End of Y axis", ymax);
//    fits.write_keyword("DY", "Y step", dy);
//    fits.write_keyword("NY", "Number of pixels in Y", img_Ny);
//    fits.write_keyword("YUNIT", "Units of Y axis", "RG");
//
//    fits.write_image(disc_emis, img_Nx, img_Ny, false);
//    fits.set_ext_name("EMIS");
//    fits.write_comment("Accretion disc emissivity");
//    fits.write_keyword("AXIS1", "Quantity along X axis of image", "Image plane X");
//    fits.write_keyword("AXIS2", "Quantity along Y axis of image", "Image plane Y");
//    fits.write_keyword("PIXVAL", "Pixel value quantity", "EMISSIVITY");
//    fits.write_keyword("PIXUNIT", "Pixel value unit", "PH/RG^2/DT");
//    fits.write_keyword("X0", "Start of X axis", x0);
//    fits.write_keyword("XMAX", "End of X axis", xmax);
//    fits.write_keyword("DX", "X step", dx);
//    fits.write_keyword("NX", "Number of pixels in X", img_Nx);
//    fits.write_keyword("XUNIT", "Units of X axis", "RG");
//    fits.write_keyword("Y0", "Start of Y axis", y0);
//    fits.write_keyword("YMAX", "End of Y axis", ymax);
//    fits.write_keyword("DY", "Y step", dy);
//    fits.write_keyword("NY", "Number of pixels in Y", img_Ny);
//    fits.write_keyword("YUNIT", "Units of Y axis", "RG");
//
//    fits.write_image(disc_x, img_Nx, img_Ny, false);
//    fits.set_ext_name("XVALS");
//    fits.write_comment("Start X-vals");
//    fits.write_keyword("AXIS1", "Quantity along X axis of image", "Image plane X");
//    fits.write_keyword("AXIS2", "Quantity along Y axis of image", "Image plane Y");
//    fits.write_keyword("PIXVAL", "Pixel value quantity", "FLUX");
//    fits.write_keyword("PIXUNIT", "Pixel value unit", "PH*EN^2/DT");
//    fits.write_keyword("X0", "Start of X axis", x0);
//    fits.write_keyword("XMAX", "End of X axis", xmax);
//    fits.write_keyword("DX", "X step", dx);
//    fits.write_keyword("NX", "Number of pixels in X", Nx);
//    fits.write_keyword("XUNIT", "Units of X axis", "RG");
//    fits.write_keyword("Y0", "Start of Y axis", y0);
//    fits.write_keyword("YMAX", "End of Y axis", ymax);
//    fits.write_keyword("DY", "Y step", dy);
//    fits.write_keyword("NY", "Number of pixels in Y", Ny);
//    fits.write_keyword("YUNIT", "Units of Y axis", "RG");
//
//    fits.write_image(disc_y, img_Nx, img_Ny, false);
//    fits.set_ext_name("YVALS");
//    fits.write_comment("Start Y-vals");
//    fits.write_keyword("AXIS1", "Quantity along X axis of image", "Image plane X");
//    fits.write_keyword("AXIS2", "Quantity along Y axis of image", "Image plane Y");
//    fits.write_keyword("PIXVAL", "Pixel value quantity", "FLUX");
//    fits.write_keyword("PIXUNIT", "Pixel value unit", "PH*EN^2/DT");
//    fits.write_keyword("X0", "Start of X axis", x0);
//    fits.write_keyword("XMAX", "End of X axis", xmax);
//    fits.write_keyword("DX", "X step", dx);
//    fits.write_keyword("NX", "Number of pixels in X", Nx);
//    fits.write_keyword("XUNIT", "Units of X axis", "RG");
//    fits.write_keyword("Y0", "Start of Y axis", y0);
//    fits.write_keyword("YMAX", "End of Y axis", ymax);
//    fits.write_keyword("DY", "Y step", dy);
//    fits.write_keyword("NY", "Number of pixels in Y", Ny);
//    fits.write_keyword("YUNIT", "Units of Y axis", "RG");
//
////    fits.write_image(disc_angle1, img_Nx, img_Ny, false);
////    fits.set_ext_name("ANGLE1");
////    fits.write_comment("Angle of ray to disc normal");
////    fits.write_keyword("AXIS1", "Quantity along X axis of image", "Image plane X");
////    fits.write_keyword("AXIS2", "Quantity along Y axis of image", "Image plane Y");
////    fits.write_keyword("PIXVAL", "Pixel value quantity", "ANGLE1");
////    fits.write_keyword("PIXUNIT", "Pixel value unit", "deg");
////    fits.write_keyword("X0", "Start of X axis", x0);
////    fits.write_keyword("XMAX", "End of X axis", xmax);
////    fits.write_keyword("DX", "X step", dx);
////    fits.write_keyword("NX", "Number of pixels in X", img_Nx);
////    fits.write_keyword("XUNIT", "Units of X axis", "RG");
////    fits.write_keyword("Y0", "Start of Y axis", y0);
////    fits.write_keyword("YMAX", "End of Y axis", ymax);
////    fits.write_keyword("DY", "Y step", dy);
////    fits.write_keyword("NY", "Number of pixels in Y", img_Ny);
////    fits.write_keyword("YUNIT", "Units of Y axis", "RG");
////
////    fits.write_image(disc_angle2, img_Nx, img_Ny, false);
////    fits.set_ext_name("ANGLE1");
////    fits.write_comment("Angle of ray to outward");
////    fits.write_keyword("AXIS1", "Quantity along X axis of image", "Image plane X");
////    fits.write_keyword("AXIS2", "Quantity along Y axis of image", "Image plane Y");
////    fits.write_keyword("PIXVAL", "Pixel value quantity", "ANGLE2");
////    fits.write_keyword("PIXUNIT", "Pixel value unit", "deg");
////    fits.write_keyword("X0", "Start of X axis", x0);
////    fits.write_keyword("XMAX", "End of X axis", xmax);
////    fits.write_keyword("DX", "X step", dx);
////    fits.write_keyword("NX", "Number of pixels in X", img_Nx);
////    fits.write_keyword("XUNIT", "Units of X axis", "RG");
////    fits.write_keyword("Y0", "Start of Y axis", y0);
////    fits.write_keyword("YMAX", "End of Y axis", ymax);
////    fits.write_keyword("DY", "Y step", dy);
////    fits.write_keyword("NY", "Number of pixels in Y", img_Ny);
////    fits.write_keyword("YUNIT", "Units of Y axis", "RG");
//
//    fits.close();
//
//    cout << "Done" << endl;
//
//    return 0;
//}


int main(int argc, char **argv)
{
    ParameterArgs par_args(argc, argv);

    // parameter configuration file
    char default_par_filename[] = "../par/imageplane_disc_image.par";
    char *par_filename;
    if(par_args.key_exists("--parfile"))
    {
        string par_filename_str = par_args.get_string_parameter("--parfile");
        par_filename = (char *) par_filename_str.c_str();
    }
    else
        par_filename = default_par_filename;

    ParameterFile par_file(par_filename);
    string out_filename = (par_args.key_exists("--outfile")) ? par_args.get_parameter<string>("--outfile")
                                                             : par_file.get_parameter<string>("outfile");
    double dist = par_file.get_parameter<double>("dist");
    double incl = (par_args.key_exists("--incl")) ? par_args.get_parameter<double>("--incl")
                                                  : par_file.get_parameter<double>("incl");
    double plane_phi0 = (par_args.key_exists("--plane_phi0")) ? par_args.get_parameter<double>("--plane_phi0")
                                                  : par_file.get_parameter<double>("plane_phi0", 0);
    double spin = (par_args.key_exists("--spin")) ? par_args.get_parameter<double>("--spin")
                                                  : par_file.get_parameter<double>("spin");

    double r_disc = (par_args.key_exists("--r_disc")) ? par_args.get_parameter<double>("--r_disc")
                                                          :  par_file.get_parameter<double>("r_disc", 200); // radial distance that the disc spans

    double r_angle_disc_dis = (par_args.key_exists("--r_angle_disc_dis")) ? par_args.get_parameter<double>("--r_angle_disc_dis")
                                                      :  par_file.get_parameter<double>("r_angle_disc_dis", 10);

    double r_torus = par_file.get_parameter<double>("r_torus");   // radius of the torus

    double thetalim = (par_args.key_exists("--thetalim")) ? par_args.get_parameter<double>("--thetalim")
                                                                          :  par_file.get_parameter<double>("thetalim", 0.001); // limiting theta value of the part of the warped disc, looking down the y axis at origin

    double semi_major_axis = par_file.get_parameter<double>("semi_major_axis");   // major axis size for ellipse geometry
    //double semi_minor_axis = par_file.get_parameter<double>("semi_minor_axis");   // minor axis size for ellipse geometry
    double semi_minor_axis = (par_args.key_exists("--semi_minor_axis")) ? par_args.get_parameter<double>("--semi_minor_axis")
                                                :  par_file.get_parameter<double>("semi_minor_axis", 2.5);

    double efficiency = par_file.get_parameter<double>("efficiency");
    double edd = (par_args.key_exists("--edd")) ? par_args.get_parameter<double>("--edd")
                                                          :  par_file.get_parameter<double>("edd", 0.3); // accretion rate

    double x0 = par_file.get_parameter<double>("x0", -1 * r_disc);
    double xmax = par_file.get_parameter<double>("xmax", r_disc);
    int Nx = par_file.get_parameter<int>("Nx");
    double y0 = par_file.get_parameter<double>("y0", x0);
    double ymax = par_file.get_parameter<double>("ymax", xmax);
    int Ny = par_file.get_parameter<int>("Ny", Nx);
    int img_Nx = par_file.get_parameter<int>("img_Nx", Nx);
    int img_Ny = par_file.get_parameter<int>("img_Ny", img_Nx);
    double q1 = par_file.get_parameter<double>("q1", 3);
    double rb1 = par_file.get_parameter<double>("rb1", 4);
    double q2 = par_file.get_parameter<double>("q2", 3);
    double rb2 = par_file.get_parameter<double>("rb2", 10);
    double q3 = par_file.get_parameter<double>("q3", 3);
    double precision = par_file.get_parameter<double>("precision", PRECISION);
    double max_tstep = par_file.get_parameter<double>("max_tstep", MAXDT);
    bool flip_image = par_file.get_parameter<bool>("flip_image", true);
    int geometry = par_file.get_parameter<int>("geometry", 0);

    double dx = exp(log(xmax/x0)/(Nx));
    double dy = exp(log(ymax/y0)/(Ny));

    int ix, iy, full_ix, full_iy;

    // the full grid is twice these dimensions as we have 4 quadrants
    int full_Nx = 2*Nx;
    int full_Ny = 2*Ny;

    cout << "*****" << endl;
    cout << "Image plane at  d = " << dist << " , incl = " << incl << endl;
    cout << "Emissivity profile q1=" << q1 << ", rb1=" << rb1 << ", q2=" << q2 << ", rb2=" << rb2 << ", q3=" << q3 << endl;
    cout << "Disc radius = " << r_disc << endl;
    cout << "Spin a = " << spin << endl;
    cout << "*****" << endl << endl;

    int disc_count = 0;

    Array2D<double> disc_flux(full_Nx, full_Ny);
    Array2D<double> disc_r(full_Nx, full_Ny);
    Array2D<double> disc_phi(full_Nx, full_Ny);
    Array2D<double> disc_enshift(full_Nx, full_Ny);
    Array2D<double> disc_time(full_Nx, full_Ny);
    Array2D<double> disc_emis(full_Nx, full_Ny);
    Array2D<int> disc_Nrays(full_Nx, full_Ny);
    Array2D<double> disc_x(full_Nx, full_Ny);
    Array2D<double> disc_y(full_Nx, full_Ny);
    Array2D<double> disc_weight(full_Nx, full_Ny);

    disc_flux.zero();
    disc_r.zero();
    disc_phi.zero();
    disc_enshift.zero();
    disc_time.zero();
    disc_emis.zero();
    disc_Nrays.zero();
    disc_x.zero();
    disc_y.zero();
    disc_weight.zero();

    double r_isco = kerr_isco<double>(spin, +1);
    cout << "ISCO at " << r_isco << endl;

    LogImagePlane<double> raytrace_source(dist, incl, x0, xmax, dx, y0, ymax, dy, spin, 0);
    int num = raytrace_source.get_count();
    std::vector<double> xValues(4*num);
    std::vector<double> yValues(4*num);
    std::vector<double> zValues(4*num);

    for(int quad=0; quad<4; quad++) {
        int quad_count = 0;

        LogImagePlane<double> raytrace_source(dist, incl, x0, xmax, dx, y0, ymax, dy, spin, quad, plane_phi0);
        //raytrace_source.set_max_tstep(max_tstep);

        //ZDestination<double>* my_destination = new ZDestination<double>(thetalim, r_disc);
        DelayedFlaredDisc<double>* my_destination = new DelayedFlaredDisc<double>(thetalim, r_disc, r_angle_disc_dis);
        //AngledDiscsDestination<double> *my_destination = new AngledDiscsDestination<double>(thetalim, thetalim, r_angle_disc_dis);
        //TorusDiscDestination<double>* my_destination = new TorusDiscDestination<double>(r_torus, r_disc, r_isco);
        //InclPortionDiscDestination<double>* my_destination = new InclPortionDiscDestination<double>(M_PI/4, M_PI/4, r_angle_disc_dis);
        //EllipseDiscDestination<double>* my_destination = new EllipseDiscDestination<double>(r_disc, r_isco, semi_major_axis, semi_minor_axis);
        //SinDiscDestination<double>* my_destination = new SinDiscDestination<double>(r_disc);
        //ShakuraDiscDestination<double>* my_destination = new ShakuraDiscDestination<double>(efficiency, edd, r_isco);


        raytrace_source.redshift_start();
        raytrace_source.run_raytrace(my_destination, 1.1 * dist, r_disc);
        raytrace_source.redshift(my_destination, -1);
        raytrace_source.range_phi();
        //raytrace_source.calculate_ray_angles(-1, true);

        //num = raytrace_source.get_count();
//        std::vector<double> xValues(4*num);
//        std::vector<double> yValues(4*num);
//        std::vector<double> zValues(4*num);

        for (int ray = 0; ray < raytrace_source.get_count(); ray++) {
            ix = ray % (Nx - 1);
            iy = ray / (Nx - 1);



            full_ix = (quad == 0 || quad == 3) ? Nx + ix : Nx - ix;
            full_iy = (quad == 0 || quad == 1) ? Ny + iy : Ny - iy;

            if(full_ix >=0  && full_ix < full_Nx  &&  full_iy >=0  &&  full_iy < full_Ny)
            {
                disc_x[full_ix][full_iy] = raytrace_source.rays[ray].alpha;
                disc_y[full_ix][full_iy] = raytrace_source.rays[ray].beta;
            }
        }

        for (int ray = 0; ray < raytrace_source.get_count(); ray++) {
            if (raytrace_source.rays[ray].steps > 0) {
                double xx, yy, zz;
                cartesian(xx, yy, zz, raytrace_source.rays[ray].r, raytrace_source.rays[ray].theta,
                          raytrace_source.rays[ray].phi, spin);
                //if (raytrace_source.rays[ray].r >= r_angle_disc_dis && raytrace_source.rays[ray].r < r_disc &&
                    //raytrace_source.rays[ray].redshift > 0 && raytrace_source.rays[ray].status == RAY_STOP_DEST) {
                if (raytrace_source.rays[ray].r >= r_isco && raytrace_source.rays[ray].r < r_disc &&
                    raytrace_source.rays[ray].redshift > 0 && raytrace_source.rays[ray].status == RAY_STOP_DEST) {
                //if (raytrace_source.rays[ray].r >= r_isco && raytrace_source.rays[ray].r < r_disc && raytrace_source.rays[ray].status == RAY_STOP_DEST) {
                        xValues[disc_count] = xx;  //  saving the positions of the rays into a file
                        yValues[disc_count] = yy;
                        zValues[disc_count] = zz;

                        double x = raytrace_source.rays[ray].alpha;
                        double y = raytrace_source.rays[ray].beta;

//                    ix = ray % (Nx - 1);
//                    iy = ray / (Nx - 1);

                        ix = ray / (Ny - 1.1);
                        iy = ray % (Ny - 1);

//                    ix = static_cast<int>( (x - x0)/dx ); //+ static_cast<int>(log(x0/x0) / log(dx));
//                    iy = static_cast<int>( (y - y0)/dy ); //+ static_cast<int>(log(y0/y0) / log(dy));

//                    ix = static_cast<int>(log((x - x0)/dx));
//                    iy = static_cast<int>((y - y0)/dy);

//                    cout << "X val image " << ix << " Y val image " << iy << endl;

                        //ix = RaytraceSource->GetXIndex(ray) + static_cast<int>( log(run_x0/x0) / log(dx));
                        //iy = RaytraceSource->GetYIndex(ray) + static_cast<int>( log(run_y0/y0) / log(dy));

                        // the x and y indices within the full grid (remember that Nx and Ny are the midpoints of the grid)
                        full_ix = (quad == 0 || quad == 3) ? Nx + ix : Nx - ix;
                        full_iy = (quad == 0 || quad == 1) ? Ny + iy : Ny - iy;
//                    full_ix = Nx + ix;
//                    full_iy = Ny + iy;

//                    cout << "X val image " << full_ix << " Y val image " << full_iy << endl;

                        double energy = 1 / raytrace_source.rays[ray].redshift;
                        //double emis = powerlaw3(raytrace_source.rays[ray].r, q1, rb1, q2, rb2, q3);

                        if (full_ix >= 0 && full_ix < full_Nx && full_iy >= 0 && full_iy < full_Ny) {
                            disc_emis[full_ix][full_iy] = powerlaw3(raytrace_source.rays[ray].r, q1, rb1, q2, rb2, q3);
                            disc_r[full_ix][full_iy] = raytrace_source.rays[ray].r;
                            disc_phi[full_ix][full_iy] = raytrace_source.rays[ray].phi;
//                        disc_x[full_ix][full_ix] = raytrace_source.rays[ray].alpha;
//                        disc_y[full_ix][full_ix] = raytrace_source.rays[ray].beta;
                            disc_weight[full_ix][full_iy] = raytrace_source.rays[ray].weight;
                            disc_enshift[full_ix][full_iy] = energy;
                            disc_flux[full_ix][full_iy] =
                                    powerlaw3(raytrace_source.rays[ray].r, q1, rb1, q2, rb2, q3) * energy * energy;
                            disc_time[full_ix][full_iy] = raytrace_source.rays[ray].t;
                            ++disc_Nrays[full_ix][full_iy];
                        }

                        ++quad_count;
                        ++disc_count;


//                    double x = raytrace_source.rays[ray].alpha;
//                    double y = raytrace_source.rays[ray].beta;
//                    int ix = static_cast<int>((x - x0) / img_dx);
//                    int iy = static_cast<int>((y - y0) / img_dy);
//
//                    if (flip_image) {
//                        iy = img_Ny - iy - 1;
//                    }
//
//                    if (ix >= 0 && ix < img_Nx && iy >= 0 && iy < img_Ny)
//                        ++disc_Nrays[ix][iy];
//
//                    double emis = powerlaw3(raytrace_source.rays[ray].r, q1, rb1, q2, rb2, q3);
//
//                    disc_flux[ix][iy] += emis / pow(raytrace_source.rays[ray].redshift, 3);
//                    disc_r[ix][iy] += raytrace_source.rays[ray].r;
//                    disc_phi[ix][iy] += raytrace_source.rays[ray].phi;
//                    disc_enshift[ix][iy] += 1. / raytrace_source.rays[ray].redshift;
//                    disc_time[ix][iy] += raytrace_source.rays[ray].t;
//                    disc_emis[ix][iy] += emis;
//                    //disc_angle1[ix][iy] += (180 / M_PI) * raytrace_source.rays[ray].angle1;
//                    //disc_angle2[ix][iy] += (180 / M_PI) * raytrace_source.rays[ray].angle2;
//                    ++disc_count;
                   // }
                }
            }
        }
        cout << quad_count << " rays hit the disk from quadrant " << quad << endl;
    }


    const char* filename = "output.csv";

    std::ofstream outputFile(filename);
    if (!outputFile) {
        std::cout << "Error opening the file." << std::endl;
        return 1;
    }

    outputFile << "x, y, z" << std::endl;
    for (int i = 0; i < 4*num; i++) {
        outputFile << xValues[i] << ", " << yValues[i] << ", " << zValues[i] << std::endl;
    }

    outputFile.close();

    cout << disc_count << " rays hit the disc" << endl;

//    for(int ix = 0; ix < full_Nx; ix++)
//        for(int iy = 0; iy < full_Ny; iy++)
//            if(disc_Nrays[ix][iy] > 0)
//                disc_flux[ix][iy] /= disc_Nrays[ix][iy];
//
//    disc_r /= disc_Nrays;
//    disc_phi /= disc_Nrays;
//    disc_enshift /= disc_Nrays;
//    disc_time /= disc_Nrays;
//    disc_emis /= disc_Nrays;
//    disc_angle1 /= disc_Nrays;
//    disc_angle2 /= disc_Nrays;

    FITSOutput<double> fits(out_filename);
    fits.create_primary();
    fits.write_comment("Raytraced images of accretion disc");
    fits.write_keyword("GENERATOR", "Simulation results were generated by this software", "imageplane_disc_image");
    fits.write_keyword("raytrace_source", "Raytrace source configuration", "imageplane");
    fits.write_keyword("DIST", "Distance to image plane", dist);
    fits.write_keyword("INCL", "Inclination of line of sight", incl);
    fits.write_keyword("SPIN", "Black hole spin", spin);
    fits.write_keyword("ISCO", "Innermost stable circular orbit", r_isco);

    fits.write_keyword("RDISC", "Maximum radius of disc", r_disc);
    fits.write_keyword("RANGLEDISC", "Distance from center to angled part of disc", r_angle_disc_dis);
    fits.write_keyword("RTORUS", "Radius of torus", r_torus);

    fits.write_keyword("Q1", "Emissivity profile inner power law index", q1);
    fits.write_keyword("RB1", "Emissivity profile inner break radius", rb1);
    fits.write_keyword("Q2", "Emissivity profile middle power law index", q2);
    fits.write_keyword("RB2", "Emissivity profile outer break radius", rb2);
    fits.write_keyword("Q3", "Emissivity profile outer power law index", q3);
    fits.write_keyword("NRAYS", "Number of rays", Nx * Ny);
    fits.write_keyword("DISCRAYS", "Rays hitting disc", disc_count);

    fits.write_image(disc_flux, full_Nx, full_Ny, false);
    fits.set_ext_name("FLUX");
    fits.write_comment("Observed flux");
    fits.write_keyword("AXIS1", "Quantity along X axis of image", "Image plane X");
    fits.write_keyword("AXIS2", "Quantity along Y axis of image", "Image plane Y");
    fits.write_keyword("PIXVAL", "Pixel value quantity", "FLUX");
    fits.write_keyword("PIXUNIT", "Pixel value unit", "PH*EN^2/DT");
    fits.write_keyword("X0", "Start of X axis", x0);
    fits.write_keyword("XMAX", "End of X axis", xmax);
    fits.write_keyword("DX", "X step", dx);
    fits.write_keyword("NX", "Number of pixels in X", Nx);
    fits.write_keyword("XUNIT", "Units of X axis", "RG");
    fits.write_keyword("Y0", "Start of Y axis", y0);
    fits.write_keyword("YMAX", "End of Y axis", ymax);
    fits.write_keyword("DY", "Y step", dy);
    fits.write_keyword("NY", "Number of pixels in Y", Ny);
    fits.write_keyword("YUNIT", "Units of Y axis", "RG");

    fits.write_image(disc_r, full_Nx, full_Ny, false);
    fits.set_ext_name("RADIUS");
    fits.write_comment("Emission radius");
    fits.write_keyword("AXIS1", "Quantity along X axis of image", "Image plane X");
    fits.write_keyword("AXIS2", "Quantity along Y axis of image", "Image plane Y");
    fits.write_keyword("PIXVAL", "Pixel value quantity", "RADIUS");
    fits.write_keyword("PIXUNIT", "Pixel value unit", "RG");
    fits.write_keyword("X0", "Start of X axis", x0);
    fits.write_keyword("XMAX", "End of X axis", xmax);
    fits.write_keyword("DX", "X step", dx);
    fits.write_keyword("NX", "Number of pixels in X", Nx);
    fits.write_keyword("XUNIT", "Units of X axis", "RG");
    fits.write_keyword("Y0", "Start of Y axis", y0);
    fits.write_keyword("YMAX", "End of Y axis", ymax);
    fits.write_keyword("DY", "Y step", dy);
    fits.write_keyword("NY", "Number of pixels in Y", Ny);
    fits.write_keyword("YUNIT", "Units of Y axis", "RG");

    fits.write_image(disc_phi, full_Nx, full_Ny, false);
    fits.set_ext_name("PHI");
    fits.write_comment("Phi co-ordinate of emission");
    fits.write_keyword("AXIS1", "Quantity along X axis of image", "Image plane X");
    fits.write_keyword("AXIS2", "Quantity along Y axis of image", "Image plane Y");
    fits.write_keyword("PIXVAL", "Pixel value quantity", "PHI");
    fits.write_keyword("PIXUNIT", "Pixel value unit", "RAD");
    fits.write_keyword("X0", "Start of X axis", x0);
    fits.write_keyword("XMAX", "End of X axis", xmax);
    fits.write_keyword("DX", "X step", dx);
    fits.write_keyword("NX", "Number of pixels in X", Nx);
    fits.write_keyword("XUNIT", "Units of X axis", "RG");
    fits.write_keyword("Y0", "Start of Y axis", y0);
    fits.write_keyword("YMAX", "End of Y axis", ymax);
    fits.write_keyword("DY", "Y step", dy);
    fits.write_keyword("NY", "Number of pixels in Y", Ny);
    fits.write_keyword("YUNIT", "Units of Y axis", "RG");

    fits.write_image(disc_enshift, full_Nx, full_Ny, false);
    fits.set_ext_name("ENSHIFT");
    fits.write_comment("Observed photon energy");
    fits.write_keyword("AXIS1", "Quantity along X axis of image", "Image plane X");
    fits.write_keyword("AXIS2", "Quantity along Y axis of image", "Image plane Y");
    fits.write_keyword("PIXVAL", "Pixel value quantity", "ENERGY");
    fits.write_keyword("PIXUNIT", "Pixel value unit", "OBS/EMIT");
    fits.write_keyword("X0", "Start of X axis", x0);
    fits.write_keyword("XMAX", "End of X axis", xmax);
    fits.write_keyword("DX", "X step", dx);
    fits.write_keyword("NX", "Number of pixels in X", Nx);
    fits.write_keyword("XUNIT", "Units of X axis", "RG");
    fits.write_keyword("Y0", "Start of Y axis", y0);
    fits.write_keyword("YMAX", "End of Y axis", ymax);
    fits.write_keyword("DY", "Y step", dy);
    fits.write_keyword("NY", "Number of pixels in Y", Ny);
    fits.write_keyword("YUNIT", "Units of Y axis", "RG");

    fits.write_image(disc_time, full_Nx, full_Ny, false);
    fits.set_ext_name("TIME");
    fits.write_comment("Ray travel time to image plane");
    fits.write_keyword("AXIS1", "Quantity along X axis of image", "Image plane X");
    fits.write_keyword("AXIS2", "Quantity along Y axis of image", "Image plane Y");
    fits.write_keyword("PIXVAL", "Pixel value quantity", "TIME");
    fits.write_keyword("PIXUNIT", "Pixel value unit", "GM/c^3");
    fits.write_keyword("X0", "Start of X axis", x0);
    fits.write_keyword("XMAX", "End of X axis", xmax);
    fits.write_keyword("DX", "X step", dx);
    fits.write_keyword("NX", "Number of pixels in X", Nx);
    fits.write_keyword("XUNIT", "Units of X axis", "RG");
    fits.write_keyword("Y0", "Start of Y axis", y0);
    fits.write_keyword("YMAX", "End of Y axis", ymax);
    fits.write_keyword("DY", "Y step", dy);
    fits.write_keyword("NY", "Number of pixels in Y", Ny);
    fits.write_keyword("YUNIT", "Units of Y axis", "RG");

    fits.write_image(disc_emis, full_Nx, full_Ny, false);
    fits.set_ext_name("EMIS");
    fits.write_comment("Accretion disc emissivity");
    fits.write_keyword("AXIS1", "Quantity along X axis of image", "Image plane X");
    fits.write_keyword("AXIS2", "Quantity along Y axis of image", "Image plane Y");
    fits.write_keyword("PIXVAL", "Pixel value quantity", "EMISSIVITY");
    fits.write_keyword("PIXUNIT", "Pixel value unit", "PH/RG^2/DT");
    fits.write_keyword("X0", "Start of X axis", x0);
    fits.write_keyword("XMAX", "End of X axis", xmax);
    fits.write_keyword("DX", "X step", dx);
    fits.write_keyword("NX", "Number of pixels in X", Nx);
    fits.write_keyword("XUNIT", "Units of X axis", "RG");
    fits.write_keyword("Y0", "Start of Y axis", y0);
    fits.write_keyword("YMAX", "End of Y axis", ymax);
    fits.write_keyword("DY", "Y step", dy);
    fits.write_keyword("NY", "Number of pixels in Y", Ny);
    fits.write_keyword("YUNIT", "Units of Y axis", "RG");

    fits.write_image(disc_x, full_Nx, full_Ny, false);
    fits.set_ext_name("XVALS");
    fits.write_comment("Start X-vals");
    fits.write_keyword("AXIS1", "Quantity along X axis of image", "Image plane X");
    fits.write_keyword("AXIS2", "Quantity along Y axis of image", "Image plane Y");
    fits.write_keyword("PIXVAL", "Pixel value quantity", "FLUX");
    fits.write_keyword("PIXUNIT", "Pixel value unit", "PH*EN^2/DT");
    fits.write_keyword("X0", "Start of X axis", x0);
    fits.write_keyword("XMAX", "End of X axis", xmax);
    fits.write_keyword("DX", "X step", dx);
    fits.write_keyword("NX", "Number of pixels in X", Nx);
    fits.write_keyword("XUNIT", "Units of X axis", "RG");
    fits.write_keyword("Y0", "Start of Y axis", y0);
    fits.write_keyword("YMAX", "End of Y axis", ymax);
    fits.write_keyword("DY", "Y step", dy);
    fits.write_keyword("NY", "Number of pixels in Y", Ny);
    fits.write_keyword("YUNIT", "Units of Y axis", "RG");

    fits.write_image(disc_y, full_Nx, full_Ny, false);
    fits.set_ext_name("YVALS");
    fits.write_comment("Start Y-vals");
    fits.write_keyword("AXIS1", "Quantity along X axis of image", "Image plane X");
    fits.write_keyword("AXIS2", "Quantity along Y axis of image", "Image plane Y");
    fits.write_keyword("PIXVAL", "Pixel value quantity", "FLUX");
    fits.write_keyword("PIXUNIT", "Pixel value unit", "PH*EN^2/DT");
    fits.write_keyword("X0", "Start of X axis", x0);
    fits.write_keyword("XMAX", "End of X axis", xmax);
    fits.write_keyword("DX", "X step", dx);
    fits.write_keyword("NX", "Number of pixels in X", Nx);
    fits.write_keyword("XUNIT", "Units of X axis", "RG");
    fits.write_keyword("Y0", "Start of Y axis", y0);
    fits.write_keyword("YMAX", "End of Y axis", ymax);
    fits.write_keyword("DY", "Y step", dy);
    fits.write_keyword("NY", "Number of pixels in Y", Ny);
    fits.write_keyword("YUNIT", "Units of Y axis", "RG");

    fits.write_image(disc_weight, full_Nx, full_Ny, false);
    fits.set_ext_name("WEIGHT");
    fits.write_comment("weights");
    fits.write_keyword("AXIS1", "Quantity along X axis of image", "Image plane X");
    fits.write_keyword("AXIS2", "Quantity along Y axis of image", "Image plane Y");
    fits.write_keyword("PIXVAL", "Pixel value quantity", "FLUX");
    fits.write_keyword("PIXUNIT", "Pixel value unit", "PH*EN^2/DT");
    fits.write_keyword("X0", "Start of X axis", x0);
    fits.write_keyword("XMAX", "End of X axis", xmax);
    fits.write_keyword("DX", "X step", dx);
    fits.write_keyword("NX", "Number of pixels in X", Nx);
    fits.write_keyword("XUNIT", "Units of X axis", "RG");
    fits.write_keyword("Y0", "Start of Y axis", y0);
    fits.write_keyword("YMAX", "End of Y axis", ymax);
    fits.write_keyword("DY", "Y step", dy);
    fits.write_keyword("NY", "Number of pixels in Y", Ny);
    fits.write_keyword("YUNIT", "Units of Y axis", "RG");

//    fits.write_image(disc_angle1, img_Nx, img_Ny, false);
//    fits.set_ext_name("ANGLE1");
//    fits.write_comment("Angle of ray to disc normal");
//    fits.write_keyword("AXIS1", "Quantity along X axis of image", "Image plane X");
//    fits.write_keyword("AXIS2", "Quantity along Y axis of image", "Image plane Y");
//    fits.write_keyword("PIXVAL", "Pixel value quantity", "ANGLE1");
//    fits.write_keyword("PIXUNIT", "Pixel value unit", "deg");
//    fits.write_keyword("X0", "Start of X axis", x0);
//    fits.write_keyword("XMAX", "End of X axis", xmax);
//    fits.write_keyword("DX", "X step", dx);
//    fits.write_keyword("NX", "Number of pixels in X", img_Nx);
//    fits.write_keyword("XUNIT", "Units of X axis", "RG");
//    fits.write_keyword("Y0", "Start of Y axis", y0);
//    fits.write_keyword("YMAX", "End of Y axis", ymax);
//    fits.write_keyword("DY", "Y step", dy);
//    fits.write_keyword("NY", "Number of pixels in Y", img_Ny);
//    fits.write_keyword("YUNIT", "Units of Y axis", "RG");
//
//    fits.write_image(disc_angle2, img_Nx, img_Ny, false);
//    fits.set_ext_name("ANGLE1");
//    fits.write_comment("Angle of ray to outward");
//    fits.write_keyword("AXIS1", "Quantity along X axis of image", "Image plane X");
//    fits.write_keyword("AXIS2", "Quantity along Y axis of image", "Image plane Y");
//    fits.write_keyword("PIXVAL", "Pixel value quantity", "ANGLE2");
//    fits.write_keyword("PIXUNIT", "Pixel value unit", "deg");
//    fits.write_keyword("X0", "Start of X axis", x0);
//    fits.write_keyword("XMAX", "End of X axis", xmax);
//    fits.write_keyword("DX", "X step", dx);
//    fits.write_keyword("NX", "Number of pixels in X", img_Nx);
//    fits.write_keyword("XUNIT", "Units of X axis", "RG");
//    fits.write_keyword("Y0", "Start of Y axis", y0);
//    fits.write_keyword("YMAX", "End of Y axis", ymax);
//    fits.write_keyword("DY", "Y step", dy);
//    fits.write_keyword("NY", "Number of pixels in Y", img_Ny);
//    fits.write_keyword("YUNIT", "Units of Y axis", "RG");

    fits.close();

    cout << "Done" << endl;

    return 0;
}