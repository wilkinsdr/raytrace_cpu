/*
 * gpu_raytracer.h
 *
 *  Created on: 4 Sep 2013
 *      Author: drw
 *
 *  General relativistic ray tracing in the Kerr spacetime
 *
 *  This class provides the basic functionality fundamental to the rays (e.g. propagation, redshift calculation)
 *  and is used as the base class for specific types of source for ray tracing in different scenarios.
 *
 *  This class is a template to set the floating point precision (e.g. float or double)
 */

#ifndef MAPPER_H_
#define MAPPER_H_

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

#include "../include/kerr.h"
#include "../include/array.h"
#include "../raytracer/raytracer.h"
#include "../include/progress_bar.h"

#include "H5Cpp.h"
using namespace H5;

class Mapper : public Raytracer<double>
{
private:
	int motion;
	int vel_mode;
	double vel;

	bool reverse;

    virtual long get_num_rays()
    {
        cerr << "I'm supposed to be overriden!" << endl;
        return 0;
    }

public:
    Mapper(int num_rays, float spin_par, double init_r0, double init_rmax, int init_Nr, int init_Ntheta, int init_Nphi, bool init_logbin_r, double init_thetamax = M_PI_2, float toler = TOL, bool reverse = false);
    Mapper(char* load_filename);
    ~Mapper( );

	double r0, rmax, theta_max;
	int Nr, Ntheta, Nphi;
	double bin_dr, bin_dtheta, bin_dphi;
	bool logbin_r;

    long num_rays;

	Array3D<double> *map_time, *map_redshift, *map_flux;
	Array3D<int> *map_Nrays;
	Array3D<double> *bin_volume;

    void run_map( double r_max = 1000, int show_progress = 1 );
    inline int map_ray(int ray, const double rlim, const double thetalim, const int steplim);
    void average_rays();

	void set_motion(double init_vel, int init_motion = 1,  int init_vel_mode = 0)
	{
		vel = init_vel;
		vel_mode = init_vel_mode;
		motion = init_motion;
	}

    void save(char* filename);

    void save_hdf(char* filename)
    {
	    H5File h5outfile(filename, H5F_ACC_TRUNC);

		DataSpace attr_dataspace = DataSpace(H5S_SCALAR);

		Attribute attr_r0 = h5outfile.createAttribute("r0", PredType::NATIVE_DOUBLE, attr_dataspace);
		attr_r0.write(PredType::NATIVE_DOUBLE, &r0);
		Attribute attr_rmax = h5outfile.createAttribute("rmax", PredType::NATIVE_DOUBLE, attr_dataspace);
		attr_rmax.write(PredType::NATIVE_DOUBLE, &rmax);
		Attribute attr_Nr = h5outfile.createAttribute("Nr", PredType::NATIVE_INT, attr_dataspace);
		attr_Nr.write(PredType::NATIVE_INT, &Nr);
		Attribute attr_dr = h5outfile.createAttribute("dr", PredType::NATIVE_DOUBLE, attr_dataspace);
		attr_dr.write(PredType::NATIVE_DOUBLE, &bin_dr);
		Attribute attr_logbin_r = h5outfile.createAttribute("logbin_r", PredType::NATIVE_INT, attr_dataspace);
		int logbin_r_int = (logbin_r) ? 1 : 0;
		attr_logbin_r.write(PredType::NATIVE_INT, &logbin_r_int);
		Attribute attr_theta_max = h5outfile.createAttribute("theta_max", PredType::NATIVE_DOUBLE, attr_dataspace);
		attr_theta_max.write(PredType::NATIVE_DOUBLE, &theta_max);
		Attribute attr_Ntheta = h5outfile.createAttribute("Ntheta", PredType::NATIVE_INT, attr_dataspace);
		attr_Ntheta.write(PredType::NATIVE_INT, &Ntheta);
		Attribute attr_dtheta = h5outfile.createAttribute("dtheta", PredType::NATIVE_DOUBLE, attr_dataspace);
		attr_dtheta.write(PredType::NATIVE_DOUBLE, &bin_dtheta);
		Attribute attr_Nphi = h5outfile.createAttribute("Nphi", PredType::NATIVE_INT, attr_dataspace);
		attr_Nphi.write(PredType::NATIVE_INT, &Nphi);
		Attribute attr_dphi = h5outfile.createAttribute("dphi", PredType::NATIVE_DOUBLE, attr_dataspace);
		attr_dphi.write(PredType::NATIVE_DOUBLE, &bin_dphi);

		map_time->write_hdf(&h5outfile, "time");
		map_redshift->write_hdf(&h5outfile, "redshift");
		map_Nrays->write_hdf(&h5outfile, "Nrays");
		bin_volume->write_hdf(&h5outfile, "volume");
    }

    void calculate_volume();

    inline int r_index(double r)
    {
        int ir = (logbin_r) ? static_cast<int>( log(r / r0) / log(bin_dr)) : static_cast<int>((r - r0) / bin_dr);

        if(ir < 0) ir = 0;
        else if  (ir >= Nr) ir = Nr - 1;

        return ir;
    }
    inline int theta_index(double theta)
    {
        int itheta = static_cast<int>(theta / bin_dtheta);

        if(itheta < 0) itheta = 0;
        else if  (itheta >= Ntheta) itheta = Ntheta - 1;

        return itheta;
    }
    inline int phi_index(double phi)
    {
        int iphi = static_cast<int>((phi + M_PI)/ bin_dphi);

        if(iphi < 0) iphi = 0;
        else if  (iphi >= Nphi) iphi = Nphi - 1;

        return iphi;
    }

    inline double bin_r(int ir)
    {
        return (logbin_r) ? r0 * pow(bin_dr, ir) : r0 + ir*bin_dr;
    }
    inline double bin_theta(int itheta)
    {
        return itheta*bin_dtheta;
    }
    inline double bin_phi(int iphi)
    {
        return iphi*bin_dphi - M_PI;
    }

};

#endif /* MAPPER_H_ */
