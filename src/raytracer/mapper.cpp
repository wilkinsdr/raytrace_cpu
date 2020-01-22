/*
 * raytracer.cu
 *
 *  Created on: 4 Sep 2013
 *      Author: drw
 */

#include "mapper.h"


Mapper::Mapper( int num_rays, float spin_par, double init_r0, double init_rmax, int init_Nr, int init_Ntheta, int init_Nphi, bool init_logbin_r, double init_thetamax, float toler, bool init_reverse )
	: Raytracer<double>(num_rays, spin_par, toler), r0(init_r0), rmax(init_rmax), Nr(init_Nr), Ntheta(init_Ntheta), Nphi(init_Nphi), logbin_r(init_logbin_r), theta_max(init_thetamax), reverse(init_reverse)
{
	bin_dr = (logbin_r) ? exp(log(rmax/r0)/(Nr-1)) : (rmax - r0)/(Nr - 1);
	bin_dtheta = (M_PI_2)/(Ntheta - 1);
	bin_dphi = (2*M_PI)/(Nphi - 1);

	map_time = new Array3D<double>(Nr, Ntheta, Nphi);
	map_redshift = new Array3D<double>(Nr, Ntheta, Nphi);
	map_Nrays = new Array3D<int>(Nr, Ntheta, Nphi);

	map_time->zero();
	map_redshift->zero();
	map_Nrays->zero();

	bin_volume = new Array3D<double>(Nr, Ntheta, Nphi);
	calculate_volume();
}


Mapper::Mapper(char* load_filename) : Raytracer<double>(1, 0)
{
	ifstream infile(load_filename, ios::binary);
	infile.read(reinterpret_cast<char *> (&r0), sizeof(double));
	infile.read(reinterpret_cast<char *> (&rmax), sizeof(double));
	infile.read(reinterpret_cast<char *> (&Nr), sizeof(int));
	infile.read(reinterpret_cast<char *> (&logbin_r), sizeof(bool));
	infile.read(reinterpret_cast<char *> (&theta_max), sizeof(double));
	infile.read(reinterpret_cast<char *> (&Ntheta), sizeof(int));
	infile.read(reinterpret_cast<char *> (&Nphi), sizeof(int));

	bin_dr = (logbin_r) ? exp(log(rmax/r0)/(Nr-1)) : (rmax - r0)/(Nr - 1);
	bin_dtheta = (M_PI_2)/(Ntheta - 1);
	bin_dphi = (2*M_PI)/(Nphi - 1);

	map_time = new Array3D<double>(Nr, Ntheta, Nphi);
	map_redshift = new Array3D<double>(Nr, Ntheta, Nphi);
	map_Nrays = new Array3D<int>(Nr, Ntheta, Nphi);

	map_time->read(&infile);
	map_redshift->read(&infile);
	map_Nrays->read(&infile);
	bin_volume->read(&infile);

	infile.close();
}


Mapper::~Mapper( )
{

}


void Mapper::run_map( double r_max, int show_progress )
{
	//
	// Runs the ray tracing algorithm once the rays have been set up.
	// Integration of geodesic equation for each ray proceeds until it reaches limiting radius or is lost through the event horizon.
	// Alternatively, Integration stop when the maximum number of steps (defined by STEPLIM in the header file) is reached to prevent
	// incredibly long loops of very small steps.
	//
	// To work around thread time limits on devices also running an X server, each kernel execution
	// only integrates for the number of steps defined in THREAD_STEPLIM. An unfinished flag is then
	// set and the kernel is executed repeatedly until each ray has finished.
	//
	// Calls the GPURaytrace kernel to perform calculation on the GPU.
	//
	// Arguments:
	//	r_max		double		Limiting outer radius for propagation in Rg (default value = 1000)
	//
	cout << "Running mapper..." << endl;

	ProgressBar prog(Raytracer<double>::nRays, "Ray", 0, (show_progress > 0));
	for(int ray=0; ray<Raytracer<double>::nRays; ray++)
	{
		//if(ray % 1 == 0) cout << "\rRay " << ray << '/' << Raytracer<double>::nRays;
		if(show_progress != 0 && (ray % show_progress) == 0) prog.show(ray+1);
		if(Raytracer<double>::m_steps[ray] == -1) continue;
		else if(Raytracer<double>::m_steps[ray] >= STEPLIM) continue;

		int n;
		n = map_ray(ray, r_max, theta_max, STEPLIM);
		Raytracer<double>::m_steps[ray] += n;
	}
	//prog.show_complete();
	prog.done();
	//cout << endl;

	average_rays();
}


inline int Mapper::map_ray(int ray, const double rlim, const double thetalim, const int steplim)
{
	//
	// propagate the photon along its geodesic until limiting r or theta reached
	//
	int steps = 0;

	int rsign_count = COUNT_MIN;
	int thetasign_count = COUNT_MIN;
	
	double rdotsq, thetadotsq;

	double step;

	double x, y, z;

	// copy variables locally to simplify the code
	double a = Raytracer<double>::spin;
	double t = Raytracer<double>::m_t[ray];
	double r = Raytracer<double>::m_r[ray];
	double theta = Raytracer<double>::m_theta[ray];
	double phi = Raytracer<double>::m_phi[ray];
	double pt = Raytracer<double>::m_pt[ray];
	double pr = Raytracer<double>::m_pr[ray];
	double ptheta = Raytracer<double>::m_ptheta[ray];
	double pphi = Raytracer<double>::m_pphi[ray];
	int rdot_sign = Raytracer<double>::m_rdot_sign[ray];
	int thetadot_sign = Raytracer<double>::m_thetadot_sign[ray];

	const double k = Raytracer<double>::m_k[ray];
	const double h = Raytracer<double>::m_h[ray];
	const double Q = Raytracer<double>::m_Q[ray];

	bool write_started = false;

	int ir, itheta, iphi;
	double V, redshift;

	int last_ir = -1;
	int last_itheta = -1;
	int last_iphi = -1;

	// integrate geodesic equations until limit reached
	// if thetalim is positive, we go until theta exceeds it, if it is negative, we go until it is less than the abs value to allow tracing back to theta=0
	while( r < rlim  && ( (thetalim > 0 && theta < thetalim) || (thetalim < 0 && theta > abs(thetalim)) || thetalim == 0 )  &&  steps < steplim )
	//while( r < rlim  && theta < thetalim  &&  steps < steplim )
	{
		++steps;

		const double rhosq = r*r + (a*cos(theta))*(a*cos(theta));
		const double delta = r*r - 2*r + a*a;
		const double sigmasq = (r*r + a*a)*(r*r + a*a) - a*a*delta*sin(theta)*sin(theta);
		const double e2nu = rhosq * delta / sigmasq;
		const double e2psi = sigmasq * sin(theta)*sin(theta) / rhosq;
		const double omega = 2*a*r / sigmasq;
		const double grr = -rhosq/delta;
		const double gthth = -rhosq;
		const double gphph = -e2psi;

		// tdot
		pt = (rhosq*(r*r + a*a) + 2*a*a*r*sin(theta)*sin(theta))*k - 2*a*r*h;
		pt = pt / ( r*r * (1 + (a*cos(theta)/r)*(a*cos(theta)/r) - 2/r)*(r*r + a*a) + 2*a*a*r*sin(theta)*sin(theta) );

		// phidot
		pphi = 2*a*r*sin(theta)*sin(theta)*k + (r*r + (a*cos(theta))*(a*cos(theta)) - 2*r)*h;
		pphi = pphi / ( (r*r + a*a)*(r*r + (a*cos(theta))*(a*cos(theta)) - 2*r)*sin(theta)*sin(theta) + 2*a*a*r*sin(theta)*sin(theta)*sin(theta)*sin(theta) );

		// thetadot
		thetadotsq = Q + (k*a*cos(theta) + h/tan(theta))*(k*a*cos(theta) - h/tan(theta));
		thetadotsq = thetadotsq / (rhosq*rhosq);

		if(thetadotsq < 0 && thetasign_count >= COUNT_MIN)
		{
			thetadot_sign *= -1;
			thetasign_count = 0;
			continue;
		}
		if (thetasign_count <= COUNT_MIN) thetasign_count++;

		// take the square roots and get the right signs
		ptheta = sqrt(abs(thetadotsq)) * thetadot_sign;

		// rdot
		rdotsq = k*pt - h*pphi - rhosq*ptheta*ptheta;
		rdotsq = rdotsq * delta/rhosq;

		if(rdotsq < 0 && rsign_count >= COUNT_MIN)
		{
			rdot_sign *= -1;
			rsign_count = 0;
			continue;
		}
		if (rsign_count <= COUNT_MIN) rsign_count++;

		pr = sqrt(abs(rdotsq)) * rdot_sign;

		step = abs( (r-(double)Raytracer<double>::horizon)/pr ) / Raytracer<double>::tolerance;
		// if the step is smaller in theta or phi (near 0/pi/2pi), use that instead
		if( step > abs( theta/ptheta ) / Raytracer<double>::tolerance ) step = abs( theta/ptheta ) / Raytracer<double>::tolerance;
		if( step > abs( phi/pphi ) / Raytracer<double>::tolerance ) step = abs( phi/pphi ) / Raytracer<double>::tolerance;
//		if( step > abs( (phi - M_PI)/pphi ) / tol ) step = abs( (phi - M_PI)/pphi ) / tol;
//		if( step > abs( (phi - 2*M_PI)/pphi ) / tol ) step = abs( (phi - 2*M_PI)/pphi ) / tol;
		// don't let the step be stupidly small
		if( step < MIN_STEP ) step = MIN_STEP;

		// calculate new position
		double dt = pt*step;
		double dr = pr*step;
		double dtheta = ptheta*step;
		double dphi = pphi*step;
		
		t += dt;
		r += dr;
		theta += dtheta;
		phi += dphi;

		if(r <= Raytracer<double>::horizon) break;

		ir = (logbin_r) ? static_cast<int>( log(r / r0) / log(bin_dr)) : static_cast<int>((r - r0) / bin_dr);
		itheta = static_cast<int>(theta / bin_dtheta);
		iphi = static_cast<int>((phi + M_PI) / bin_dphi);

		if(ir != last_ir && itheta != last_itheta && iphi != last_iphi
			&& ir > 0 && ir < Nr && itheta > 0 && itheta < Ntheta && iphi > 0 && iphi < Nphi)
		{
			if(motion == 1)
			{
				if (vel_mode == 0) V = vel;
				else if (vel_mode == 1) V = vel * (r / rmax);
				else if (vel_mode == 2) V = vel * sqrt(r / rmax);
			}
			else
				V = 1 / (a + r * sin(theta) * sqrt(r * sin(theta)));    // project the radius parallel to the equatorial plane

			redshift = Raytracer<double>::ray_redshift(V, false, false, r, theta, phi, k, h, Q, rdot_sign, thetadot_sign, Raytracer<double>::m_emit[ray], motion);
			(*map_time)[ir][itheta][iphi] += t;
			(*map_redshift)[ir][itheta][iphi] += redshift;
			++(*map_Nrays)[ir][itheta][iphi];
		}

	}

	Raytracer<double>::m_t[ray] = t;
	Raytracer<double>::m_r[ray] = r;
	Raytracer<double>::m_theta[ray] = theta;
	Raytracer<double>::m_phi[ray] = phi;
	Raytracer<double>::m_pt[ray] = pt;
	Raytracer<double>::m_pr[ray] = pr;
	Raytracer<double>::m_ptheta[ray] = ptheta;
	Raytracer<double>::m_pphi[ray] = pphi;
	Raytracer<double>::m_rdot_sign[ray] = rdot_sign;
	Raytracer<double>::m_thetadot_sign[ray] = thetadot_sign;

	return steps;
}


void Mapper::save(char* filename)
{
	ofstream outfile(filename, ios::binary);
	outfile.write(reinterpret_cast<char *> (&r0), sizeof(double));
	outfile.write(reinterpret_cast<char *> (&rmax), sizeof(double));
	outfile.write(reinterpret_cast<char *> (&Nr), sizeof(int));
	outfile.write(reinterpret_cast<char *> (&logbin_r), sizeof(bool));
	outfile.write(reinterpret_cast<char *> (&theta_max), sizeof(double));
	outfile.write(reinterpret_cast<char *> (&Ntheta), sizeof(int));
	outfile.write(reinterpret_cast<char *> (&Nphi), sizeof(int));

	map_time->write(&outfile);
	map_redshift->write(&outfile);
	map_Nrays->write(&outfile);
	bin_volume->write(&outfile);
	outfile.close();
}


void Mapper::average_rays()
{
	(*map_time) /= (*map_Nrays);
	(*map_redshift) /= (*map_Nrays);
}


void Mapper::calculate_volume()
{
	const double a = Raytracer<double>::spin;

	for(int ir=0; ir<Nr; ir++)
	{
		const double r = (logbin_r) ? r0 * pow(bin_dr, ir) : r0 + bin_dr*ir;
		const double this_bin_dr = (logbin_r) ? r*(bin_dr - 1) : bin_dr;

		for(int itheta=0; itheta<Ntheta; itheta++)
		{
			const double theta = itheta * bin_dtheta;

			const double rhosq = r * r + (a * cos(theta)) * (a * cos(theta));
			const double delta = r * r - 2 * r + a * a;
			const double sigmasq = (r * r + a * a) * (r * r + a * a) - a * a * delta * sin(theta) * sin(theta);
			const double e2nu = rhosq * delta / sigmasq;
			const double e2psi = sigmasq * sin(theta) * sin(theta) / rhosq;
			const double omega = 2 * a * r / sigmasq;
			const double grr = -rhosq / delta;
			const double gthth = -rhosq;
			const double gphph = -e2psi;

			for (int iphi = 0; iphi < Nphi; iphi++)
				(*bin_volume)[ir][itheta][iphi] = sqrt(-1 * grr * gthth * gphph) * this_bin_dr * bin_dtheta * bin_dphi;
		}
	}
}

