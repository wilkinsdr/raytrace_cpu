/*
 * raytracer.cu
 *
 *  Created on: 4 Sep 2013
 *      Author: drw
 */

#include "raytracer.h"
#include "raytrace_destination.h"
#include "vector"


template <typename T>
Raytracer<T>::Raytracer( int num_rays, T spin_par, T init_precision, T init_max_phistep, T init_max_tstep )
	: nRays( num_rays )
	,  spin(spin_par)
	, precision(init_precision)
    , theta_precision(THETA_PRECISION)
	, max_phistep(init_max_phistep)
	, max_tstep(init_max_tstep)
    , maxtstep_rlim(MAXDT_RLIM)
{
	//
	// Constructor function - allocates host and device memory for each ray to store ray position, momentum,
	// integration steps taken, redshift and constants of motion. Allocates only revice memory for the sign
	// of rdot and thetadot as well as the emitted energy used for redshift calculations.
	//
	// Rays are set up on a 2D grid of GPU threads (for easy variation of 2 parameters between rays)
	//
	// Arguments:
	//	num_rays	int		Number of rays required
	//	spin_par	T		Dimensionless spin parameter of black hole
	//	toler		T		Tolerance used to set step size in numerical integration of geodesics
	//

	// calculate horizon and store in GPU shared memory
	horizon = kerr_horizon<T>(spin);
	cout << "Event horizon at " << horizon << endl;

    int mb = 1<<20;

	cout << "Allocating memory (" << nRays*sizeof(Ray<T>) / mb << "MB)" << endl;
    rays = new Ray<T>[nRays];

    for(int ray=0; ray<nRays; ray++)
    {
	    rays[ray].steps = -1;
	    //rays[ray].status = 0;
    }
}

template <typename T>
Raytracer<T>::~Raytracer( )
{
	//
	// Destructor - frees host and device memory used for ray tracing variables
	//
	cout << "Cleaning up raytracer" << endl;

	delete[] rays;
}

// T t = rays[ray].t;
//        T r = rays[ray].r;
//        T theta = rays[ray].theta;
//        T phi = rays[ray].phi;
//        T pt = rays[ray].pt;
//        T pr = rays[ray].pr;
//        T ptheta = rays[ray].ptheta;
//        T pphi = rays[ray].pphi;
//        int rdot_sign = rays[ray].rdot_sign;
//        int thetadot_sign = rays[ray].thetadot_sign;
//
//        const T k = rays[ray].k;
//        const T h = rays[ray].h;
//        const T Q = rays[ray].Q;

//template<typename T>
//void Raytracer<T>::ray_sort(const int threads)
//{
//    rays_1 = new Ray<T>[nRays];
//
//    for (int i = 0; i < nRays; ++i) {
//        int index = i % threads;
//        rays_new_order[index].t = rays[i].t;
//
//    }
//
//    int ind = 0;
//    for (int ray_index = 0; ray_index < nRays / threads; ++ray_index) {
//        for (int index = 0; index < threads; ++index) {
//            rays[ind] = ray_groups[index][ray_index];
//            ++ind;
//        }
//    }
//}

//
//    for (int i = 0; i < threads; ++i) {
//        const int startIndex = i * raysPerThread;
//        for (int j = 0; j < raysPerThread; ++j) {
//            const int sourceIndex = j + startIndex;
//            const int targetIndex = sourceIndex + raysPerThread;
//            rays[targetIndex] = raysPerThreadArray[i][j];
//        }
//    }
//}

template <typename T>
void Raytracer<T>::run_raytrace(RayDestination<T>* destination, T r_max, T rad_disc, int show_progress, int status, T rad_in, T z_min, T bound, TextOutput* outfile, int write_step, T write_rmax, T write_rmin, bool write_cartesian)
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
	//	r_max		T		Limiting outer radius for propagation in Rg (default value = 1000)
	//
	cout << "Running raytracer..." << endl;

	//ProgressBar prog(nRays, "Ray", 0, (show_progress > 0));
    //show_progress = abs(show_progress);
//    #pragma omp parallel for
	for(int ray=0; ray<nRays; ray++)
	{
       //if(show_progress != 0 && (ray % show_progress) == 0) prog.show(ray+1);

//        int thread_num = omp_get_thread_num();
//        cout << "ray number " << ray << " running on thread " << thread_num << endl;

        if(rays[ray].steps < 0) continue;
		else if(rays[ray].steps >= STEPLIM) continue;

		//int n;
		//n = propagate(destination, ray, r_max, status, rad_disc, rad_in, z_min, bound,STEPLIM, outfile, write_step, write_rmax, write_rmin, write_cartesian);
		T rlim = r_max;
        T ray_status = status;
        T r_disc = rad_disc;
        T r_in = rad_in;
        T boundary = bound;
        T steplim = STEPLIM;



        //
        // propagate the photon along its geodesic until limiting r or theta reached
        //
        int steps = 0;

        int rsign_count = COUNT_MIN;
        int thetasign_count = COUNT_MIN;

        T rhosq, delta;
        T rdotsq, thetadotsq;

        T step;

        T x, y, z;

        // copy variables locally to simplify the code
        T a = spin;
        T t = rays[ray].t;
        T r = rays[ray].r;
        T theta = rays[ray].theta;
        T phi = rays[ray].phi;
        T pt = rays[ray].pt;
        T pr = rays[ray].pr;
        T ptheta = rays[ray].ptheta;
        T pphi = rays[ray].pphi;
        int rdot_sign = rays[ray].rdot_sign;
        int thetadot_sign = rays[ray].thetadot_sign;

        const T k = rays[ray].k;
        const T h = rays[ray].h;
        const T Q = rays[ray].Q;

        bool write_started = false;

        bool rdotsign_unlocked = false;
        bool thetadotsign_unlocked = false;

        bool under;

        // integrate geodesic equations until limit reached
        // if thetalim is positive, we go until theta exceeds it, if it is negative, we go until it is less than the abs value to allow tracing back to theta=0
        //while( r < rlim  && ( (thetalim > 0 && theta < thetalim) || (thetalim < 0 && theta > abs(thetalim)) || thetalim == 0 )  &&  steps < steplim )
        while(r < rlim  &&  steps < steplim)
        {
            ++steps;

            //const T thetalim = destination->get_thetalim(r, theta, phi, a);
            //T under = destination->under_function(r, theta);

            if(destination->stopping_fn(r, theta, phi, a)) {
                rays[ray].status = RAY_STOP_DEST;
                break;
            }

            rhosq = r*r + (a*cos(theta))*(a*cos(theta));
            delta = r*r - 2*r + a*a;

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

            if(rdotsign_unlocked && rdotsq <= 0 && rsign_count >= COUNT_MIN)
            {
                rdot_sign *= -1;
                rsign_count = 0;
                //continue;
            }
            else
            {
                rsign_count++;
                rdotsign_unlocked = true;
            }

            pr = sqrt(abs(rdotsq)) * rdot_sign;

            
	    step = abs((r - (T) horizon) / pr) / precision;

            if(step > abs(theta / ptheta) / precision)
            {
                step = abs(theta / ptheta) / theta_precision;
            }
            if(max_tstep > 0 && r < maxtstep_rlim && step > abs(max_tstep / pt))
            {
                step = abs(max_tstep / pt);
            }
            if(max_phistep > 0 && step > abs(max_phistep / pphi))
            {
                step = abs(max_phistep / pphi);
            }


            // don't let the step be stupidly small
            if(step < MIN_STEP) step = MIN_STEP;

            // make sure we don't go past rlim
            if(rlim > 0 && r + pr * step > rlim) step = abs((rlim - r) / pr);

            if(theta < 0)
            {
                theta *= -1;
                phi += M_PI;
                thetadot_sign *= -1;
            }

            if(step == 0)
            {
                rays[ray].status = 6;
                break;
            }
//        if(r <= horizon)
//        {
//            rays[ray].status = RAY_STOP_HORIZON;
//            break;
//        }
            if(rlim > 0 && r >= rlim)
            {
                rays[ray].status = RAY_STOP_RLIM;
                break;
            }
//        if(r_disc > 0 && thetalim > 0)
//        {
//            if(((!under) && theta >= thetalim) || (under && theta <= (M_PI - thetalim)))
//            {
//                if(r < r_disc && r >= r_in)
//                {
//                    rays[ray].status = RAY_STOP_DEST;
//                    break;
//                }
//            }
//        }
//        else if(thetalim > 0 && theta >= thetalim)
//        {
//            rays[ray].status = RAY_STOP_DEST;
//            break;
//        }
//        if(thetalim < 0 && thetalim > -4 && theta <= abs(thetalim))
//        {
//            rays[ray].status = RAY_STOP_DEST;
//            break;
//        }
//        if(zlim > 0 && r_disc > 0)
//        {
//            T z = r * cos(theta);
//            if(r < r_disc && r >= r_in)
//                if(abs(z) <= zlim)
//                {
//                    rays[ray].status = RAY_STOP_ZLIM;
//                    break;
//                }
//        }
//        else if(zlim > 0)
//        {
//            T z = r * cos(theta);
//            if(abs(z) <= zlim)
//            {
//                rays[ray].status = RAY_STOP_ZLIM;
//                break;
//            }
//        }
//        if(boundary > 0 && r <= boundary) {
//            rays[ray].status = RAY_STOP_BOUND;
//            break;
//        }


            // same for thetalim (but only if in range of r that would hit disc)

            //if(M_PI/4.0 > 0 && theta + ptheta * step > M_PI/4.0) step = abs((M_PI/4.0 - theta) / ptheta);

            // step condition
            //if(thetalim > 0 && theta + ptheta * step > thetalim) step = abs((thetalim - theta) / ptheta);


            if (theta + ptheta * step > M_PI_2) step = abs((M_PI_2 - theta) / ptheta);

            step = destination->step_function(r, theta, phi, step, ptheta);

            //step = step / 48;

            if(step < MIN_STEP) step = MIN_STEP;

            //cout << step << endl;
          // if (step > destination->step_function(r, theta, phi, step, ptheta, a)) {
	  // }

            // the code below is if we wanted to be able to trace rays below the disc (not implemented)
//        if((!under) && thetalim > 0 && theta + ptheta * step > thetalim)
//        {
//            // would the new step size put us in the range of r for the disc?
//            if(r_disc <= 0 || ((r + pr * abs((thetalim - theta) / ptheta)) < r_disc &&
//                               (r + pr * abs((thetalim - theta) / ptheta)) > r_in))
//                step = abs((thetalim - theta) / ptheta);
//        }
//        else if(r_disc > 0 && under && thetalim > 0 && theta + ptheta * step < (M_PI - thetalim))
//        {
//            if((r + pr * abs((thetalim - theta) / ptheta)) < r_disc &&
//               (r + pr * abs((thetalim - theta) / ptheta)) > r_in)
//                step = abs((thetalim - theta) / ptheta);
//        }

            // the code below is if we wanted to run rays backwards (not implemented)
            //if(reverse) step *= -1;

            // throw away the ray if tdot goes negative (inside the ergosphere - these are not physical)
            if(pt <= 0)
            {
                rays[ray].status = -2;
                rays[ray].steps = -1;
                break;
            }

            // calculate new position
//            if (flag){
//                t += pt*step;
//                r += pr*step;
//                phi += pphi*step;
//            } else {
                t += pt*step;
                r += pr*step;
                theta += ptheta*step;
                phi += pphi*step;
//            }
//            t += pt*step;
//            r += pr*step;
//            theta += ptheta*step;
//            phi += pphi*step;

            //aff += step;

            if(r <= horizon) break;

            if(outfile != 0 && (steps % write_step) == 0 )
            {
                if((write_rmax < 0 || r < write_rmax) && (write_rmin < 0 || r > write_rmin) )
                {
                    write_started = true;
                    if(write_cartesian)
                    {
                        cartesian<T>(x, y, z, r, theta, phi, a);
                        (*outfile) << t << x << y << z << endl;
                    }
                    else
                    {
                        (*outfile) << t << r << theta << phi << endl;
                    }
                }
                else if(write_started)
                {
                    break;
                }
            }
        }

        rays[ray].t = t;
        rays[ray].r = r;
        rays[ray].theta = theta;
        rays[ray].phi = phi;
        rays[ray].pt = pt;
        rays[ray].pr = pr;
        rays[ray].ptheta = ptheta;
        rays[ray].pphi = pphi;
        rays[ray].rdot_sign = rdot_sign;
        rays[ray].rdot_sign = thetadot_sign;

        if(steps > 0) rays[ray].steps += steps;





        //rays[ray].steps += n;

		if(outfile != 0)
			outfile->newline(2);
	}
//    prog.done();
}

template <typename T>
inline int Raytracer<T>::propagate(RayDestination<T>* destination, int ray, const T rlim, int ray_status, const T r_disc, const T r_in, const T zlim, const T boundary, const int steplim, TextOutput* outfile, int write_step, T write_rmax, T write_rmin, bool write_cartesian )
{
	//
	// propagate the photon along its geodesic until limiting r or theta reached
	//
	int steps = 0;

	int rsign_count = COUNT_MIN;
	int thetasign_count = COUNT_MIN;

	T rhosq, delta;
	T rdotsq, thetadotsq;

	T step;

	T x, y, z;

	// copy variables locally to simplify the code
	T a = spin;
	T t = rays[ray].t;
	T r = rays[ray].r;
	T theta = rays[ray].theta;
	T phi = rays[ray].phi;
	T pt = rays[ray].pt;
	T pr = rays[ray].pr;
	T ptheta = rays[ray].ptheta;
	T pphi = rays[ray].pphi;
	int rdot_sign = rays[ray].rdot_sign;
	int thetadot_sign = rays[ray].thetadot_sign;

	const T k = rays[ray].k;
	const T h = rays[ray].h;
	const T Q = rays[ray].Q;

	bool write_started = false;

	bool rdotsign_unlocked = false;
	bool thetadotsign_unlocked = false;

    bool under;

	// integrate geodesic equations until limit reached
	// if thetalim is positive, we go until theta exceeds it, if it is negative, we go until it is less than the abs value to allow tracing back to theta=0
    //while( r < rlim  && ( (thetalim > 0 && theta < thetalim) || (thetalim < 0 && theta > abs(thetalim)) || thetalim == 0 )  &&  steps < steplim )
    while(r < rlim  &&  steps < steplim)
	{
        ++steps;

        //const T thetalim = destination->get_thetalim(r, theta, phi, a);
        //T under = destination->under_function(r, theta);


        if(destination->stopping_fn(r, theta, phi, a)) {
            rays[ray].status = RAY_STOP_DEST;
            break;
        }

		rhosq = r*r + (a*cos(theta))*(a*cos(theta));
		delta = r*r - 2*r + a*a;

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

		if(rdotsign_unlocked && rdotsq <= 0 && rsign_count >= COUNT_MIN)
		{
			rdot_sign *= -1;
			rsign_count = 0;
			//continue;
		}
		else
		{
			rsign_count++;
			rdotsign_unlocked = true;
		}

		pr = sqrt(abs(rdotsq)) * rdot_sign;

        step = abs((r - (T) horizon) / pr) / precision;
        if(step > abs(theta / ptheta) / precision)
        {
            step = abs(theta / ptheta) / theta_precision;
        }
        if(max_tstep > 0 && r < maxtstep_rlim && step > abs(max_tstep / pt))
        {
            step = abs(max_tstep / pt);
        }
        if(max_phistep > 0 && step > abs(max_phistep / pphi))
        {
            step = abs(max_phistep / pphi);
        }
        // don't let the step be stupidly small
        if(step < MIN_STEP) step = MIN_STEP;

        // make sure we don't go past rlim
         if(rlim > 0 && r + pr * step > rlim) step = abs((rlim - r) / pr);

        if(theta < 0)
        {
            theta *= -1;
            phi += M_PI;
            thetadot_sign *= -1;
        }

        if(step == 0)
        {
            rays[ray].status = 6;
            break;
        }
//        if(r <= horizon)
//        {
//            rays[ray].status = RAY_STOP_HORIZON;
//            break;
//        }
        if(rlim > 0 && r >= rlim)
        {
            rays[ray].status = RAY_STOP_RLIM;
            break;
        }
//        if(r_disc > 0 && thetalim > 0)
//        {
//            if(((!under) && theta >= thetalim) || (under && theta <= (M_PI - thetalim)))
//            {
//                if(r < r_disc && r >= r_in)
//                {
//                    rays[ray].status = RAY_STOP_DEST;
//                    break;
//                }
//            }
//        }
//        else if(thetalim > 0 && theta >= thetalim)
//        {
//            rays[ray].status = RAY_STOP_DEST;
//            break;
//        }
//        if(thetalim < 0 && thetalim > -4 && theta <= abs(thetalim))
//        {
//            rays[ray].status = RAY_STOP_DEST;
//            break;
//        }
//        if(zlim > 0 && r_disc > 0)
//        {
//            T z = r * cos(theta);
//            if(r < r_disc && r >= r_in)
//                if(abs(z) <= zlim)
//                {
//                    rays[ray].status = RAY_STOP_ZLIM;
//                    break;
//                }
//        }
//        else if(zlim > 0)
//        {
//            T z = r * cos(theta);
//            if(abs(z) <= zlim)
//            {
//                rays[ray].status = RAY_STOP_ZLIM;
//                break;
//            }
//        }
//        if(boundary > 0 && r <= boundary) {
//            rays[ray].status = RAY_STOP_BOUND;
//            break;
//        }


        // same for thetalim (but only if in range of r that would hit disc)

        //if(M_PI/4.0 > 0 && theta + ptheta * step > M_PI/4.0) step = abs((M_PI/4.0 - theta) / ptheta);

        // step condition
        //if(thetalim > 0 && theta + ptheta * step > thetalim) step = abs((thetalim - theta) / ptheta);
//        if (step > destination->step_function(r, theta, phi, step, ptheta, a)) {
//	     step = destination->step_function(r, theta, phi, step, ptheta, a);
//	}
	// destination->step_function(r, theta, phi, step, ptheta, a);

        // the code below is if we wanted to be able to trace rays below the disc (not implemented)
//        if((!under) && thetalim > 0 && theta + ptheta * step > thetalim)
//        {
//            // would the new step size put us in the range of r for the disc?
//            if(r_disc <= 0 || ((r + pr * abs((thetalim - theta) / ptheta)) < r_disc &&
//                               (r + pr * abs((thetalim - theta) / ptheta)) > r_in))
//                step = abs((thetalim - theta) / ptheta);
//        }
//        else if(r_disc > 0 && under && thetalim > 0 && theta + ptheta * step < (M_PI - thetalim))
//        {
//            if((r + pr * abs((thetalim - theta) / ptheta)) < r_disc &&
//               (r + pr * abs((thetalim - theta) / ptheta)) > r_in)
//                step = abs((thetalim - theta) / ptheta);
//        }

        // the code below is if we wanted to run rays backwards (not implemented)
        //if(reverse) step *= -1;

		// throw away the ray if tdot goes negative (inside the ergosphere - these are not physical)
		if(pt <= 0)
		{
			//rays[ray].status =sort_rays -2;
            rays[ray].steps = -1;
			break;
		}

		// calculate new position
		t += pt*step;
		r += pr*step;
		theta += ptheta*step;
		phi += pphi*step;

		//aff += step;

		if(r <= horizon) break;

		if(outfile != 0 && (steps % write_step) == 0 )
		{
			if((write_rmax < 0 || r < write_rmax) && (write_rmin < 0 || r > write_rmin) )
			{
				write_started = true;
				if(write_cartesian)
				{
                    cartesian<T>(x, y, z, r, theta, phi, a);
					(*outfile) << t << x << y << z << endl;
				}
				else
				{
					(*outfile) << t << r << theta << phi << endl;
				}
			}
			else if(write_started)
			{
				break;
			}
		}
	}

	rays[ray].t = t;
    rays[ray].r = r;
    rays[ray].theta = theta;
    rays[ray].phi = phi;
    rays[ray].pt = pt;
    rays[ray].pr = pr;
    rays[ray].ptheta = ptheta;
    rays[ray].pphi = pphi;
    rays[ray].rdot_sign = rdot_sign;
    rays[ray].rdot_sign = thetadot_sign;

	if(steps > 0) rays[ray].steps += steps;

	return steps;
}

template <typename T>
void Raytracer<T>::redshift_start(T V, bool reverse, bool projradius )
{
	//
	// Calculates the initial energy of each ray at emission for use in redshift calculations.
	// Emitter is orbitting the black hole rotation axis at angular velocity V, or if V = -1, the
	// angular velocity is calculated for a circular orbit at the current r co-ordinate of the ray.
	//
	// If ray is being propagated backwards in time (for image planes), the spin parameter is revered
	// (back to its true value as it will already have been reversed for the propagation) and the
	// spatial components of the ray's 4-momentum will be reversed so the ray is travelling in the
	// correct direction wrt the emitting material.
	//
	// This function should be called prior to running the raytrace if redshifts are required.
	//
	// Arguments:
	//	V		T		Angular velocity of emitter (set to -1 for circular orbit at current r of ray)
	//	reverse	bool	Whether the ray is being propagated backwards in time (for image planes, default value = false)
	//
	cout << "Calculating initial energies" << endl;

	for(int ray=0; ray<nRays; ray++)
	{
		T p[4];

		const T a = (reverse) ? -1*spin : spin;

		// metric coefficients
		const T rhosq = rays[ray].r*rays[ray].r + (a*cos(rays[ray].theta))*(a*cos(rays[ray].theta));
		const T delta = rays[ray].r*rays[ray].r - 2*rays[ray].r + a*a;
		const T sigmasq = (rays[ray].r*rays[ray].r + a*a)*(rays[ray].r*rays[ray].r + a*a) - a*a*delta*sin(rays[ray].theta)*sin(rays[ray].theta);

		const T e2nu = rhosq * delta / sigmasq;
		const T e2psi = sigmasq * sin(rays[ray].theta)*sin(rays[ray].theta) / rhosq;
		const T omega = 2*a*rays[ray].r / sigmasq;

		T g[16];
		for(int i=0; i<16; i++)
			g[i] = 0;

		// g[i][j] -> g[i*4 + j]
		g[0*4 + 0] = e2nu - omega*omega*e2psi;
		g[0*4 + 3] = omega*e2psi;
		g[3*4 + 0] = g[0*4 + 3];
		g[1*4 + 1] = -rhosq/delta;
		g[2*4 + 2] = -rhosq;
		g[3*4 + 3] = -e2psi;

		// if V==-1, calculate orbital velocity for a geodesic circular orbit in equatorial plane
		if(V == -1 && projradius)
			V = 1 / (a + rays[ray].r*sin(rays[ray].theta)*sqrt(rays[ray].r*sin(rays[ray].theta)));	// project the radius parallel to the equatorial plane
		else if(V == -1)
			V = 1 / (a + rays[ray].r*sqrt(rays[ray].r));

		// if(reverse) V *= -1;


		// timelike basis vector
		const T et[] = { (1/sqrt(e2nu))/sqrt(1 - (V - omega)*(V - omega)*e2psi/e2nu)
                            , 0 , 0 ,
                            (1/sqrt(e2nu))*V / sqrt(1 - (V - omega)*(V - omega)*e2psi/e2nu) };



        // photon momentum
        momentum_from_consts<T>(p[0], p[1], p[2], p[3], rays[ray].k, rays[ray].h, rays[ray].Q, rays[ray].rdot_sign,
                                rays[ray].thetadot_sign, rays[ray].r, rays[ray].theta, rays[ray].phi, spin);

		// if we're propagating backwards, reverse the direction of the photon momentum
		if(reverse) p[1] *= -1; p[2] *= -1; p[3] *= -1;

		// evaluate dot product to get energy
		rays[ray].emit = 0;
		for(int i=0; i<4; i++)
			for(int j=0; j<4; j++)
				rays[ray].emit += g[i*4 + j] * et[i]* p[j];
	}
}


template <typename T>
void Raytracer<T>::redshift(RayDestination<T>* destination, T V, bool reverse, bool projradius, int motion )
{
	//
	// Calculates the redshift of the ray (emitted / received energy).
	// Receiver is orbitting the black hole rotation axis at angular velocity V, or if V = -1, the
	// angular velocity is calculated for a circular orbit at the current r co-ordinate of the ray.
	//
	// If ray is being propagated backwards in time (for image planes), the spin parameter is revered
	// (back to its true value as it will already have been reversed for the propagation) and the
	// spatial components of the ray's 4-momentum will be reversed so the ray is travelling in the
	// correct direction wrt the receiving material.
	//
	// redshift_start( ) needs to have been called before rays were traced. Use this function after
	// propagation.
	//
	// Arguments:
	//	V		T		Angular velocity of receiver (set to -1 for circular orbit at current r of ray)
	//	reverse	bool	Whether the ray is being propagated backwards in time (for image planes, default value = false)
	//	projradius bool	Whether to use r*sin(theta) as the radial co-ordinate in angular velocity calculation to use radius projected parallel to equatorial plane
	//
	cout << "Calculating ray redshifts..." << endl;

	for(int ray=0; ray<nRays; ray++)
	{
		rays[ray].redshift = ray_redshift(destination, V, reverse, projradius, rays[ray].r, rays[ray].theta, rays[ray].phi, rays[ray].k, rays[ray].h, rays[ray].Q, rays[ray].rdot_sign, rays[ray].thetadot_sign, rays[ray].emit, motion);
	}
}


template <typename T>
inline T Raytracer<T>::ray_redshift(RayDestination<T>* destination, T V, bool reverse, bool projradius, T r, T theta, T phi, T k, T h, T Q, int rdot_sign, int thetadot_sign, T emit, int motion )
{
	// calculate the redshift of a single ray

	T p[4];

    const T a = (reverse) ? -1*spin : spin;

	// metric coefficients
	const T rhosq = r*r + (a*cos(theta))*(a*cos(theta));
	const T delta = r*r - 2*r + a*a;
	const T sigmasq = (r*r + a*a)*(r*r + a*a) - a*a*delta*sin(theta)*sin(theta);

	const T e2nu = rhosq * delta / sigmasq;
	const T e2psi = sigmasq * sin(theta)*sin(theta) / rhosq;
	const T omega = 2*a*r / sigmasq;

	T g[16];
	for(int i=0; i<16; i++)
		g[i] = 0;

	// g[i][j] -> g[i*4 + j]
	g[0*4 + 0] = e2nu - omega*omega*e2psi;
	g[0*4 + 3] = omega*e2psi;
	g[3*4 + 0] = g[0*4 + 3];
	g[1*4 + 1] = -rhosq/delta;
	g[2*4 + 2] = -rhosq;
	g[3*4 + 3] = -e2psi;

	T et[] = {0, 0, 0, 0};

    destination->velocity_fn(et[0], et[1], et[2], et[3], r, theta, phi, a, h, k);

	// photon momentum
    momentum_from_consts<T>(p[0], p[1], p[2], p[3], k, h, Q, rdot_sign, thetadot_sign, r, theta, phi, spin);

	// if we're propagating backwards, reverse the direction of the photon momentum
	if(reverse)
	{
	    p[1] *= -1; p[2] *= -1; p[3] *= -1;
	}

	// evaluate dot product to get energy
	T recv = 0;
	for(int i=0; i<4; i++)
		for(int j=0; j<4; j++)
			recv += g[i*4 + j] * et[i] * p[j];

	return (reverse) ? recv / emit : emit / recv;
}


template <typename T>
void Raytracer<T>::range_phi(T min, T max )
{
	//
	// Puts the azimuthal angle co-ordinate (phi) in the required range
	//
	// Arguments:
	//	min		T		Lower bound of range (default value = -1 * M_PI)
	//	max		T		Upper bound of range (default value = M_PI)
	//
	cout << "Putting phi co-ordinate into range [" << min << "," << max << "]" << endl;
	for( int ray=0; ray<nRays; ray++ )
	{
		// check phi isn't something horrible so we don't enter an infinite loop
		if( abs(rays[ray].phi) > 1000 || rays[ray].phi != rays[ray].phi || !(rays[ray].steps>0) ) continue;

		while( rays[ray].phi >= max ) rays[ray].phi -= 2*M_PI;
		while( rays[ray].phi < min ) rays[ray].phi += 2*M_PI;
	}
}


template <typename T>
inline void Raytracer<T>::calculate_constants(int ray, T alpha, T beta, T V, T E)
{
	//
	// Compute the constants of motion for a ray emitted at polar angles alpha and beta in the frame
	// of a source at (rays[ray].t,rays[ray].r,rays[ray].theta,rays[ray].phi) orbiting azimuthally at angular velocity V
	//
	const T rhosq = rays[ray].r*rays[ray].r + (spin*cos(rays[ray].theta))*(spin*cos(rays[ray].theta));
	const T delta = rays[ray].r*rays[ray].r - 2*rays[ray].r + spin*spin;
	const T sigmasq = (rays[ray].r*rays[ray].r + spin*spin)*(rays[ray].r*rays[ray].r + spin*spin) - spin*spin*delta*sin(rays[ray].theta)*sin(rays[ray].theta);

	// metric coefficients
	const T e2nu = rhosq * delta / sigmasq;
	const T e2psi = sigmasq * sin(rays[ray].theta)*sin(rays[ray].theta) / rhosq;
	const T omega = 2*spin*rays[ray].r / sigmasq;

	// tetrad basis vector components
	const T et0 = (1/sqrt(e2nu))/sqrt(1 - (V - omega)*(V - omega)*e2psi/e2nu);
	const T et3 = (1/sqrt(e2nu))*V / sqrt(1 - (V - omega)*(V - omega)*e2psi/e2nu);
	//
	const T e10 = (V - omega)*sqrt(e2psi/e2nu) / sqrt(e2nu - (V-omega)*(V-omega)*e2psi);
	const T e13 = (1/sqrt(e2nu*e2psi))*(e2nu + V*omega*e2psi - omega*omega*e2psi) / sqrt(e2nu - (V-omega)*(V-omega)*e2psi);
	//
	const T e22 = -1/sqrt(rhosq);
	//
	const T e31 = sqrt(delta/rhosq);

	// photon 4-momentum in source frame
	const T rdotprime[] = { E, E*sin(alpha)*cos(beta), E*sin(alpha)*sin(beta), E*cos(alpha) };

	const T tdot = rdotprime[0]*et0 + rdotprime[1]*e10;
	const T phidot = rdotprime[0]*et3 + rdotprime[1]*e13;
	 const T rdot = rdotprime[3]*e31;
	const T thetadot = rdotprime[2]*e22;

	// find the corresponding values of k, h and Q using the geodesic equations
	rays[ray].k = (1 - 2*rays[ray].r/rhosq)*tdot + (2*spin*rays[ray].r*sin(rays[ray].theta)*sin(rays[ray].theta)/rhosq)*phidot;

	rays[ray].h = phidot * ( (rays[ray].r*rays[ray].r + spin*spin)*(rays[ray].r*rays[ray].r + spin*spin*cos(rays[ray].theta)*cos(rays[ray].theta) - 2*rays[ray].r)*sin(rays[ray].theta)*sin(rays[ray].theta) + 2*spin*spin*rays[ray].r*sin(rays[ray].theta)*sin(rays[ray].theta)*sin(rays[ray].theta)*sin(rays[ray].theta) );
	rays[ray].h = rays[ray].h - 2*spin*rays[ray].r*rays[ray].k*sin(rays[ray].theta)*sin(rays[ray].theta);
	rays[ray].h = rays[ray].h / ( rays[ray].r*rays[ray].r + spin*spin*cos(rays[ray].theta)*cos(rays[ray].theta) - 2*rays[ray].r );

	rays[ray].Q = rhosq*rhosq*thetadot*thetadot - (spin*rays[ray].k*cos(rays[ray].theta) + rays[ray].h/tan(rays[ray].theta))*(spin*rays[ray].k*cos(rays[ray].theta) - rays[ray].h/tan(rays[ray].theta));

	rays[ray].rdot_sign = (rdot >= 0) ? 1 : -1;
	rays[ray].thetadot_sign = (thetadot > 0) ? 1 : -1;

	//if(abs(rdot) < (1e-2 * e31)) rays[ray].steps = -1;
//	if(abs(rays[ray].r*phidot/rdot) > 1e3) rays[ray].steps = -1;
}

template <typename T>
inline void Raytracer<T>::calculate_constants_from_p(int ray, T pt, T pr, T ptheta, T pphi)
{
	//
	// calculate constants of motion from the 4-momentum and location of a photon
	//
	const T a = spin;
	const T r = rays[ray].r;
	const T theta = rays[ray].theta;

	const T rhosq = r*r + (a*cos(theta))*(a*cos(theta));

	T k = (1 - 2*r/rhosq)*pr + (2*a*r*sin(theta)*sin(theta)/rhosq)*pphi;

	T h = pphi * ( (r*r + a*a)*(r*r + a*a*cos(theta)*cos(theta) - 2*r)*sin(theta)*sin(theta) + 2*a*a*r*sin(theta)*sin(theta)*sin(theta)*sin(theta) );
	h = h - 2*a*r*k*sin(theta)*sin(theta);
	h = h / ( r*r + a*a*cos(theta)*cos(theta) - 2*r );

	T Q = rhosq*rhosq*ptheta*ptheta - (a*k*cos(theta) + h/tan(theta))*(a*k*cos(theta) - h/tan(theta));
	
	rays[ray].k = k;
	rays[ray].h = h;
	rays[ray].Q = Q;
}


template<typename T>
void Raytracer<T>::calculate_momentum( )
{
	//
	// Calculate photon momentum from constants of motion at a location
	//
	T a = spin;

	for(int ray = 0; ray < nRays; ray++)
	{
		const T t = rays[ray].t;
		const T r = rays[ray].r;
		const T theta = rays[ray].theta;
		const T phi = rays[ray].phi;
		const int rdot_sign = rays[ray].rdot_sign;
		const int thetadot_sign = rays[ray].thetadot_sign;

		const T k = rays[ray].k;
		const T h = rays[ray].h;
		const T Q = rays[ray].Q;

		const T rhosq = r*r + (a*cos(theta))*(a*cos(theta));
		const T delta = r*r - 2*r + a*a;

		// tdot
		rays[ray].pt = (rhosq*(r*r + a*a) + 2*a*a*r*sin(theta)*sin(theta))*k - 2*a*r*h;
		rays[ray].pt = rays[ray].pt / ( r*r * (1 + (a*cos(theta)/r)*(a*cos(theta)/r) - 2/r)*(r*r + a*a) + 2*a*a*r*sin(theta)*sin(theta) );

		// phidot
		rays[ray].pphi = 2*a*r*sin(theta)*sin(theta)*k + (r*r + (a*cos(theta))*(a*cos(theta)) - 2*r)*h;
		rays[ray].pphi = rays[ray].pphi / ( (r*r + a*a)*(r*r + (a*cos(theta))*(a*cos(theta)) - 2*r)*sin(theta)*sin(theta) + 2*a*a*r*sin(theta)*sin(theta)*sin(theta)*sin(theta) );

		// thetadot
		T thetadotsq = Q + (k*a*cos(theta) + h/tan(theta))*(k*a*cos(theta) - h/tan(theta));
		thetadotsq = thetadotsq / (rhosq*rhosq);

		// take the square roots and get the right signs
		rays[ray].ptheta = sqrt(abs(thetadotsq)) * thetadot_sign;

		// rdot
		T rdotsq = k*rays[ray].pt - h*rays[ray].pphi - rhosq*rays[ray].ptheta*rays[ray].ptheta;
		rdotsq = rdotsq * delta/rhosq;

		rays[ray].pr = sqrt(abs(rdotsq)) * rdot_sign;
	}
}

template class Raytracer<double>;
template class Raytracer<float>;
