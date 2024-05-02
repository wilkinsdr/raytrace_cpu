/*
 * imagePlane.cpp
 *
 *  Created on: 25 Apr 2016
 *      Author: drw
 */

#include "imageplane.h"

template <typename T>
ImagePlane<T>::ImagePlane(T dist, T inc, T x0, T xmax, T dx, T y0, T ymax, T dy, T spin, T phi, T precision )
	: Raytracer<T>( (((xmax - x0) / dx) + 1) * (((ymax - y0) / dy) + 1) , spin , precision ),
	        Nx(((xmax - x0) / dx) + 1),
	        Ny(((ymax - y0) / dy) + 1),
	        D(dist),
	        incl(inc),
	        phi0(phi),
            m_x0(x0), m_xmax(xmax), m_dx(dx),
            m_y0(y0), m_ymax(ymax), m_dy(dy)
{
	cout << "Setting up image plane with (" << Nx << 'x' << Ny << ") rays" << endl;
	//init_image_plane(D, incl * M_PI / 180, phi0, x0, xmax, dx, y0, ymax, dy);
    init_image_plane(D, incl * M_PI / 180, phi0 * M_PI / 180, x0, xmax, dx, y0, ymax, dy);
}

//template <typename T>
//ImagePlane<T>::~ImagePlane()
//{
//
//}

//template<typename T>
//void ImagePlane<T>::ray_sort(const int threads, T rays[])
//{
//    const int num_rays = (Nx*Ny);
//    int ray_groups[threads][num_rays / threads];
//
//    for (int i = 0; i < num_rays; ++i) {
//        int index = i % threads;
//        int ray_index = i / threads;
//        ray_groups[index][ray_index] = rays[i];
//    }
//
//    int ind = 0;
//    for (int ray_index = 0; ray_index < num_rays / threads; ++ray_index) {
//        for (int index = 0; index < threads; ++index) {
//            rays[ind] = ray_groups[index][ray_index];
//            ++ind;
//        }
//    }
//}

template <typename T>
void ImagePlane<T>::init_image_plane(T D, T incl, T phi0,
									T x0, T xmax, T dx,
                                    T y0, T ymax, T dy)
{
	const int Nx = ((xmax - x0) / dx) + 1;
	const int Ny = ((ymax - y0) / dy) + 1;
	
	const double a = Raytracer<T>::spin;

//    int threads = 1;
//    #pragma omp parallel
//    threads = omp_get_num_threads();
    //cout << "Threads " << threads << endl;
    //const int threads = 8;
	for(int i=0; i<Nx; i++)
    {
        //const T x = x0 + ((i % threads)*threads + static_cast<int>(i / threads)) * dx;
        const T x = x0 + i * dx;
        //cout << "Ray " << i << " at x =" << x << endl;

		for (int j = 0; j < Ny; j++)
		{
			const int ix = i * Ny + j;

			const T y = y0 + j*dy;

			// initialise position of photon
			const T t = 0;
			const T r = sqrt(D*D + x*x + y*y);
			const T theta = acos( (D*cos(incl) + y*sin(incl)) / r );
			const T phi = phi0 + atan2( x, D*sin(incl) - y*cos(incl) );

            //cout << "phi " << phi << endl;
			// and the momentum
			const T pr = D/r;
			const T ptheta = sin(acos(D/r))/r;
			const T pphi = x*sin(incl)/(x*x+(D*sin(incl)-y*cos(incl))*(D*sin(incl)-y*cos(incl)));

			// metric coefficients
			const T rhosq = r*r + (a*cos(theta))*(a*cos(theta));
			const T delta = r*r - 2*r + a*a;
			const T sigmasq = (r*r + a*a)*(r*r + a*a) - a*a*delta*sin(theta)*sin(theta);

			const T e2nu = rhosq * delta / sigmasq;
			const T e2psi = sigmasq * sin(theta)*sin(theta) / rhosq;
			const T omega = 2*a*r / sigmasq;

			const T g00 = e2nu - omega*omega*e2psi;
			const T g03 = omega*e2psi;
			const T g11 = -rhosq/delta;
			const T g22 = -rhosq;
			const T g33 = -e2psi;

			// solve quadratic in pt to make this a null vector
			const T A = g00;
			const T B = 2*g03*pphi;
			const T C = g11*pr*pr + g22*ptheta*ptheta + g33*pphi*pphi;

			// take the appropriate (positive) root
			T pt = ( -B + sqrt(B*B - 4*A*C) ) / (2*A);
			if(pt < 0) pt = ( -B - sqrt(B*B - 4*A*C) ) / (2*A);

			Raytracer<T>::rays[ix].t = t;
			Raytracer<T>::rays[ix].r = r;
			Raytracer<T>::rays[ix].theta = theta;
			Raytracer<T>::rays[ix].phi = phi;
			Raytracer<T>::rays[ix].pt = pt;
			Raytracer<T>::rays[ix].pr = pr;
			Raytracer<T>::rays[ix].ptheta = ptheta;
			Raytracer<T>::rays[ix].pphi = pphi;

//            cout << "r " << Raytracer<T>::rays[ix].r << endl;
//            cout << "theta " << Raytracer<T>::rays[ix].theta << endl;
//            cout << "phi " << Raytracer<T>::rays[ix].phi << endl;
			// calculate constants of motion
			Raytracer<T>::calculate_constants_from_p(ix, pt, pr, ptheta, pphi);
			Raytracer<T>::rays[ix].rdot_sign = -1;
			Raytracer<T>::rays[ix].thetadot_sign = 1;

			Raytracer<T>::rays[ix].k = 1;

			const T b = sqrt(x*x + y*y);
			T beta = asin(y/b);
			if(x < 0) beta=M_PI-beta;

			T h = -1.*b*sin(incl)*cos(beta);
			T ltheta = b*sin(beta);
			T Q = (ltheta*ltheta) - (a*cos(theta))*(a*cos(theta))+((h/tan(theta)))*((h/tan(theta)));

			Raytracer<T>::rays[ix].h = h;
			Raytracer<T>::rays[ix].Q = Q;

			Raytracer<T>::rays[ix].thetadot_sign = (ltheta>=0) ? 1 : -1;
			
			Raytracer<T>::rays[ix].steps = 0;

            Raytracer<T>::rays[ix].alpha = x;
            Raytracer<T>::rays[ix].beta = y;
		}
	}
}

// template<typename T>
//    void ImagePlane<T>::init_log_imageplane(T *t, T *r, T *theta, T *phi,
//                                        T *pt, T *pr, T *ptheta, T *pphi,
//                                        T *k, T *h, T *Q, int *rdot_sign, int *thetadot_sign,
//                                        long *steps, int *status,
//                                        T *weight,
//                                        T D, T incl, T phi0,
//                                        T x0, T xmax, T dx,
//                                        T y0, T ymax, T dy,
//                                        T a,
//                                        int quad,
//                                        T excise)
//    {
// //        const int i = blockIdx.x * blockDim.x + threadIdx.x;
// //        const int j = blockIdx.y * blockDim.y + threadIdx.y;
// //
// //        // array index for this thread
// //        const int ix = j * blockDim.x * gridDim.x + i;
//
//        for(int i=0; i<Nx; i++)
//        {
//            //const T x = x0 + ((i % threads)*threads + static_cast<int>(i / threads)) * dx;
//            const T x = x0 + i * dx;
//            //cout << "Ray " << i << " at x =" << x << endl;
//
//            for (int j = 0; j < Ny; j++)
//            {
//                const int ix = i * Ny + j;
//
//                const T y = y0 + j*dy;
//
// //        T x = x0 * pow(dx, i);
// //        T y = y0 * pow(dy, j);
//
// //                if(x >= xmax || y >= ymax || x < x0 || y < y0)
// //                {
// //                    steps[ix] = -1;
// //                    weight[ix] = 0;
// //                    return;
// //                }
// //                if(x == 0 || y == 0)
// //                {
// //                    steps[ix] = -1;
// //                    status[ix] = -1;
// //                    weight[ix] = 0;
// //                    return;
// //                }
// //                if((x * x + y * y) < excise * excise)
// //                {
// //                    steps[ix] = -1;
// //                    status[ix] = -5;
// //                    weight[ix] = 0;
// //                    return;
// //                }
//
//                // flip signs of x and y according to the quadrant we're in
//                // numbered anti-clockwise from quadrant 0 (x+,y+)
//                if(quad == 1 || quad == 2) x *= -1;
//                if(quad == 2 || quad == 3) y *= -1;
//
// //                t[ix] = 0;
// //                r[ix] = sqrt(D * D + x * x + y * y);
// //                theta[ix] = acosf((D * cos(incl) + y * sin(incl)) / r[ix]);
// //                phi[ix] = phi0 + atan2f(x, D * sin(incl) - y * cos(incl));
//
//                // initialise position of photon
//                const T t = 0;
//                const T r = sqrt(D*D + x*x + y*y);
//                const T theta = acos( (D*cos(incl) + y*sin(incl)) / r );
//                const T phi = phi0 + atan2( x, D*sin(incl) - y*cos(incl) );
//
//                // and the momentum
//                const T pr = D/r;
//                const T ptheta = sin(acos(D/r))/r;
//                const T pphi = x*sin(incl)/(x*x+(D*sin(incl)-y*cos(incl))*(D*sin(incl)-y*cos(incl)));
//
//                // metric coefficients
//                const T rhosq = r*r + (a*cos(theta))*(a*cos(theta));
//                const T delta = r*r - 2*r + a*a;
//                const T sigmasq = (r*r + a*a)*(r*r + a*a) - a*a*delta*sin(theta)*sin(theta);
//
//                const T e2nu = rhosq * delta / sigmasq;
//                const T e2psi = sigmasq * sin(theta)*sin(theta) / rhosq;
//                const T omega = 2*a*r / sigmasq;
//
//                const T g00 = e2nu - omega*omega*e2psi;
//                const T g03 = omega*e2psi;
//                const T g11 = -rhosq/delta;
//                const T g22 = -rhosq;
//                const T g33 = -e2psi;
//
//                // solve quadratic in pt to make this a null vector
//                const T A = g00;
//                const T B = 2*g03*pphi;
//                const T C = g11*pr*pr + g22*ptheta*ptheta + g33*pphi*pphi;
//
//                // take the appropriate (positive) root
//                T pt = ( -B + sqrt(B*B - 4*A*C) ) / (2*A);
//                if(pt < 0) pt = ( -B - sqrt(B*B - 4*A*C) ) / (2*A);
//
//                Raytracer<T>::rays[ix].t = t;
//                Raytracer<T>::rays[ix].r = r;
//                Raytracer<T>::rays[ix].theta = theta;
//                Raytracer<T>::rays[ix].phi = phi;
//                Raytracer<T>::rays[ix].pt = pt;
//                Raytracer<T>::rays[ix].pr = pr;
//                Raytracer<T>::rays[ix].ptheta = ptheta;
//                Raytracer<T>::rays[ix].pphi = pphi;
//
//                // calculate constants of motion
//                Raytracer<T>::calculate_constants_from_p(ix, pt, pr, ptheta, pphi);
//                Raytracer<T>::rays[ix].rdot_sign = -1;
//                Raytracer<T>::rays[ix].thetadot_sign = 1;
//
//                Raytracer<T>::rays[ix].k = 1;
//
//                const T b = sqrt(x*x + y*y);
//                T beta = asin(y/b);
//                if(x < 0) beta=M_PI-beta;
//
//                T h = -1.*b*sin(incl)*cos(beta);
//                T ltheta = b*sin(beta);
//                T Q = (ltheta*ltheta) - (a*cos(theta))*(a*cos(theta))+((h/tan(theta)))*((h/tan(theta)));
//
//                Raytracer<T>::rays[ix].h = h;
//                Raytracer<T>::rays[ix].Q = Q;
//
//                Raytracer<T>::rays[ix].thetadot_sign = (ltheta>=0) ? 1 : -1;
//
// //                k[ix] = 1;
// //
// //                T b = sqrt(x * x + y * y);
// //                T beta = asinf(y / b);
// //                if(x < 0) beta = M_PI - beta;
// //
// //                h[ix] = -1. * b * sin(incl) * cos(beta);
// //
// //                T ltheta = b * sin(beta);
// //                Q[ix] = (ltheta * ltheta) - (a * cos(theta[ix])) * (a * cos(theta[ix])) +
// //                        ((h[ix] / tan(theta[ix]))) * ((h[ix] / tan(theta[ix])));
// //
// //                rdot_sign[ix] = -1;
// //                thetadot_sign[ix] = (ltheta >= 0) ? 1 : -1;
//
//                // the weighting of this ray is equal to the area of the grid element
//                //weight[ix] = x0 * (pow(dx,i+1) - pow(dx,i)) * y0 * (pow(dy,j+1) - pow(dy,j));
//                weight[ix] = abs(x) * (dx - 1) * abs(y) * (dy - 1);
//
// //                steps[ix] = 0;
// //                status[ix] = 0;
//                Raytracer<T>::rays[ix].steps = 0;
//                Raytracer<T>::rays[ix].status = 0;
//
//                Raytracer<T>::rays[ix].alpha = x;
//                Raytracer<T>::rays[ix].beta = y;
//    }
//}

template <typename T>
void ImagePlane<T>::redshift_start( )
{
	//
	// Call the redshift_start function of the base class using the source's angular velocity
	//
    Raytracer<T>::redshift_start(0, true);
}

template <typename T>
void ImagePlane<T>::redshift(RayDestination<T>* destination, bool projradius )
{
	//
	// Call the redshift_start function of the base class using the angular velocity for a circular orbit at the ray's end point
	// for rays incident on the accretion disc
	//
	Raytracer<T>::redshift(destination, -1, true, projradius);
}

template class ImagePlane<double>;
template class ImagePlane<float>;
