//
// Created by Will Surgent on 8/8/23.
//

#include "log_imageplane.h"

template <typename T>
LogImagePlane<T>::LogImagePlane(T dist, T incl, T x0, T xmax, T dx, T y0, T ymax, T dy, T spin, int quad, T phi0, T excise, T precision, T max_tstep)
        : Raytracer<T>( (log(xmax/x0) / log(dx)) * (log(ymax/y0) / log(dy)), -1 * spin , precision ), //static_cast<int>( log(xmax / x0) / log(dx)) * static_cast<int>( log(ymax / y0) / log(dy)), spin , precision),
            Nx((log(xmax / x0) / log(dx))),
            Ny((log(ymax / y0) / log(dy))),
            D(dist),
            incl(incl),
            phi0(phi0),
            m_x0(x0), m_xmax(xmax), m_dx(dx),
            m_y0(y0), m_ymax(ymax), m_dy(dy),
            m_quad(quad)
{
    cout << endl << "Image plane at incl = " << incl << ", distance = " << dist << ", phi0 = " << phi0 << endl;
    cout << "Tracing rays with logarithmic spacing in range x = [" << x0 << ':' << m_dx << ':' << xmax << "], y = ["
         << y0 << ':' << m_dy << ':' << ymax << "], quadrant " << quad << endl << endl;
    //init_log_imageplane(D, incl * M_PI / 180, x0, xmax, dx, y0, ymax, dy, spin, quad, phi0); //0, PRECISION);
    init_log_imageplane(D, incl * M_PI / 180, x0, xmax, dx, y0, ymax, dy, spin, quad, phi0 * M_PI / 180); //0, PRECISION);
}

 template<typename T>
    void LogImagePlane<T>::init_log_imageplane(T dist, T incl, T x0, T xmax, T dx, T y0, T ymax, T dy, T a, int quad, T phi0,
                                               T excise, T precision)
    {
//        const int Nx = ((xmax - x0) / dx) + 1;
//        const int Ny = ((ymax - y0) / dy) + 1;

        const int Nx = (log(xmax/x0) / log(dx));
        const int Ny = (log(ymax/y0) / log(dy));

//        const int Nx = (log(xmax - x0) / log(dx));
//        const int Ny = (log(ymax - y0) / log(dy));

        for(int i=0; i<Nx; i++) {
            //double x = x0 + i * dx;
            double x = x0 * pow(dx, i);
            if (quad == 1 || quad == 2) x *= -1;

            for (int j = 0; j < Ny; j++) {
                const int ix = i * Ny + j;
                //const int ix = pow(Ny, i) + j;

                //double y = y0 + j * dy;
                double y = y0 * pow(dy, j);

                //        T x = x0 * pow(dx, i);
                //        T y = y0 * pow(dy, j);

                //                if(x >= xmax || y >= ymax || x < x0 || y < y0)
                //                {
                //                    steps[ix] = -1;
                //                    weight[ix] = 0;
                //                    return;
                //                }
                //                if(x == 0 || y == 0)
                //                {
                //                    steps[ix] = -1;
                //                    status[ix] = -1;
                //                    weight[ix] = 0;
                //                    return;
                //                }
                //                if((x * x + y * y) < excise * excise)
                //                {
                //                    steps[ix] = -1;
                //                    status[ix] = -5;
                //                    weight[ix] = 0;
                //                    return;
                //                }

                // flip signs of x and y according to the quadrant we're in
                // numbered anti-clockwise from quadrant 0 (x+,y+)
                if (quad == 2 || quad == 3) y *= -1;

                //                t[ix] = 0;
                //                r[ix] = sqrt(D * D + x * x + y * y);
                //                theta[ix] = acosf((D * cos(incl) + y * sin(incl)) / r[ix]);
                //                phi[ix] = phi0 + atan2f(x, D * sin(incl) - y * cos(incl));

                // initialise position of photon
                const T t = 0;
                const T r = sqrt(D * D + x * x + y * y);
                const T theta = acos((D * cos(incl) + y * sin(incl)) / r);
                const T phi = phi0 + atan2(x, D * sin(incl) - y * cos(incl));

                // and the momentum
                const T pr = D / r;
                const T ptheta = sin(acos(D / r)) / r;
                const T pphi =
                        x * sin(incl) / (x * x + (D * sin(incl) - y * cos(incl)) * (D * sin(incl) - y * cos(incl)));

                // metric coefficients
                const T rhosq = r * r + (a * cos(theta)) * (a * cos(theta));
                const T delta = r * r - 2 * r + a * a;
                const T sigmasq = (r * r + a * a) * (r * r + a * a) - a * a * delta * sin(theta) * sin(theta);

                const T e2nu = rhosq * delta / sigmasq;
                const T e2psi = sigmasq * sin(theta) * sin(theta) / rhosq;
                const T omega = 2 * a * r / sigmasq;

                const T g00 = e2nu - omega * omega * e2psi;
                const T g03 = omega * e2psi;
                const T g11 = -rhosq / delta;
                const T g22 = -rhosq;
                const T g33 = -e2psi;

                // solve quadratic in pt to make this a null vector
                const T A = g00;
                const T B = 2 * g03 * pphi;
                const T C = g11 * pr * pr + g22 * ptheta * ptheta + g33 * pphi * pphi;

                // take the appropriate (positive) root
                T pt = (-B + sqrt(B * B - 4 * A * C)) / (2 * A);
                if (pt < 0) pt = (-B - sqrt(B * B - 4 * A * C)) / (2 * A);

                Raytracer<T>::rays[ix].t = t;
                Raytracer<T>::rays[ix].r = r;
                Raytracer<T>::rays[ix].theta = theta;
                Raytracer<T>::rays[ix].phi = phi;
                Raytracer<T>::rays[ix].pt = pt;
                Raytracer<T>::rays[ix].pr = pr;
                Raytracer<T>::rays[ix].ptheta = ptheta;
                Raytracer<T>::rays[ix].pphi = pphi;

                // calculate constants of motion
                Raytracer<T>::calculate_constants_from_p(ix, pt, pr, ptheta, pphi);
                Raytracer<T>::rays[ix].rdot_sign = -1;
                Raytracer<T>::rays[ix].thetadot_sign = 1;

                Raytracer<T>::rays[ix].k = 1;

                const T b = sqrt(x * x + y * y);
                T beta = asin(y / b);
                if (x < 0) beta = M_PI - beta;

                T h = -1. * b * sin(incl) * cos(beta);
                T ltheta = b * sin(beta);
                T Q = (ltheta * ltheta) - (a * cos(theta)) * (a * cos(theta)) + ((h / tan(theta))) * ((h / tan(theta)));

                Raytracer<T>::rays[ix].h = h;
                Raytracer<T>::rays[ix].Q = Q;

                Raytracer<T>::rays[ix].thetadot_sign = (ltheta >= 0) ? 1 : -1;

                //                k[ix] = 1;
                //
                //                T b = sqrt(x * x + y * y);
                //                T beta = asinf(y / b);
                //                if(x < 0) beta = M_PI - beta;
                //
                //                h[ix] = -1. * b * sin(incl) * cos(beta);
                //
                //                T ltheta = b * sin(beta);
                //                Q[ix] = (ltheta * ltheta) - (a * cos(theta[ix])) * (a * cos(theta[ix])) +
                //                        ((h[ix] / tan(theta[ix]))) * ((h[ix] / tan(theta[ix])));
                //
                //                rdot_sign[ix] = -1;
                //                thetadot_sign[ix] = (ltheta >= 0) ? 1 : -1;

                // the weighting of this ray is equal to the area of the grid element
                //weight[ix] = x0 * (pow(dx,i+1) - pow(dx,i)) * y0 * (pow(dy,j+1) - pow(dy,j));

                //weight[ix] = abs(x) * (dx - 1) * abs(y) * (dy - 1);
                Raytracer<T>::rays[ix].weight = abs(x) * (dx - 1) * abs(y) * (dy - 1);

                //                steps[ix] = 0;
                //                status[ix] = 0;
                Raytracer<T>::rays[ix].steps = 0;
                Raytracer<T>::rays[ix].status = 0;

                // cout << ix << " " << i << " " << j << " " << x << " " << y << endl;
                Raytracer<T>::rays[ix].alpha = x;
                Raytracer<T>::rays[ix].beta = y;
            }
        }
}

template <typename T>
void LogImagePlane<T>::redshift_start( )
{
    //
    // Call the redshift_start function of the base class using the source's angular velocity
    //
    Raytracer<T>::redshift_start(0, true);
}

template <typename T>
void LogImagePlane<T>::redshift(RayDestination<T>* destination, bool projradius )
{
    //
    // Call the redshift_start function of the base class using the angular velocity for a circular orbit at the ray's end point
    // for rays incident on the accretion disc
    //
    Raytracer<T>::redshift(destination, -1, true, projradius);
}

template class LogImagePlane<double>;
template class LogImagePlane<float>;
