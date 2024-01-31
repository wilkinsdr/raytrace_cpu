//
// Created by Will Surgent on 5/10/23.
//

#ifndef RAYTRACE_CPU_RAYTRACE_DESTINATION_H
#define RAYTRACE_CPU_RAYTRACE_DESTINATION_H

#include<cmath>
#include "../include/kerr.h"

template <typename T>
class RayDestination {
public:
    virtual bool stopping_fn(T r, T theta, T phi, T spin)
    { }
    virtual void velocity_fn(T& vt, T& vr, T& vtheta, T& vphi, T r, T theta, T phi, T spin, T h, T k)
    { }
    virtual double step_function(T r, T theta, T phi, T step, T ptheta, T pr, T pphi, T r_disc, T spin)
    { }
};

template <typename T>
class ZDestination:
        public RayDestination<T> {
private:
    T thetalim;
    T rlim;
public:
    ZDestination(T thetalim_val, T r_lim_val) {
        thetalim = thetalim_val;
        rlim = r_lim_val;
    }

    bool stopping_fn(T r, T theta, T phi, T spin) {
        return theta >= thetalim; //&& r >= rlim;
    }

    void velocity_fn(T& vt, T& vr, T& vtheta, T& vphi, T r, T theta, T phi, T spin, T h, T k) {
         const T rhosq = r*r + (spin*cos(theta))*(spin*cos(theta));
         const T delta = r*r - 2*r + spin*spin;
         const T sigmasq = (r*r + spin*spin)*(r*r + spin*spin) - spin*spin*delta*sin(theta)*sin(theta);

         const T e2nu = rhosq * delta / sigmasq;
         const T e2psi = sigmasq * sin(theta)*sin(theta) / rhosq;
         const T omega = 2*spin*r / sigmasq;

         const T V = 1 / (spin + r * sqrt(r));

         vt = (1 / sqrt(e2nu)) / sqrt(1 - (V - omega) * (V - omega) * e2psi / e2nu);
         vphi = (1 / sqrt(e2nu)) * V / sqrt(1 - (V - omega) * (V - omega) * e2psi / e2nu);
    }

    T step_function(T r, T theta, T phi, T step, T ptheta) {
        if (thetalim > 0 && theta + ptheta * step > thetalim) {
            return abs((thetalim - theta) / ptheta);
        }
    }
};

template <typename T>
class AngledDiscsDestination:
        public RayDestination<T> {
private:
    T thetalim_1;
    T thetalim_2;
    T r_inner;
public:
    AngledDiscsDestination(T thetalim_val_1, T thetalim_val_2, T r_inner_val) {
        thetalim_1 = thetalim_val_1;
        thetalim_2 = thetalim_val_2;
        r_inner = r_inner_val;
    }

    bool stopping_fn(T r, T theta, T phi, T spin) {
        double x, y, z;
        cartesian(x, y, z, r, theta, phi, spin);
        if (0 < phi && phi <= M_PI) {
            if (r < r_inner) {
                return theta >= M_PI_2;
            } else {
                return y*sin(thetalim_1) + z*cos(thetalim_1) <= 0;
            }
        } else {
            if (r < r_inner) {
                return theta >= M_PI_2;
            } else {
                return y*sin(thetalim_2) + z*cos(thetalim_2) <= 0;
            }
        }
    }

    void velocity_fn(T& vt, T& vr, T& vtheta, T& vphi, T r, T theta, T phi, T spin, T h, T k) {
        const T rhosq = r*r + (spin*cos(theta))*(spin*cos(theta));
        const T delta = r*r - 2*r + spin*spin;
        const T sigmasq = (r*r + spin*spin)*(r*r + spin*spin) - spin*spin*delta*sin(theta)*sin(theta);

        const T e2nu = rhosq * delta / sigmasq;
        const T e2psi = sigmasq * sin(theta)*sin(theta) / rhosq;
        const T omega = 2*spin*r / sigmasq;

        const T V = 1 / (spin + r * sqrt(r));

        const T a = spin;
        T pt = (rhosq*(r*r + a*a) + 2*a*a*r*sin(theta)*sin(theta))*k - 2*a*r*h;
		pt = pt / ( r*r * (1 + (a*cos(theta)/r)*(a*cos(theta)/r) - 2/r)*(r*r + a*a) + 2*a*a*r*sin(theta)*sin(theta) );

        if (0 < phi && phi <= M_PI) {
            if (r < r_inner) {
                vt = (1 / sqrt(e2nu)) / sqrt(1 - (V - omega) * (V - omega) * e2psi / e2nu);
                vr = 0;
                vtheta = 0;
                vphi = (1 / sqrt(e2nu)) * V / sqrt(1 - (V - omega) * (V - omega) * e2psi / e2nu);
            } else {
                T pphi = atan((sin(theta)*sin(phi)*cos(thetalim_1) - cos(theta)* sin(thetalim_1)) / (sin(theta)*cos(phi)));
                T vpphi = V*pt;

                vt = (1 / sqrt(e2nu)) / sqrt(1 - (V - omega) * (V - omega) * e2psi / e2nu);
                vr = 0;
                vtheta = (vpphi*cos(pphi)*sin(thetalim_1)) / sin(theta);
                vphi = (vpphi*cos(thetalim_1)*((1/cos(pphi)*1/cos(pphi)))) / ((1/cos(phi)*1/cos(phi)));
            }
        } else {
            if (r < r_inner) {
                vt = (1 / sqrt(e2nu)) / sqrt(1 - (V - omega) * (V - omega) * e2psi / e2nu);
                vr = 0;
                vtheta = 0;
                vphi = (1 / sqrt(e2nu)) * V / sqrt(1 - (V - omega) * (V - omega) * e2psi / e2nu);
            } else {
                T pphi = atan((sin(theta)*sin(phi)*cos(thetalim_2)- cos(theta)* sin(thetalim_2))/ (sin(theta)*cos(phi)));
                T vpphi = V*pt;

                vt = (1 / sqrt(e2nu)) / sqrt(1 - (V - omega) * (V - omega) * e2psi / e2nu);
                vr = 0;
                vtheta = (vpphi*cos(pphi)*sin(thetalim_2))/ sin(theta);
                vphi = (vpphi*cos(thetalim_2)*((1/cos(pphi)*1/cos(pphi)))) / ((1/cos(phi)*1/cos(phi)));
            }
        }
    }

    double step_function(T r, T theta, T phi, T step, T ptheta, T pr, T pphi, T r_disc, T spin) {
//        T thetalim = M_PI_2;
//        if (thetalim > 0 && theta + ptheta * step > thetalim) {
//            return abs((thetalim - theta) / ptheta);
//        }
//   }

        if (r < r_inner) {
            if (theta + ptheta * step > M_PI_2) {
                return abs((M_PI_2 - theta) / ptheta);
            } else {
                return step;
            }
        } else {
            double theta_new = acos(sin(theta)*sin(phi)*sin(thetalim_1) + cos(theta)*cos(thetalim_1));
            double dtheta = (((cos(thetalim_1)*sin(theta)) - (sin(thetalim_1)*cos(theta)*sin(phi))) /
                    sqrt(1 - ((sin(thetalim_1)*sin(theta)*sin(phi) + cos(thetalim_1)*cos(theta))*(sin(thetalim_1)*sin(theta)*sin(phi) + cos(thetalim_1)*cos(theta)))))*ptheta
                    - ((sin(thetalim_1)*sin(theta)*cos(phi))/
                    sqrt(1 - ((sin(thetalim_1)*sin(theta)*sin(phi) + cos(thetalim_1)*cos(theta))*(sin(thetalim_1)*sin(theta)*sin(phi) + cos(thetalim_1)*cos(theta)))))*pphi;
//            if (theta_new + dtheta * step > M_PI_2) {
//                return abs((theta_new - M_PI_2) / dtheta);
//            } else {
//                return step;
//            }
            return abs((M_PI_2 - theta_new) / dtheta);



//           if (theta + ptheta * step > thetalim_1) {    // change to use thetalim_1 and _2
//               return abs((thetalim_1 - theta) / ptheta);
//           } else {
//               return step;
//           }
//       }
//            double x, y, z;
//            cartesian(x, y, z, r, theta, phi, spin);
//            double dx = (r / sqrt(r*r + spin*spin))*pr*sin(theta)*cos(phi) + sqrt(r*r + spin*spin)*cos(theta)*cos(phi)*ptheta - sqrt(r*r + spin*spin)*sin(theta)*sin(phi)*pphi; //pr * sin(theta) * cos(phi) + r * ptheta * cos(theta) * cos(phi) - r * sin(theta) * pphi * sin(phi);
//            double dy = (r / sqrt(r*r + spin*spin))*pr*sin(theta)*sin(phi) + sqrt(r*r + spin*spin)*cos(theta)*sin(phi)*ptheta + sqrt(r*r + spin*spin)*sin(theta)*cos(phi)*pphi;//pr * sin(theta) * cos(phi) + r * ptheta * cos(theta) * cos(phi) + r * sin(theta) * pphi * cos(phi);
//            double dz = pr * cos(theta) - r * ptheta * sin(theta);

            //double xlim = sqrt(r_disc*r_disc - (y*y*cos(thetalim_1)*cos(thetalim_1)) + 2*y*z*cos(thetalim_1)*sin(thetalim_1) - (z*z*sin(thetalim_1)*sin(thetalim_1)));
            //double ylim = sqrt(abs((r_disc*r_disc*(1/cos(thetalim_1))*(1/cos(thetalim_1))) - (x*x*(1/cos(thetalim_1))*(1/cos(thetalim_1))))) + z*tan(thetalim_1);
            //double ylim = (2*z*cos(thetalim_1)*sin(thetalim_1) + sqrt(abs((2*z*cos(thetalim_1)*sin(thetalim_1))*(2*z*cos(thetalim_1)*sin(thetalim_1)) - 4*(cos(thetalim_1)*cos(thetalim_1)*(x*x + z*z*sin(thetalim_1)*sin(thetalim_1) - r_disc*r_disc))))) / (2*cos(thetalim_1)*cos(thetalim_1));
            //double zlim = (y*(1/tan(thetalim_1))) - (1/sin(thetalim_1))*(1/sin(thetalim_1))*sqrt((sin(thetalim_1)*sin(thetalim_1)*(r*r - x*x)));
            //double zlim = (2*y*cos(thetalim_1)*sin(thetalim_1) + sqrt(abs((2*y*cos(thetalim_1)*sin(thetalim_1))*(2*y*cos(thetalim_1)*sin(thetalim_1)) - 4*(sin(thetalim_1)*sin(thetalim_1)*(x*x + y*y*cos(thetalim_1)*cos(thetalim_1) - r_disc*r_disc))))) / (2*sin(thetalim_1)*sin(thetalim_1));

//            double xlim = r*cos(theta);
//            double ylim = r*cos(theta)*sin(thetalim_1);
//            double zlim = -r*sin(theta)*sin(thetalim_1);

//            double rlim = r_disc;
//            double thetalim = phi * sin(thetalim_1);

            //double step_r = abs((r - rlim) / pr);
            //double step_y = abs((theta - thetalim) / ptheta);

//            double step_r, step_theta;
//            if (r + pr * step > r_disc) {
//                step_r = abs((r - rlim) / pr);
//            } else {
//                return step;
//            }
//
//            if (theta + ptheta * step > thetalim) {
//                step_theta =  abs((theta - thetalim) / ptheta);
//            } else {
//                return step;
//            }
            //double step_z = abs((phi - philim) / pphi);

//            double step_x = abs((x - xlim) / dx);
//            double step_y = abs((y - ylim) / dy);
//            double step_z = abs((z - zlim) / dz);
//
//            if (step_x < step_y && step_x < step_z) {
//                return step_x;
//            } else if (step_y < step_x && step_y < step_z) {
//                return step_y;
//            } else {
//                return step_z;
//            }

//            if ((theta - thetalim) == 0) {
//                step_theta_zero_flag = true;
//                return step_r;
//            }
//
//            if (step_r < step_theta) {
//                return step_r;
//            } else {
//                return step_theta;
//            }
        }
    }
};

template <typename T>
class InclPortionDiscDestination:
        public RayDestination<T> {
private:
    T thetalim_1;
    T thetalim_2;
    T r_inner;
public:
    InclPortionDiscDestination(T thetalim_val_1, T thetalim_val_2, T r_inner_val) {
        thetalim_1 = thetalim_val_1;
        thetalim_2 = thetalim_val_2;
        r_inner = r_inner_val;
    }

    bool stopping_fn(T r, T theta, T phi, T spin) {
        double x, y, z;
        cartesian(x, y, z, r, theta, phi, spin);
        if (r > r_inner) {
            return y*sin(thetalim_1) + z*cos(thetalim_1) <= 0;
        }
//        double x, y, z;
//        cartesian(x, y, z, r, theta, phi, spin);
//        if (0 < phi && phi <= M_PI) {
//            if (r < r_inner) {
//                return false;
//            } else {
//                return y*sin(thetalim_1) + z*cos(thetalim_1) <= 0;
//            }
//        } else {
//            if (r < r_inner) {
//                return false;
//            } else {
//                return y*sin(thetalim_2) + z*cos(thetalim_2) <= 0;
//            }
//        }
    }

    void velocity_fn(T& vt, T& vr, T& vtheta, T& vphi, T r, T theta, T phi, T spin, T h, T k) {
        const T rhosq = r*r + (spin*cos(theta))*(spin*cos(theta));
        const T delta = r*r - 2*r + spin*spin;
        const T sigmasq = (r*r + spin*spin)*(r*r + spin*spin) - spin*spin*delta*sin(theta)*sin(theta);

        const T e2nu = rhosq * delta / sigmasq;
        const T e2psi = sigmasq * sin(theta)*sin(theta) / rhosq;
        const T omega = 2*spin*r / sigmasq;

        const T V = 1 / (spin + r * sqrt(r));

        const T a = spin;
        T pt = (rhosq*(r*r + a*a) + 2*a*a*r*sin(theta)*sin(theta))*k - 2*a*r*h;
        pt = pt / ( r*r * (1 + (a*cos(theta)/r)*(a*cos(theta)/r) - 2/r)*(r*r + a*a) + 2*a*a*r*sin(theta)*sin(theta) );

        if (0 < phi && phi <= M_PI) {
            if (r < r_inner) {
                vt = (1 / sqrt(e2nu)) / sqrt(1 - (V - omega) * (V - omega) * e2psi / e2nu);
                vr = 0;
                vtheta = 0;
                vphi = (1 / sqrt(e2nu)) * V / sqrt(1 - (V - omega) * (V - omega) * e2psi / e2nu);
            } else {
                T pphi = atan((sin(theta)*sin(phi)*cos(thetalim_1) - cos(theta)* sin(thetalim_1)) / (sin(theta)*cos(phi)));
                T vpphi = V*pt;

                vt = (1 / sqrt(e2nu)) / sqrt(1 - (V - omega) * (V - omega) * e2psi / e2nu);
                vr = 0;
                vtheta = (vpphi*cos(pphi)*sin(thetalim_1)) / sin(theta);
                vphi = (vpphi*cos(thetalim_1)*((1/cos(pphi)*1/cos(pphi)))) / ((1/cos(phi)*1/cos(phi)));
            }
        } else {
            if (r < r_inner) {
                vt = (1 / sqrt(e2nu)) / sqrt(1 - (V - omega) * (V - omega) * e2psi / e2nu);
                vr = 0;
                vtheta = 0;
                vphi = (1 / sqrt(e2nu)) * V / sqrt(1 - (V - omega) * (V - omega) * e2psi / e2nu);
            } else {
                T pphi = atan((sin(theta)*sin(phi)*cos(thetalim_2)- cos(theta)* sin(thetalim_2))/ (sin(theta)*cos(phi)));
                T vpphi = V*pt;

                vt = (1 / sqrt(e2nu)) / sqrt(1 - (V - omega) * (V - omega) * e2psi / e2nu);
                vr = 0;
                vtheta = (vpphi*cos(pphi)*sin(thetalim_2))/ sin(theta);
                vphi = (vpphi*cos(thetalim_2)*((1/cos(pphi)*1/cos(pphi)))) / ((1/cos(phi)*1/cos(phi)));
            }
        }
    }

    double step_function(T r, T theta, T phi, T step, T ptheta, T pr, T pphi, T r_disc, T spin) {
        return 0;
    }
};

template <typename T>
class TorusDiscDestination:
        public RayDestination<T> {
private:
    T r_torus;
    T r_disc;
    T r_inner;
public:
    TorusDiscDestination(T r_torus_val, T r_disc_val, T r_inner_val) {
        r_torus = r_torus_val;
        r_disc = r_disc_val;
        r_inner = r_inner_val;
    }

    bool stopping_fn(T r, T theta, T phi, T spin) {
       //T r_cen = (r_disc - r_torus);
       T r_cen = (r_inner + r_torus);
        double x, y, z;
        cartesian(x, y, z, r, theta, phi, spin);
        if (r <= r_inner) {
	     //return theta <= M_PI_2;
	     return z <= 0;
	    } else if (r > r_inner && r <= (r_inner + 2*r_torus)) {
             return z <= sqrt((r_torus * r_torus) - (r * sin(theta) - r_cen) * (r * sin(theta) - r_cen));
        } else {
             //return theta <= M_PI_2;
             return z <= 0;
        }
        //if (r <= (r_cen - r_torus)) {
        //    return theta >= M_PI_2;
        //} else {
        //    return z <= sqrt((r_torus * r_torus) - (r * sin(theta) - r_cen) * (r * sin(theta) - r_cen));
        //}
    }

    void velocity_fn(T& vt, T& vr, T& vtheta, T& vphi, T r, T theta, T phi, T spin, T h, T k) {
         const T rhosq = r*r + (spin*cos(theta))*(spin*cos(theta));
         const T delta = r*r - 2*r + spin*spin;
         const T sigmasq = (r*r + spin*spin)*(r*r + spin*spin) - spin*spin*delta*sin(theta)*sin(theta);

         const T e2nu = rhosq * delta / sigmasq;
         const T e2psi = sigmasq * sin(theta)*sin(theta) / rhosq;
         const T omega = 2*spin*r / sigmasq;

         const T V = 1 / (spin + (r*sin(theta)) * sqrt(r*sin(theta)));

         vt = (1 / sqrt(e2nu)) / sqrt(1 - (V - omega) * (V - omega) * e2psi / e2nu);
         vphi = (1 / sqrt(e2nu)) * V / sqrt(1 - (V - omega) * (V - omega) * e2psi / e2nu);
    }

    T step_function(T r, T theta, T phi, T step, T ptheta, T spin) {
        T r_cen = (r_disc - r_torus);
        double x, y, z;
        cartesian(x, y, z, r, theta, phi, spin);
        //T thetalim = asin((-(r_torus*r_torus)/(2*r_cen*r) + r/(2*r_cen) + r_cen/(2*r)));
        //T thetalim = asin(0.5 * sqrt(((-r_torus*r_torus + r*r + r_cen*r_cen) * (-r_torus*r_torus + r*r + r_cen*r_cen)) / (r*r*r_cen*r_cen)));
        //T thetalim = 0.5 * acos((r_torus*r_torus)/(r*r) + (r_torus*r_torus)/(r_cen*r_cen) - (r_torus*r_torus*r_torus*r_torus)/(2*r*r*r_cen*r_cen) - (r*r)/(2*r_cen*r_cen) - (r_cen*r_cen)/(2*r*r));
        T thetalim = M_PI_2;
        if (thetalim > 0 && theta + ptheta * step > thetalim) {
            return abs((thetalim - theta) / ptheta);
        }
    }
};

template <typename T>
class SinDiscDestination:
        public RayDestination<T> {
private:
    T r_disc;
public:
    SinDiscDestination(T r_disc_val) {
        r_disc = r_disc_val;
    }

    bool stopping_fn(T r, T theta, T phi, T spin) {
        double x, y, z;
        cartesian(x, y, z, r, theta, phi, spin);
        return z <= sin(sqrt(x*x + y*y));
    }

    void velocity_fn(T& vt, T& vr, T& vtheta, T& vphi, T r, T theta, T phi, T spin, T h, T k) {
         const T rhosq = r*r + (spin*cos(theta))*(spin*cos(theta));
         const T delta = r*r - 2*r + spin*spin;
         const T sigmasq = (r*r + spin*spin)*(r*r + spin*spin) - spin*spin*delta*sin(theta)*sin(theta);

         const T e2nu = rhosq * delta / sigmasq;
         const T e2psi = sigmasq * sin(theta)*sin(theta) / rhosq;
         const T omega = 2*spin*r / sigmasq;

         const T V = 1 / (spin + (r*sin(theta)) * sqrt(r*sin(theta)));

         vt = (1 / sqrt(e2nu)) / sqrt(1 - (V - omega) * (V - omega) * e2psi / e2nu);
         vphi = (1 / sqrt(e2nu)) * V / sqrt(1 - (V - omega) * (V - omega) * e2psi / e2nu);
    }

    T step_function(T r, T theta, T phi, T step, T ptheta, T spin) {
        T thetalim = M_PI_2;
        if (thetalim > 0 && theta + ptheta * step > thetalim) {
            return abs((thetalim - theta) / ptheta);
        }
    }
};

template <typename T>
class EllipseDiscDestination:
        public RayDestination<T> {
private:
    T r_disc;
    T r_inner;
    T major_axis;
    T minor_axis;
public:
    EllipseDiscDestination(T r_disc_val, T r_inner_val, T major_axis_val, T minor_axis_val) {
        r_disc = r_disc_val;
        r_inner = r_inner_val;
        major_axis = major_axis_val;
        minor_axis = minor_axis_val;
    }

    bool stopping_fn(T r, T theta, T phi, T spin) {
        double x, y, z;
        cartesian(x, y, z, r, theta, phi, spin);
        T h = (r_inner + major_axis);
        T val = sqrt(x*x + y*y);

	if (r >= r_inner && r <= (r_inner + 2*major_axis)) {
	     return z <= sqrt(minor_axis*minor_axis*(1 - ((val-h)*(val-h))/(major_axis*major_axis)));	
	} else if (r > r_inner && r > (r_inner + 2*major_axis)) {
	     return theta >= M_PI_2;
    }
    }

    void velocity_fn(T& vt, T& vr, T& vtheta, T& vphi, T r, T theta, T phi, T spin, T h, T k) {
         const T rhosq = r*r + (spin*cos(theta))*(spin*cos(theta));
         const T delta = r*r - 2*r + spin*spin;
         const T sigmasq = (r*r + spin*spin)*(r*r + spin*spin) - spin*spin*delta*sin(theta)*sin(theta);

         const T e2nu = rhosq * delta / sigmasq;
         const T e2psi = sigmasq * sin(theta)*sin(theta) / rhosq;
         const T omega = 2*spin*r / sigmasq;

         const T V = 1 / (spin + (r*sin(theta)) * sqrt(r*sin(theta)));

         vt = (1 / sqrt(e2nu)) / sqrt(1 - (V - omega) * (V - omega) * e2psi / e2nu);
         vphi = (1 / sqrt(e2nu)) * V / sqrt(1 - (V - omega) * (V - omega) * e2psi / e2nu);
    }

    T step_function(T r, T theta, T phi, T step, T ptheta, T pr, T pphi, T r_disc, T spin) {
        T thetalim = M_PI_2;
        if (theta + ptheta * step > thetalim) {
            return abs((thetalim - theta) / ptheta);
        } else {
            return step;
        }
    }
};


template <typename T>
class ShakuraDiscDestination:
        public RayDestination<T> {
private:
    T efficiency;
    T edd_frac;
    T r_isco;
public:
    ShakuraDiscDestination(T efficiency_val, T edd_frac_val, T r_isco_val) {
        efficiency = efficiency_val;
        edd_frac = edd_frac_val;
        r_isco = r_isco_val;
    }

    bool stopping_fn(T r, T theta, T phi, T spin) {
        double x, y, z;
        cartesian(x, y, z, r, theta, phi, spin);
        return z <= (3/2)*(1/efficiency)*(edd_frac)*(1 - sqrt(r_isco/(r*sin(theta))));
    }

    void velocity_fn(T& vt, T& vr, T& vtheta, T& vphi, T r, T theta, T phi, T spin, T h, T k) {
        const T rhosq = r*r + (spin*cos(theta))*(spin*cos(theta));
        const T delta = r*r - 2*r + spin*spin;
        const T sigmasq = (r*r + spin*spin)*(r*r + spin*spin) - spin*spin*delta*sin(theta)*sin(theta);

        const T e2nu = rhosq * delta / sigmasq;
        const T e2psi = sigmasq * sin(theta)*sin(theta) / rhosq;
        const T omega = 2*spin*r / sigmasq;

        const T V = 1 / (spin + (r*sin(theta)) * sqrt(r*sin(theta)));

        vt = (1 / sqrt(e2nu)) / sqrt(1 - (V - omega) * (V - omega) * e2psi / e2nu);
        vphi = (1 / sqrt(e2nu)) * V / sqrt(1 - (V - omega) * (V - omega) * e2psi / e2nu);
    }

    double step_function(T r, T theta, T phi, T step, T ptheta, T pr, T pphi, T r_disc, T spin) {
        //double x, y, z;
        //cartesian(x, y, z, r, theta, phi, spin);
        //double dz = pr* cos(theta) - r*sin(theta)*ptheta;
        //double zlim = 3.3;
        //if (z + dz * step > zlim) step = abs((zlim - z) / dz);
        double thetalim;
        thetalim = acos((-3*0.35)*(1-sqrt(r_isco/r)/(2*M_PI*r)));
        if (theta + ptheta * step > thetalim) step = abs((thetalim - theta) / ptheta);
        return step;
    }
};


#endif //RAYTRACE_CPU_RAYTRACE_DESTINATION_H
