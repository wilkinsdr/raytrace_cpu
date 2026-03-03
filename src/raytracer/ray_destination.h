/*
 * ray_destination.h
 *
 *  Pluggable ray termination criteria and disc velocity fields for the RK4 ray tracer.
 *
 *  Subclass RayDestination<T> and implement:
 *    reached()  — return true to stop propagation at the current position.
 *    velocity() — return the angular velocity Ω = dφ/dt of disc material at the
 *                 termination point, for use by redshift(dest*). Return -1 to fall
 *                 back to the Keplerian circular-orbit velocity (default).
 */

#ifndef RAY_DESTINATION_H_
#define RAY_DESTINATION_H_

#include <cmath>

// Abstract base class — subclass this to define a custom ray termination surface
// and an associated velocity field.
//
// reached() is called after every RK4 position update; return true to stop the ray.
//
// velocity() returns the scalar angular velocity Ω = dφ/dt of material at (r,θ,φ).
// Return -1 for equatorial Keplerian (default). Used by four_velocity() default.
//
// four_velocity() populates et[4] = {ut, ur, uθ, uφ} (contravariant Boyer-Lindquist)
// for the material at (r,θ,φ). The default computes the circular-orbit 4-velocity from
// velocity(), using the Keplerian value if velocity() returns -1. Override to supply a
// fully general velocity field (e.g. warped disc, outflow) — the vector must satisfy
// g_μν et^μ et^ν = +1 (timelike, +--- signature convention used throughout the codebase).
// spin is the black-hole spin parameter a/M.
template <typename T>
class RayDestination {
public:
    virtual ~RayDestination() = default;
    virtual bool reached(T r, T theta, T phi) const = 0;
    // Direction-aware overload: thetadot_sign = +1 means theta is increasing (ray approaching
    // the equatorial plane from above); -1 means decreasing (approaching from below).
    // Default ignores direction and delegates to the positional overload.
    virtual bool reached(T r, T theta, T phi, int thetadot_sign) const {
        return reached(r, theta, phi);
    }
    virtual T velocity(T r, T theta, T phi) const { return -1; }
    virtual void four_velocity(T r, T theta, T phi, T spin, T et[4]) const
    {
        T V = velocity(r, theta, phi);

        const T rhosq   = r*r + (spin*cos(theta))*(spin*cos(theta));
        const T delta   = r*r - 2*r + spin*spin;
        const T sigmasq = (r*r + spin*spin)*(r*r + spin*spin)
                          - spin*spin*delta*sin(theta)*sin(theta);
        const T e2nu    = rhosq * delta / sigmasq;
        const T e2psi   = sigmasq * sin(theta)*sin(theta) / rhosq;
        const T omega   = 2*spin*r / sigmasq;

        if (V == -1) V = 1 / (spin + r * sqrt(r));   // equatorial Keplerian

        const T gamma_factor = 1 / sqrt(1 - (V - omega)*(V - omega)*e2psi/e2nu);
        et[0] = gamma_factor / sqrt(e2nu);
        et[1] = 0;
        et[2] = 0;
        et[3] = gamma_factor * V / sqrt(e2nu);
    }
};

// Concrete implementation equivalent to the existing theta_lim behaviour.
//   theta_lim > 0: stop when theta >= theta_lim  (e.g. M_PI/2 for equatorial disc)
//   theta_lim < 0: stop when theta <= |theta_lim| (reverse integration)
//   theta_lim = 0: never stop on theta (equivalent to passing thetalim=0 to propagate_rk4)
template <typename T>
class FlatDiscDestination : public RayDestination<T> {
    T theta_lim;
public:
    explicit FlatDiscDestination(T theta_lim = M_PI_2) : theta_lim(theta_lim) {}
    bool reached(T r, T theta, T phi) const override {
        if (theta_lim > 0) return theta >= theta_lim;
        if (theta_lim < 0) return theta <= -theta_lim;
        return false;
    }
};

// Like FlatDiscDestination but only triggers when r >= r_isco.
// Rays that cross theta_lim inside the ISCO pass through without stopping, allowing them
// to continue to the horizon.  Rays that have orbited the black hole and approach the
// disc from below (thetadot_sign < 0) ARE accepted — these produce the photon ring image.
// A ray inside the ISCO has no stable radial turning point and falls to the horizon, so
// the r >= r_isco guard alone is sufficient to exclude spurious inside-ISCO crossings.
template <typename T>
class DiscWithISCODestination : public RayDestination<T> {
    T theta_lim;
    T r_isco;
public:
    explicit DiscWithISCODestination(T r_isco, T theta_lim = M_PI_2)
        : theta_lim(theta_lim), r_isco(r_isco) {}
    bool reached(T r, T theta, T phi) const override {
        if (r < r_isco) return false;
        if (theta_lim > 0) return theta >= theta_lim;
        if (theta_lim < 0) return theta <= -theta_lim;
        return false;
    }
    bool reached(T r, T theta, T phi, int thetadot_sign) const override {
        if (r < r_isco) return false;
        // Only stop on the correct approach direction: thetadot_sign > 0 means theta is
        // increasing (approaching from above for theta_lim > 0); < 0 means from below.
        if (theta_lim > 0 && thetadot_sign <= 0) return false;
        if (theta_lim < 0 && thetadot_sign >= 0) return false;
        if (theta_lim > 0) return theta >= theta_lim;
        if (theta_lim < 0) return theta <= -theta_lim;
        return false;
    }

};

#endif /* RAY_DESTINATION_H_ */
