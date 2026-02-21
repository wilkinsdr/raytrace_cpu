/*
 * ray_destination.h
 *
 *  Pluggable ray termination criteria for the RK4 ray tracer.
 *
 *  Subclass RayDestination<T> and implement reached() to define a custom
 *  termination surface. The propagator calls reached(r, theta, phi) after
 *  every RK4 position update; return true to stop the ray at that point.
 */

#ifndef RAY_DESTINATION_H_
#define RAY_DESTINATION_H_

#include <cmath>

// Abstract base class — subclass this to define a custom ray termination surface.
// reached() receives the current Boyer-Lindquist coordinates of the ray.
// Return true to stop propagation at the current position.
template <typename T>
class RayDestination {
public:
    virtual ~RayDestination() = default;
    virtual bool reached(T r, T theta, T phi) const = 0;
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

#endif /* RAY_DESTINATION_H_ */
