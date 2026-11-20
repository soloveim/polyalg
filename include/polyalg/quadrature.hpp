#pragma once
// Legendre polynomials and Gauss-Legendre quadrature.
//
//  * legendre(n)          : P_n as a Polynomial<double> (three-term recurrence)
//  * legendre_eval(n, x)  : P_n(x) and P_n'(x) via the recurrence (stable, O(n))
//  * legendre_roots(n)    : the n roots of P_n in (-1, 1), ascending (Newton)
//  * gauss_legendre(n)    : nodes x_i (= roots of P_n) and weights
//                           w_i = 2 / ((1 - x_i^2) * P_n'(x_i)^2)
//
// An n-point rule integrates every polynomial of degree <= 2n-1 exactly.

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "polynomial.hpp"

namespace polyalg {

/// Legendre polynomial P_n(x) in monomial form:
///   (k+1) P_{k+1} = (2k+1) x P_k - k P_{k-1}.
/// Note: monomial coefficients become ill-conditioned for large n (>~ 30);
/// use legendre_eval() for numerics.
inline Polynomial<double> legendre(unsigned n) {
    using P = Polynomial<double>;
    P p0{1.0};
    if (n == 0) return p0;
    P p1{0.0, 1.0};
    const P x = P::x();
    for (unsigned k = 1; k < n; ++k) {
        const double kd = static_cast<double>(k);
        P p2 = (x * p1 * (2.0 * kd + 1.0) - p0 * kd) * (1.0 / (kd + 1.0));
        p0 = std::move(p1);
        p1 = std::move(p2);
    }
    return p1;
}

/// Returns {P_n(x), P_n'(x)} using the recurrence, valid for x in [-1, 1].
inline std::pair<double, double> legendre_eval(unsigned n, double x) {
    if (n == 0) return {1.0, 0.0};
    double p0 = 1.0, p1 = x;
    for (unsigned k = 1; k < n; ++k) {
        const double kd = static_cast<double>(k);
        const double p2 = ((2.0 * kd + 1.0) * x * p1 - kd * p0) / (kd + 1.0);
        p0 = p1;
        p1 = p2;
    }
    const double nd = static_cast<double>(n);
    double dp;
    if (std::abs(x) == 1.0) {                       // endpoint values
        dp = 0.5 * nd * (nd + 1.0);
        if (x < 0.0 && n % 2 == 0) dp = -dp;        // P_n'(-1) = (-1)^{n+1} n(n+1)/2
    } else {
        dp = nd * (x * p1 - p0) / (x * x - 1.0);
    }
    return {p1, dp};
}

/// Roots of P_n, ascending. Newton's method from the classical
/// cos(pi (i - 1/4) / (n + 1/2)) starting values; symmetry is enforced exactly.
inline std::vector<double> legendre_roots(unsigned n) {
    if (n == 0) throw std::invalid_argument("legendre_roots: n must be >= 1");
    const double pi = std::acos(-1.0);
    std::vector<double> r(n);
    const unsigned half = (n + 1) / 2;
    for (unsigned i = 0; i < half; ++i) {
        double x = std::cos(pi * (i + 0.75) / (n + 0.5));
        for (int it = 0; it < 100; ++it) {
            auto [p, dp] = legendre_eval(n, x);
            const double dx = p / dp;
            x -= dx;
            if (std::abs(dx) <= 4.0 * std::numeric_limits<double>::epsilon()) break;
        }
        if (n % 2 == 1 && i == half - 1) x = 0.0;   // exact middle root
        r[i] = -x;
        r[n - 1 - i] = x;
    }
    // for odd n the middle element was written twice with the same value 0
    return r;
}

struct GaussLegendre {
    std::vector<double> nodes;     // ascending, in (-1, 1)
    std::vector<double> weights;   // positive, sum to 2

    /// Approximate the integral of f over [a, b] (affine change of variable).
    template <typename F>
    double integrate(F&& f, double a = -1.0, double b = 1.0) const {
        const double mid = 0.5 * (a + b), half = 0.5 * (b - a);
        double s = 0.0;
        for (std::size_t i = 0; i < nodes.size(); ++i)
            s += weights[i] * f(mid + half * nodes[i]);
        return half * s;
    }
};

/// n-point Gauss-Legendre rule on [-1, 1].
inline GaussLegendre gauss_legendre(unsigned n) {
    GaussLegendre g;
    g.nodes = legendre_roots(n);
    g.weights.resize(n);
    for (unsigned i = 0; i < n; ++i) {
        const double x = g.nodes[i];
        const double dp = legendre_eval(n, x).second;
        g.weights[i] = 2.0 / ((1.0 - x * x) * dp * dp);
    }
    return g;
}

} // namespace polyalg
