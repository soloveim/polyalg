#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <polyalg/quadrature.hpp>

static int failures = 0;
#define CHECK(cond) do { if (!(cond)) { ++failures; \
    std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK failed: " #cond "\n"; } } while (0)

using namespace polyalg;
static bool close(double a, double b, double eps = 1e-12) { return std::abs(a - b) < eps; }

int main() {
    // Legendre polynomials in closed form
    CHECK(legendre(0) == (Polynomial<double>{1.0}));
    CHECK(legendre(1) == (Polynomial<double>{0.0, 1.0}));
    {
        const auto p2 = legendre(2);              // (3x^2 - 1)/2
        CHECK(close(p2[0], -0.5) && close(p2[1], 0.0) && close(p2[2], 1.5));
        const auto p3 = legendre(3);              // (5x^3 - 3x)/2
        CHECK(close(p3[1], -1.5) && close(p3[3], 2.5) && p3.degree() == 3);
    }
    // P_n(1) = 1, recurrence evaluation agrees with the polynomial
    for (unsigned n = 0; n <= 12; ++n) {
        CHECK(close(legendre(n)(1.0), 1.0));
        for (double x : {-0.9, -0.3, 0.2, 0.77}) {
            CHECK(close(legendre_eval(n, x).first, legendre(n)(x)));
            CHECK(close(legendre_eval(n, x).second, legendre(n).derivative()(x), 1e-9));
        }
    }

    // known nodes / weights
    {
        auto g2 = gauss_legendre(2);
        CHECK(close(g2.nodes[0], -1.0 / std::sqrt(3.0)) && close(g2.nodes[1], 1.0 / std::sqrt(3.0)));
        CHECK(close(g2.weights[0], 1.0) && close(g2.weights[1], 1.0));
        auto g3 = gauss_legendre(3);
        CHECK(close(g3.nodes[1], 0.0) && close(g3.nodes[2], std::sqrt(0.6)));
        CHECK(close(g3.weights[0], 5.0 / 9.0) && close(g3.weights[1], 8.0 / 9.0));
    }

    for (unsigned n = 1; n <= 20; ++n) {
        const auto g = gauss_legendre(n);
        CHECK(g.nodes.size() == n && g.weights.size() == n);
        CHECK(std::is_sorted(g.nodes.begin(), g.nodes.end()));
        double wsum = 0;
        for (unsigned i = 0; i < n; ++i) {
            wsum += g.weights[i];
            CHECK(g.weights[i] > 0 && std::abs(g.nodes[i]) < 1.0);
            CHECK(close(g.nodes[i], -g.nodes[n - 1 - i]));       // symmetry
        }
        CHECK(close(wsum, 2.0));

        // exactness for monomials of degree <= 2n-1 on [-1, 1]
        for (unsigned k = 0; k <= 2 * n - 1; ++k) {
            const double exact = (k % 2 == 0) ? 2.0 / (k + 1) : 0.0;
            CHECK(close(g.integrate([k](double x) { return std::pow(x, k); }), exact, 1e-11));
        }
        // nodes are roots of P_n (for moderate n, where monomial form is fine)
        if (n <= 10) for (double x : g.nodes) CHECK(std::abs(legendre(n)(x)) < 1e-10);
    }

    // orthogonality: int P_m P_n = 2/(2n+1) delta_mn
    {
        const auto g = gauss_legendre(8);
        for (unsigned m = 0; m <= 5; ++m)
            for (unsigned n = 0; n <= 5; ++n) {
                const double v = g.integrate([&](double x) {
                    return legendre_eval(m, x).first * legendre_eval(n, x).first; });
                CHECK(close(v, m == n ? 2.0 / (2.0 * n + 1.0) : 0.0));
            }
    }

    // agreement with exact polynomial integration on [a, b]
    {
        Polynomial<double> p{1.0, -2.0, 0.5, 3.0, -1.0, 0.25};    // degree 5
        const auto g = gauss_legendre(3);                         // exact up to degree 5
        CHECK(close(g.integrate([&](double x) { return p(x); }, 0.3, 2.1), p.integrate(0.3, 2.1)));
    }

    // smooth integrands and error handling
    {
        const auto g = gauss_legendre(10);
        const double pi = std::acos(-1.0);
        CHECK(close(g.integrate([](double x) { return std::sin(x); }, 0.0, pi), 2.0, 1e-12));
        bool threw = false;
        try { (void)gauss_legendre(0); } catch (const std::invalid_argument&) { threw = true; }
        CHECK(threw);
    }

    // large n stays accurate through recurrence-based evaluation
    {
        const auto g = gauss_legendre(100);
        double wsum = 0; for (double w : g.weights) wsum += w;
        CHECK(close(wsum, 2.0, 1e-12));
        CHECK(close(g.integrate([](double x) { return x * x; }), 2.0 / 3.0, 1e-12));
    }

    if (failures == 0) std::cout << "All quadrature tests passed.\n";
    return failures == 0 ? 0 : 1;
}
