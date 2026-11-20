#include <cmath>
#include <iomanip>
#include <iostream>
#include <polyalg/quadrature.hpp>

int main() {
    using namespace polyalg;

    std::cout << "Legendre polynomials:\n";
    for (unsigned n = 0; n <= 5; ++n)
        std::cout << "  P_" << n << "(x) = " << legendre(n) << "\n";

    const unsigned n = 5;
    std::cout << std::setprecision(15);
    std::cout << "\nRoots of P_" << n << " (Gauss-Legendre nodes) and weights:\n";
    const auto rule = gauss_legendre(n);
    const auto P = legendre(n);
    for (unsigned i = 0; i < n; ++i)
        std::cout << "  x = " << std::setw(20) << rule.nodes[i]
                  << "   w = " << std::setw(18) << rule.weights[i]
                  << "   P_n(x) = " << P(rule.nodes[i]) << "\n";

    const double pi = std::acos(-1.0);
    std::cout << "\nIntegrals with a 8-point rule:\n";
    const auto g8 = gauss_legendre(8);
    std::cout << "  int_0^pi sin(x) dx = " << g8.integrate([](double x) { return std::sin(x); }, 0.0, pi)
              << "  (exact 2)\n";
    std::cout << "  int_0^1  exp(x) dx = " << g8.integrate([](double x) { return std::exp(x); }, 0.0, 1.0)
              << "  (exact " << std::exp(1.0) - 1.0 << ")\n";
    return 0;
}
