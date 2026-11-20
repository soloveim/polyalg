#include <iostream>
#include <polyalg/polynomial.hpp>

int main() {
    using P = polyalg::Polynomial<double>;

    P p{-1, 0, 1};        // x^2 - 1
    P q{1, 1};            // x + 1
    P r{-2, 0, 1, 3};     // 3x^3 + x^2 - 2

    std::cout << "p        = " << p << "\n";
    std::cout << "q        = " << q << "\n";
    std::cout << "p + q    = " << p + q << "\n";
    std::cout << "p * q    = " << p * q << "\n";
    std::cout << "p / q    = " << p / q << "   (remainder " << p % q << ")\n";
    std::cout << "r / q    = " << r / q << "   (remainder " << r % q << ")\n";
    std::cout << "q^5      = " << q.pow(5) << "\n";
    std::cout << "p(q(x))  = " << p.compose(q) << "\n";
    std::cout << "r'       = " << r.derivative() << "\n";
    std::cout << "int r    = " << r.antiderivative() << "\n";
    std::cout << "int_0^1 r = " << r.integrate(0.0, 1.0) << "\n";
    std::cout << "r(2)     = " << r(2.0) << "\n";
    std::cout << "gcd(p*q, q^2*r) = " << P::gcd(p * q, q * q * r, 1e-9) << "\n";
    return 0;
}
