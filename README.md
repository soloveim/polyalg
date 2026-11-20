# polyalg

Header-only C++17 library for dense univariate polynomial algebra and
Gauss-Legendre quadrature.

## Polynomial algebra (`polyalg/polynomial.hpp`)
Construction, +, -, *, scalar multiply, division with remainder, `pow`,
`compose`, derivative, antiderivative, definite integral, Horner evaluation,
monic GCD, pretty printing.

    polyalg::Polynomial<double> p{1, 2, 3};   // 3x^2 + 2x + 1  (ascending coefficients)

## Legendre polynomials and Gaussian quadrature (`polyalg/quadrature.hpp`)

    #include <polyalg/quadrature.hpp>
    auto P5    = polyalg::legendre(5);            // Polynomial<double>
    auto roots = polyalg::legendre_roots(5);      // nodes, ascending
    auto rule  = polyalg::gauss_legendre(8);      // rule.nodes, rule.weights
    double I   = rule.integrate([](double x){ return std::sin(x); }, 0.0, 3.14159);

* `legendre(n)`: P_n via the three-term recurrence (monomial form; ill-conditioned for n > ~30).
* `legendre_eval(n, x)`: stable evaluation of P_n(x) and P_n'(x), used for all numerics.
* `legendre_roots(n)`: roots of P_n by Newton's method with Chebyshev-type starting values.
* `gauss_legendre(n)`: nodes = roots of P_n, weights w_i = 2 / ((1 - x_i^2) P_n'(x_i)^2);
  exact for polynomials up to degree 2n-1.

## Build

    cmake -S . -B build
    cmake --build build
    ctest --test-dir build --output-on-failure
    ./build/polyalg_demo
    ./build/quadrature_demo
