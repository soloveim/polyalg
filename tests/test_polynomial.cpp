#include <cmath>
#include <iostream>
#include <stdexcept>
#include <polyalg/polynomial.hpp>

static int failures = 0;
#define CHECK(cond) do { if (!(cond)) { ++failures; \
    std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK failed: " #cond "\n"; } } while (0)

using PI = polyalg::Polynomial<long long>;
using PD = polyalg::Polynomial<double>;

static bool close(double a, double b, double eps = 1e-9) { return std::abs(a - b) < eps; }

int main() {
    // normalization & degree
    CHECK(PI{}.degree() == -1);
    CHECK((PI{1, 2, 0, 0}).degree() == 1);
    CHECK(PI{0}.is_zero());

    // add / sub / mul
    PI a{1, 2, 3}, b{4, 5};
    CHECK(a + b == (PI{5, 7, 3}));
    CHECK(a - a == PI{});
    CHECK(a * b == (PI{4, 13, 22, 15}));
    CHECK(a * 0LL == PI{});
    CHECK(-a == (PI{-1, -2, -3}));

    // evaluation
    CHECK(a(2LL) == 1 + 4 + 12);
    CHECK(close(PD{1, 2, 3}(0.5), 1 + 1 + 0.75));

    // pow & compose
    CHECK((PI{1, 1}).pow(3) == (PI{1, 3, 3, 1}));
    CHECK(a.pow(0) == PI{1});
    CHECK((PI{0, 0, 1}).compose(PI{1, 1}) == (PI{1, 2, 1}));   // (x+1)^2

    // calculus
    CHECK(a.derivative() == (PI{2, 6}));
    CHECK(PI{5}.derivative() == PI{});
    PD c{0, 0, 3};  // 3x^2
    CHECK(close(c.integrate(0.0, 2.0), 8.0));

    // divmod: (x^2 - 1) / (x + 1)
    {
        PD p{-1, 0, 1}, d{1, 1};
        auto [q, r] = p.divmod(d);
        CHECK(q == (PD{-1, 1}));
        CHECK(r.is_zero());
    }
    // divmod identity with remainder
    {
        PD p{-2, 0, 1, 3}, d{1, 1};
        auto [q, r] = p.divmod(d);
        CHECK(r.degree() < d.degree());
        PD back = q * d + r;
        for (double x : {-2.0, 0.3, 1.7}) CHECK(close(back(x), p(x)));
    }
    // division by zero
    {
        bool threw = false;
        try { (void)PD{1, 2}.divmod(PD{}); } catch (const std::domain_error&) { threw = true; }
        CHECK(threw);
    }

    // gcd
    {
        PD common{-1, 1};                    // x - 1
        PD g = PD::gcd(common * PD{2, 1}, common * PD{3, 1}, 1e-9);
        CHECK(g.degree() == 1);
        CHECK(close(g[0], -1.0) && close(g[1], 1.0));
    }

    // printing
    CHECK((PI{1, -2, 3}).to_string() == "3x^2 - 2x + 1");
    CHECK((PI{0, -1}).to_string() == "-x");
    CHECK(PI{}.to_string() == "0");

    if (failures == 0) std::cout << "All tests passed.\n";
    return failures == 0 ? 0 : 1;
}
