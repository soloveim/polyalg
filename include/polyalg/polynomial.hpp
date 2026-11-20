#pragma once
// Dense univariate polynomial algebra, header-only.
//
// Coefficients are stored in ascending order: c[i] is the coefficient of x^i.
// The representation is always normalized: no trailing zero coefficients;
// the zero polynomial is the empty vector (degree() == -1).
//
// T may be any arithmetic-like type supporting + - * / (division only for
// divmod/gcd), comparison with T{}, and std::abs via ADL (for gcd/trimmed).

#include <algorithm>
#include <cstddef>
#include <initializer_list>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <cmath>
#include <cstdlib>

namespace polyalg {

template <typename T>
class Polynomial {
public:
    using value_type = T;

    // ---- construction ----------------------------------------------------
    Polynomial() = default;
    Polynomial(T constant) { if (constant != T{}) c_.push_back(constant); }
    Polynomial(std::initializer_list<T> coeffs) : c_(coeffs) { normalize(); }
    explicit Polynomial(std::vector<T> coeffs) : c_(std::move(coeffs)) { normalize(); }

    /// The monomial coeff * x^power.
    static Polynomial monomial(std::size_t power, T coeff = T{1}) {
        std::vector<T> v(power + 1, T{});
        v[power] = coeff;
        return Polynomial(std::move(v));
    }
    /// The polynomial x.
    static Polynomial x() { return monomial(1); }

    // ---- inspection ------------------------------------------------------
    std::ptrdiff_t degree() const { return static_cast<std::ptrdiff_t>(c_.size()) - 1; }
    bool is_zero() const { return c_.empty(); }
    std::size_t size() const { return c_.size(); }
    const std::vector<T>& coefficients() const { return c_; }
    T operator[](std::size_t i) const { return i < c_.size() ? c_[i] : T{}; }
    T leading() const { return c_.empty() ? T{} : c_.back(); }

    // ---- evaluation (Horner) --------------------------------------------
    template <typename U>
    auto operator()(const U& x) const -> decltype(T{} * x) {
        decltype(T{} * x) r{};
        for (auto it = c_.rbegin(); it != c_.rend(); ++it) r = r * x + *it;
        return r;
    }

    // ---- arithmetic ------------------------------------------------------
    Polynomial operator-() const {
        Polynomial r(*this);
        for (auto& v : r.c_) v = -v;
        return r;
    }

    Polynomial& operator+=(const Polynomial& o) {
        if (o.c_.size() > c_.size()) c_.resize(o.c_.size(), T{});
        for (std::size_t i = 0; i < o.c_.size(); ++i) c_[i] += o.c_[i];
        normalize();
        return *this;
    }
    Polynomial& operator-=(const Polynomial& o) {
        if (o.c_.size() > c_.size()) c_.resize(o.c_.size(), T{});
        for (std::size_t i = 0; i < o.c_.size(); ++i) c_[i] -= o.c_[i];
        normalize();
        return *this;
    }
    Polynomial& operator*=(const Polynomial& o) { return *this = *this * o; }
    Polynomial& operator*=(const T& s) {
        for (auto& v : c_) v *= s;
        normalize();
        return *this;
    }

    friend Polynomial operator+(Polynomial a, const Polynomial& b) { return a += b; }
    friend Polynomial operator-(Polynomial a, const Polynomial& b) { return a -= b; }
    friend Polynomial operator*(const Polynomial& a, const Polynomial& b) {
        if (a.is_zero() || b.is_zero()) return {};
        std::vector<T> r(a.c_.size() + b.c_.size() - 1, T{});
        for (std::size_t i = 0; i < a.c_.size(); ++i)
            for (std::size_t j = 0; j < b.c_.size(); ++j)
                r[i + j] += a.c_[i] * b.c_[j];
        return Polynomial(std::move(r));
    }
    friend Polynomial operator*(Polynomial a, const T& s) { return a *= s; }
    friend Polynomial operator*(const T& s, Polynomial a) { return a *= s; }

    /// Quotient and remainder: *this = q * d + r, deg r < deg d.
    /// Requires an exact division on T (a field, or floating point).
    std::pair<Polynomial, Polynomial> divmod(const Polynomial& d) const {
        if (d.is_zero()) throw std::domain_error("Polynomial division by zero");
        if (degree() < d.degree()) return {Polynomial{}, *this};
        Polynomial r(*this);
        std::vector<T> q(static_cast<std::size_t>(degree() - d.degree()) + 1, T{});
        const T lead = d.leading();
        for (std::ptrdiff_t k = degree() - d.degree(); k >= 0; --k) {
            const std::size_t top = static_cast<std::size_t>(k) + d.c_.size() - 1;
            if (top >= r.c_.size()) continue;
            T f = r.c_[top] / lead;
            q[static_cast<std::size_t>(k)] = f;
            for (std::size_t j = 0; j < d.c_.size(); ++j)
                r.c_[static_cast<std::size_t>(k) + j] -= f * d.c_[j];
            r.c_[top] = T{};          // kill rounding residue at the pivot
            r.normalize();
        }
        return {Polynomial(std::move(q)), std::move(r)};
    }
    friend Polynomial operator/(const Polynomial& a, const Polynomial& b) { return a.divmod(b).first; }
    friend Polynomial operator%(const Polynomial& a, const Polynomial& b) { return a.divmod(b).second; }

    /// Integer power by repeated squaring.
    Polynomial pow(unsigned n) const {
        Polynomial result{T{1}}, base(*this);
        while (n) {
            if (n & 1u) result *= base;
            n >>= 1u;
            if (n) base *= base;
        }
        return result;
    }

    /// Composition: (*this)(inner(x)), evaluated by Horner's scheme.
    Polynomial compose(const Polynomial& inner) const {
        Polynomial r;
        for (auto it = c_.rbegin(); it != c_.rend(); ++it) r = r * inner + Polynomial(*it);
        return r;
    }

    // ---- calculus --------------------------------------------------------
    Polynomial derivative() const {
        if (c_.size() <= 1) return {};
        std::vector<T> r(c_.size() - 1);
        for (std::size_t i = 1; i < c_.size(); ++i) r[i - 1] = c_[i] * static_cast<T>(i);
        return Polynomial(std::move(r));
    }
    /// Antiderivative with zero constant term (needs division on T).
    Polynomial antiderivative() const {
        if (is_zero()) return {};
        std::vector<T> r(c_.size() + 1, T{});
        for (std::size_t i = 0; i < c_.size(); ++i) r[i + 1] = c_[i] / static_cast<T>(i + 1);
        return Polynomial(std::move(r));
    }
    T integrate(const T& a, const T& b) const {
        const Polynomial F = antiderivative();
        return F(b) - F(a);
    }

    // ---- utilities -------------------------------------------------------
    /// Copy with leading coefficients of magnitude <= eps removed.
    Polynomial trimmed(T eps) const {
        using std::abs;
        Polynomial r(*this);
        while (!r.c_.empty() && !(abs(r.c_.back()) > eps)) r.c_.pop_back();
        return r;
    }
    Polynomial monic() const {
        if (is_zero()) return {};
        Polynomial r(*this);
        const T l = leading();
        for (auto& v : r.c_) v = v / l;
        return r;
    }

    /// Monic greatest common divisor via Euclid. For floating point,
    /// remainders are truncated at `eps` to absorb rounding noise.
    static Polynomial gcd(Polynomial a, Polynomial b, T eps = T{}) {
        a = a.trimmed(eps);
        b = b.trimmed(eps);
        while (!b.is_zero()) {
            Polynomial r = (a % b).trimmed(eps);
            a = std::move(b);
            b = std::move(r);
        }
        return a.monic();
    }

    // ---- comparison / output --------------------------------------------
    friend bool operator==(const Polynomial& a, const Polynomial& b) { return a.c_ == b.c_; }
    friend bool operator!=(const Polynomial& a, const Polynomial& b) { return !(a == b); }

    std::string to_string(const std::string& var = "x") const {
        if (c_.empty()) return "0";
        std::ostringstream os;
        bool first = true;
        for (std::size_t k = c_.size(); k-- > 0;) {
            const T& v = c_[k];
            if (v == T{}) continue;
            const bool neg = v < T{};
            T mag = neg ? -v : v;
            if (first) { if (neg) os << "-"; }
            else os << (neg ? " - " : " + ");
            first = false;
            if (k == 0 || mag != T{1}) os << mag;
            if (k >= 1) os << var;
            if (k >= 2) os << "^" << k;
        }
        return os.str();
    }
    friend std::ostream& operator<<(std::ostream& os, const Polynomial& p) { return os << p.to_string(); }

private:
    void normalize() { while (!c_.empty() && c_.back() == T{}) c_.pop_back(); }
    std::vector<T> c_;
};

using Polynomiald = Polynomial<double>;
using Polynomiali = Polynomial<long long>;

} // namespace polyalg
