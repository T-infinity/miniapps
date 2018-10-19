#pragma once
#include <vector>
#include <array>
#include <iostream>
#include <cmath>
#include <complex>

namespace Linearize {

    template<typename T, size_t n>
    class Ddata {
    public:
        T value;
        std::array<T, n> derivatives;

        explicit Ddata(){
            value = 0.0;
            derivatives.fill(0.0);
        }

        explicit Ddata(T f) : value(f) {
            derivatives.fill(0.0);
        }

        explicit Ddata(T f_in, const std::array<T, n>& d_in) : value(f_in), derivatives(d_in) {}

        Ddata<T, n>& operator=(T rhs) {
            value = rhs;
            derivatives.fill(0.0);
            return *this;
        }
        void operator+=(const Ddata<T, n> &rhs) {
            value += rhs.value;
            for (size_t i = 0; i < n; i++)
                derivatives[i] += rhs.derivatives[i];
        }
        void operator+=(T rhs) {
            value += rhs;
        }
        void operator-=(const Ddata<T, n> &rhs) {
            value -= rhs.value;
            for (size_t i = 0; i < n; i++)
                derivatives[i] -= rhs.derivatives[i];
        }
        void operator-=(T rhs) {
            value -= rhs;
        }
        void operator*=(const Ddata<T, n> &rhs) {
            value *= rhs.value;
            for (size_t i = 0; i < n; i++)
                derivatives[i] = value * rhs.derivatives[i] + derivatives[i] * rhs.value;
        }
        void operator*=(T rhs) {
            value *= rhs;
            for (size_t i = 0; i < n; i++)
                derivatives[i] *= rhs;
        }
        void operator/=(const Ddata<T, n> &rhs) {
            value /= rhs.value;
            T denom = 1.0 / rhs.value / rhs.value;
            for (size_t i = 0; i < n; i++)
                derivatives[i] = (derivatives[i] * rhs.value - rhs.derivatives[i] * value) * denom;
        }
        void operator/=(T rhs) {
            value /= rhs;
            for (size_t i = 0; i < n; i++)
                derivatives[i] /= rhs;
        }
        Ddata<T, n> operator-() const {
            auto ddata = *this;
            auto local_copy = ddata.derivatives;
            for (auto& entry : local_copy)
                entry = -entry;
            return Ddata<T, n>(-ddata.value, local_copy);
        }

        bool operator>(const Ddata& b) const {
            return (value > b.value);
        }
        bool operator>(T b) const {
            return (value > b);
        }

        bool operator<(const Ddata& a) const {
            return (value < a.value);
        }
        bool operator<(T b) const {
            return (value < b);
        }

        bool operator>=(const Ddata& b) const {
            return (value >= b.value);
        }
        bool operator>=(T b) const {
            return (value >= b);
        }

        bool operator<=(const Ddata& a) const {
            return (value <= a.value);
        }
        bool operator<=(T b) const {
            return (value <= b);
        }

        bool operator==(const Ddata& a) const {
            return (value == a.value);
        }
        bool operator==(T b) const {
            return (value == b);
        }

        Ddata<T, n> operator+(T rhs) const {
            auto ddata = *this;
            ddata.value += rhs;
            return ddata;
        }
        Ddata<T, n> operator+(const Ddata<T, n> &rhs) const {
            auto lhs = *this;
            lhs.value += rhs.value;
            for (size_t i = 0; i < n; i++)
                lhs.derivatives[i] += rhs.derivatives[i];
            return lhs;
        }

        Ddata<T, n> operator-(T rhs) const {
            auto ddata = *this;
            ddata.value -= rhs;
            return ddata;
        }
        Ddata<T, n> operator-(const Ddata<T, n> &rhs) const {
            auto lhs = *this;
            lhs.value -= rhs.value;
            for (size_t i = 0; i < n; i++)
                lhs.derivatives[i] -= rhs.derivatives[i];
            return lhs;
        }

        Ddata<T, n> operator*(T rhs) const {
            auto ddata = *this;
            ddata.value *= rhs;
            for (size_t i = 0; i < n; i++)
                ddata.derivatives[i] *= rhs;
            return ddata;
        }
        Ddata<T, n> operator*(const Ddata<T, n> &rhs) const {
            auto lhs = *this;
            for (size_t i = 0; i < n; i++)
                lhs.derivatives[i] = lhs.value * rhs.derivatives[i] + lhs.derivatives[i] * rhs.value;
            lhs.value *= rhs.value;
            return lhs;
        }

        Ddata<T, n> operator/(T rhs) const {
            auto ddata = *this;
            T denom = 1.0 / rhs / rhs;
            for (size_t i = 0; i < n; i++)
                ddata.derivatives[i] = (ddata.derivatives[i] * rhs) * denom;
            ddata.value = ddata.value / rhs;
            return ddata;
        }
        Ddata<T, n> operator/(const Ddata<T, n> &rhs) const {
            auto lhs = *this;
            T denom = 1.0 / (rhs.value * rhs.value);
            for (size_t i = 0; i < n; i++)
                lhs.derivatives[i] = (lhs.derivatives[i] * rhs.value - rhs.derivatives[i] * lhs.value) * denom;
            lhs.value /= rhs.value;
            return lhs;
        }

        friend std::ostream& operator<<(std::ostream& os, const Ddata& ddt) {  
            os << ddt.value << ": ";
            for(size_t i = 0; i < n; i++)
                os << ddt.derivatives[i] << " ";
            return os;  
        }  

        void print() {
            printf("f: %e\n", value);
            printf("d:");
            for (T d_value : derivatives)
                printf(" %e", d_value);
            printf("\n");
        }
    };

    template<typename T, size_t n>
    std::array<Ddata<T, n>, n>
    Identity(const std::array<T, n>& func) {
        std::array<Ddata<T, n>, n> ddt;
        for (size_t i = 0; i < n; i++) {
            ddt[i].value = func[i];
            ddt[i].derivatives.fill(0.0);
            ddt[i].derivatives[i] = 1.0;
        }
        return ddt;
    }

    template<typename T, size_t NumEqns, size_t NumVars>
    std::array<T, NumEqns*NumVars> ExtractBlockJacobian(std::array<Ddata<T, NumVars>, NumEqns> ddt) {
        std::array<T, NumEqns*NumVars> block_jacobian;
        for (size_t i = 0; i < NumEqns; ++i)
            for (size_t j = 0; j < NumVars; ++j)
                block_jacobian[i*NumVars + j] = ddt[i].derivatives[j];
        return block_jacobian;
    }
}

template<typename T, typename T2, size_t n>
inline Linearize::Ddata<T, n> operator+(T2 scale, const Linearize::Ddata<T, n>& rhs){
    return Linearize::Ddata<T, n>(T(rhs.value + scale), rhs.derivatives);
}

template<typename T, typename T2, size_t n>
inline Linearize::Ddata<T, n> operator-(T2 scale, const Linearize::Ddata<T, n> &rhs) {
    auto d = rhs.derivatives;
    for (T& value : d)
        value = -value;
    return Linearize::Ddata<T, n>(T(scale - rhs.value), d);
}

//FIXME: manually stamped these out due to issue with Tensor.h
template<size_t n>
inline Linearize::Ddata<double, n> operator*(double scale, const Linearize::Ddata<double, n> &rhs) {
    auto local_copy = rhs.derivatives;
    for (auto &entry:local_copy)
        entry *= scale;
    return Linearize::Ddata<double, n>(rhs.value * scale, local_copy);
}

//FIXME: manually stamped these out due to issue with Tensor.h
template<size_t n>
inline Linearize::Ddata<std::complex<double>, n>
operator*(const std::complex<double> &scale, const Linearize::Ddata<std::complex<double>, n> &rhs) {
    auto local_copy = rhs.derivatives;
    for (auto &entry:local_copy)
        entry *= scale;
    return Linearize::Ddata<std::complex<double>, n>(rhs.value * scale, local_copy);
}

template<typename T, typename T2, size_t n>
inline Linearize::Ddata<T, n> operator/(T2 scale, const Linearize::Ddata<T, n> &rhs) {
    auto local_copy = rhs.derivatives;
    T denom = 1.0 / rhs.value / rhs.value;
    for (auto &entry:local_copy)
        entry *= -scale * denom;
    return Linearize::Ddata<T, n>(T(scale / rhs.value), local_copy);
}

template<typename T, size_t n>
inline Linearize::Ddata<T, n> max(const Linearize::Ddata<T, n>& a, const T& b) {
    if (a.value > b)
        return a;
    else
        return Linearize::Ddata<T, n>(b);
}
template<typename T, size_t n>
inline Linearize::Ddata<T, n> max(const Linearize::Ddata<T, n>& a, const Linearize::Ddata<T, n>& b) {
    if ( a.value > b.value)
        return a;
    else
        return b;
}
template<typename T, size_t n>
inline Linearize::Ddata<T, n> max(T a, const Linearize::Ddata<T, n>& b) {
    if ( a > b.value)
        return Linearize::Ddata<T, n>(a);
    else
        return b;
}
template<typename T, size_t n>
inline Linearize::Ddata<T, n> min(const Linearize::Ddata<T, n>& a, const Linearize::Ddata<T, n>& b) {
    if ( a.value < b.value)
        return a;
    else
        return b;
}
template<typename T, typename T2, size_t n>
inline Linearize::Ddata<T, n> min(T2 a, const Linearize::Ddata<T, n>& b) {
    if ( a < b.value)
        return Linearize::Ddata<T, n>(T(a));
    else
        return b;
}
template<typename T, typename T2, size_t n>
inline Linearize::Ddata<T, n> min(const Linearize::Ddata<T, n>& a, T2 b) {
    if (a.value < b)
        return a;
    else
        return Linearize::Ddata<T, n>(T(b));
}
template<typename T, size_t n>
inline Linearize::Ddata<T, n> fabs(const Linearize::Ddata<T, n>& a) {
    if ( a.value < 0.0 )
        return -a;
    else
        return a;
}
template<typename T, size_t n>
inline Linearize::Ddata<T, n> sign(const Linearize::Ddata<T, n>& a, const Linearize::Ddata<T, n>& b) {
    if ( b.value >= 0.0 )
        return fabs(a);
    else
        return -fabs(a);
}
template<typename T, typename T2, size_t n>
inline Linearize::Ddata<T, n> pow(const Linearize::Ddata<T, n>& rhs, T2 r_in) {
    auto r = T(r_in);
    auto f = pow(rhs.value,r);
    std::array<T, n> d;
    for (size_t i = 0; i < n; i++)
        d[i] = r * pow(rhs.value,(r-1.0)) * rhs.derivatives[i];
    return Linearize::Ddata<T, n>(f, d);
}
template<typename T, size_t n>
inline Linearize::Ddata<T, n> pow(const Linearize::Ddata<T, n>& rhs, const Linearize::Ddata<T, n>& r) {
    T f = pow(rhs.value,r.value);
    std::array<T, n> d;
    for (size_t i = 0; i < n; i++)
        d[i] = r.value * pow(rhs.value, (r.value - 1.0)) * rhs.derivatives[i]
               + std::log(rhs.value) * pow(rhs.value, r.value) * r.derivatives[i];
    return Linearize::Ddata<T, n>(f, d);
}
template<typename T, size_t n>
inline Linearize::Ddata<T, n> exp(const Linearize::Ddata<T, n>& rhs) {
    T f = std::exp(rhs.value);
    std::array<T, n> d;
    for (size_t i = 0; i < n; i++)
        d[i] = f * rhs.derivatives[i];
    return Linearize::Ddata<T, n>(f, d);
}
template<typename T, size_t n>
inline Linearize::Ddata<T, n> log(const Linearize::Ddata<T, n>& rhs) {
    T f = std::log(rhs.value);
    std::array<T, n> d;
    for (size_t i = 0; i < n; i++) {
        d[i] = (1.0 / rhs.value) * rhs.derivatives[i];
    }
    return Linearize::Ddata<T, n>(f, d);
}
template<typename T, size_t n>
inline Linearize::Ddata<T, n> tanh(const Linearize::Ddata<T, n>& rhs) {
    T f = std::tanh(rhs.value);
    std::array<T, n> d;
    for (size_t i = 0; i < n; i++)
        d[i] = (1.0 - pow(f, 2)) * rhs.derivatives[i];
    return Linearize::Ddata<T, n>(f, d);
}
template<typename T, size_t n>
inline Linearize::Ddata<T, n> sin(const Linearize::Ddata<T, n>& rhs) {
    T f = std::sin(rhs.value);
    std::array<T, n> d;
    for (size_t i = 0; i < n; i++)
        d[i] = std::cos(rhs.value) * rhs.derivatives[i];
    return Linearize::Ddata<T, n>(f, d);
}
template<typename T, size_t n>
inline Linearize::Ddata<T, n> cos(const Linearize::Ddata<T, n>& rhs) {
    T f = std::cos(rhs.value);
    std::array<T, n> d;
    for (size_t i = 0; i < n; i++)
        d[i] = -std::sin(rhs.value) * rhs.derivatives[i];
    return Linearize::Ddata<T, n>(f, d);
}
template<typename T, size_t n>
inline Linearize::Ddata<T, n> asin(const Linearize::Ddata<T, n>& rhs) {
    T f = std::asin(rhs.value);
    std::array<T, n> d;
    for (size_t i = 0; i < n; i++)
        d[i] = (1.0 / std::sqrt(1.0 - pow(rhs.value, 2))) * rhs.derivatives[i];
    return Linearize::Ddata<T, n>(f, d);
}
template<typename T, size_t n>
inline Linearize::Ddata<T, n> acos(const Linearize::Ddata<T, n>& rhs) {
    T f = std::acos(rhs.value);
    std::array<T, n> d;
    for (size_t i = 0; i < n; i++)
        d[i] = -(1.0 / std::sqrt(1.0 - pow(rhs.value, 2))) * rhs.derivatives[i];
    return Linearize::Ddata<T, n>(f, d);
}
template<typename T, size_t n>
inline Linearize::Ddata<T, n> sqrt(const Linearize::Ddata<T, n>& rhs) {
    T f = std::sqrt(rhs.value);
    std::array<T, n> d;
    if (f > 0.0) {
        for (size_t i = 0; i < n; i++)
            d[i] = (0.5 / f) * rhs.derivatives[i];
    } else {
        d.fill(0.0);
    }
    return Linearize::Ddata<T, n>(f, d);
}

template <typename T, size_t n>
inline bool operator>( const T& a, const Linearize::Ddata<T, n>& b ) {
    return ( a > b.value );
}
template <typename T, size_t n>
inline bool operator<( const T& a, const Linearize::Ddata<T, n>& b ) {
    return ( a < b.value );
}
template <typename T, size_t n>
inline bool operator>=( const T& a, const Linearize::Ddata<T, n>& b ) {
    return ( a >= b.value );
}
template <typename T, size_t n>
inline bool operator<=( const T& a, const Linearize::Ddata<T, n>& b ) {
    return ( a <= b.value );
}
template <typename T, size_t n>
inline bool operator==( const T& a, const Linearize::Ddata<T, n>& b ) {
    return ( a == b.value );
}

