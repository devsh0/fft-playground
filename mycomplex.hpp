#pragma once
#include <exception>
#include <iostream>
#include <cmath>
#include <type_traits>

template<typename T>
class MyComplex {
private:
    T m_real;
    T m_imag;

public:
    MyComplex (T real, T imag) : m_real (real),
                                 m_imag (imag) {}
    MyComplex () = delete;

    explicit constexpr operator bool () {
        return std::is_arithmetic_v<T>;
    }

    bool constexpr operator! () {
        return !std::is_arithmetic_v<T>;
    }

    MyComplex operator+ (const MyComplex& other) const {
        return {m_real + other.m_real, m_imag + other.m_imag};
    }

    MyComplex operator- (const MyComplex& other) const {
        return {m_real - other.m_real, m_imag - other.m_imag};
    }

    MyComplex operator* (T scalar) const {
        return {m_real * scalar, m_imag * scalar};
    }

    MyComplex operator* (const MyComplex& other) const {
        MyComplex outer = other * m_real;
        MyComplex inner = {-(m_imag * other.m_imag), m_imag * other.m_real};
        return {outer + inner};
    }

    MyComplex operator/ (const MyComplex& other) const {
        auto conj = other.conjugate ();
        T denom = (other.real () * other.real ()) + (other.imag () * other.imag ());
        auto nume = *this * conj;
        return {nume.m_real / denom, nume.m_imag / denom};
    }

    MyComplex conjugate () const {
        return {m_real, -m_imag};
    }

    MyComplex exp () const {
        if (m_real == 0 && m_imag == 0)
            return {1, 0};
        double factor = std::exp(m_real);
        int sign = m_imag >= 0 ? 1 : -1;
        return MyComplex{cos(m_imag), sign * sin(fabs(m_imag))} * factor;
    }

    T real () const {
        return m_real;
    }

    T imag () const {
        return m_imag;
    }
};
