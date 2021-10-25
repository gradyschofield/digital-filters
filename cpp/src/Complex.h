//
// Created by Grady Schofield on 10/17/21.
//

#ifndef CPP_COMPLEX_H
#define CPP_COMPLEX_H

#include<cmath>

#if 0

#include<complex> // this has a very slow multiplication operator

#else

template<typename T>
class complex{
    T re;
    T im;

public:
    complex()
        : re(0), im(0)
    {
    }

    complex(T re)
        : re(re), im(0)
    {
    }

    complex(T re, T im)
        : re(re), im(im)
    {
    }

    complex(complex<T> const & c)
        : re(c.real()), im(c.imag())
    {
    }

    T real() const {
        return re;
    }

    T imag() const {
        return im;
    }

    complex<T> operator-() const {
        return complex<T>{-re, -im};
    }

    complex<T> & operator+=(complex<T> const & c) {
        re += c.real();
        im += c.imag();
        return *this;
    }

    complex<T> & operator-=(complex<T> const & c) {
        re -= c.real();
        im -= c.imag();
        return *this;
    }

    complex<T> & operator*=(complex<T> const & c2) {
        T a = re;
        T b = im;
        T c = c2.real();
        T d = c2.imag();
        *this = complex<T>{a*c-b*d, a*d+b*c};
        return *this;
    }

    complex<T> & operator/=(complex<T> const & c2) {
        T c = c2.real();
        T d = c2.imag();
        *this *= complex<T>(c, -d);
        T inv = 1.0 / (c*c + d*d);
        re *= inv;
        im *= inv;
        return *this;
    }
};

template<typename T>
complex<T> conj(complex<T> const & c) {
    return complex<T>{c.real(), -c.imag()};
}

template<typename T>
T norm(complex<T> const & c) {
    return c.real()*c.real() + c.imag()*c.imag();
}

template<typename T>
T abs(complex<T> const & c) {
    return sqrt(norm(c));
}

template<typename T>
complex<T> exp(complex<T> const & c) {
    T r = exp(c.real());
    return complex<T>{r*cos(c.imag()), r*sin(c.imag())};
}


template<typename T>
complex<T> operator+(complex<T> const & c1, complex<T> const & c2) {
    return complex<T>{c1.real()+c2.real(), c1.imag()+c2.imag()};
}

template<typename T>
complex<T> operator-(complex<T> const & c1, complex<T> const & c2) {
    return complex<T>{c1.real()-c2.real(), c1.imag()-c2.imag()};
}

template<typename T>
complex<T> operator+(T r, complex<T> const & c) {
    return complex<T>{r+c.real(), c.imag()};
}

template<typename T>
complex<T> operator-(T r, complex<T> const & c) {
    return complex<T>{r-c.real(), -c.imag()};
}

template<typename T>
complex<T> operator+(complex<T> const & c, T r) {
    return complex<T>{r+c.real(), c.imag()};
}

template<typename T>
complex<T> operator-(complex<T> const & c, T r) {
    return complex<T>{c.real() - r, c.imag()};
}

template<typename T>
complex<T> operator*(complex<T> const & c1, complex<T> const & c2) {
    T a = c1.real();
    T b = c1.imag();
    T c = c2.real();
    T d = c2.imag();
    return complex<T>{a*c-b*d, a*d+b*c};
}

template<typename T>
complex<T> operator/(complex<T> const & c1, complex<T> const & c2) {
    complex<T> ret(c1);
    ret *= conj(c2);
    ret /= norm(c2);
    return ret;
}

template<typename T>
complex<T> operator*(complex<T> const & c, T r) {
    return complex<T>{c.real()*r, c.imag()*r};
}

template<typename T>
complex<T> operator*(T r, complex<T> const & c) {
    return complex<T>{c.real()*r, c.imag()*r};
}

template<typename T>
complex<T> operator/(complex<T> const & c, T r) {
    T inv = 1.0 / r;
    return complex<T>{c.real()*inv, c.imag()*inv};
}

template<typename T>
complex<T> operator/(T r, complex<T> const & c) {
    complex<T> ret{r, 0};
    ret /= c;
    return ret;
}

#endif

#endif //CPP_COMPLEX_H
