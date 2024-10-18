#include "polynomial.h"
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
using namespace std;

// Default constructor
Polynomial::Polynomial() : coeffs(1, 0.0) {}

// Constructor with coefficients vector
Polynomial::Polynomial(const vector<double>& coefficients) : coeffs(coefficients) {
    while (!coeffs.empty() && coeffs.back() == 0.0) {
        coeffs.pop_back(); // Remove trailing zeros
    }
}

// Copy constructor
Polynomial::Polynomial(const Polynomial& other) : coeffs(other.coeffs) {}

// Destructor
Polynomial::~Polynomial() {}

// Assignment operator
Polynomial& Polynomial::operator=(const Polynomial& other) {
    if (this != &other) {
        coeffs = other.coeffs;
    }
    return *this;
}

// Addition operator
Polynomial Polynomial::operator+(const Polynomial& other) const {
    vector<double> result(max(coeffs.size(), other.coeffs.size()), 0.0);
    for (size_t i = 0; i < result.size(); ++i) {
        if (i < coeffs.size()) result[i] += coeffs[i];
        if (i < other.coeffs.size()) result[i] += other.coeffs[i];
    }
    return Polynomial(result);
}

// Subtraction operator
Polynomial Polynomial::operator-(const Polynomial& other) const {
    vector<double> result(max(coeffs.size(), other.coeffs.size()), 0.0);
    for (size_t i = 0; i < result.size(); ++i) {
        if (i < coeffs.size()) result[i] += coeffs[i];
        if (i < other.coeffs.size()) result[i] -= other.coeffs[i];
    }
    return Polynomial(result);
}

// Multiplication operator
Polynomial Polynomial::operator*(const Polynomial& other) const {
    vector<double> result(coeffs.size() + other.coeffs.size() - 1, 0.0);
    for (size_t i = 0; i < coeffs.size(); ++i) {
        for (size_t j = 0; j < other.coeffs.size(); ++j) {
            result[i + j] += coeffs[i] * other.coeffs[j];
        }
    }
    return Polynomial(result);
}

// Output operator
ostream& operator<<(ostream& out, const Polynomial& poly) {
    bool first = true;
    for (int i = poly.coeffs.size() - 1; i >= 0; --i) {
        if (poly.coeffs[i] != 0.0) {
            if (!first) {
                out << (poly.coeffs[i] > 0 ? " + " : " - ");
            } else {
                if (poly.coeffs[i] < 0) out << "-";
                first = false;
            }
            if (abs(poly.coeffs[i]) != 1 || i == 0) out << abs(poly.coeffs[i]);
            if (i > 0) out << "x";
            if (i > 1) out << "^" << i;
        }
    }
    if (first) out << "0"; // if polynomial is zero
    return out;
}

// Degree of polynomial
int Polynomial::degree() const {
    return coeffs.size() - 1;
}

// Evaluate polynomial at real x
double Polynomial::evaluate(double x) const {
    double result = 0.0;
    for (int i = coeffs.size() - 1; i >= 0; --i) {
        result = result * x + coeffs[i];
    }
    return result;
}

// Derivative of polynomial
Polynomial Polynomial::derivative() const {
    if (coeffs.size() <= 1) return Polynomial(); // derivative of a constant is 0
    vector<double> result(coeffs.size() - 1, 0.0);
    for (size_t i = 1; i < coeffs.size(); ++i) {
        result[i - 1] = coeffs[i] * static_cast<double>(i);
    }
    return Polynomial(result);
}

// Newton's method to find complex roots
complex<double> Polynomial::getRoot(complex<double> guess, double tolerance, int maxIter) {
    Polynomial deriv = this->derivative();
    complex<double> x = guess;

    for (int i = 0; i < maxIter; ++i) {
        complex<double> fx(this->evaluate(x.real())); // Real part of evaluation for real coefficients
        complex<double> fpx(deriv.evaluate(x.real()));

        if (abs(fx) < tolerance) return x; // If the solution is within tolerance, return it
        if (abs(fpx) < tolerance) break;   // Avoid division by zero

        x = x - fx / fpx; // Newton's iteration
    }
    return x; // Return last approximation if max iterations reached
}

// Set coefficients
void Polynomial::setCoefficients(const vector<double>& coefficients) {
    coeffs = coefficients;
    while (!coeffs.empty() && coeffs.back() == 0.0) {
        coeffs.pop_back(); // Remove trailing zeros
    }
}

// Get coefficient of a specific degree
double Polynomial::getCoefficient(int degree) const {
    if (degree >= 0 && degree < coeffs.size()) {
        return coeffs[degree];
    }
    return 0.0; // Return 0 if degree is outside range
}
