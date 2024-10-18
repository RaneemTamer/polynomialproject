#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <vector>
#include <complex>
#include <iostream>

using namespace std;

class Polynomial {
private:
    vector<double> coeffs; // Coefficients remain real (double)

public:
    // Constructors
    Polynomial();
    Polynomial(const vector<double>& coefficients);
    Polynomial(const Polynomial& other);
    ~Polynomial();

    // Operators
    Polynomial& operator=(const Polynomial& other);
    Polynomial operator+(const Polynomial& other) const;
    Polynomial operator-(const Polynomial& other) const;
    Polynomial operator*(const Polynomial& other) const;

    // Output operator
    friend ostream& operator<<(ostream& out, const Polynomial& poly);

    // Polynomial operations
    int degree() const;
    double evaluate(double x) const;
    Polynomial derivative() const;

    // Root finding using Newton's method with complex roots
    complex<double> getRoot(complex<double> guess, double tolerance, int maxIter);

    // Setters and Getters
    void setCoefficients(const vector<double>& coefficients);
    double getCoefficient(int degree) const;
};

#endif // POLYNOMIAL_H
