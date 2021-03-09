#pragma once

#include "mycomplex.hpp"
#include "signalgen.h"

void test_initialization () {
    MyComplex<double> num (1, 2);
    MyComplex<std::string> nonnum ("", "");
    bool success = num && !nonnum;
    if (success)
        std::cout << "Initialization test succeeded :)\n";
    else
        std::cerr << "Initialization test failed :(\n";
}

void test_addition () {
    MyComplex<int> num1 (1, 2);
    MyComplex<int> num2 (2, 3);
    auto sum = num1 + num2;
    bool success = sum.real () == 3 && sum.imag () == 5;
    if (success)
        std::cout << "Addition test succeeded :)\n";
    else
        std::cerr << "Addition test failed :(\n";
}

void test_subtraction () {
    MyComplex<int> num1 (3, 5);
    MyComplex<int> num2 (1, 2);
    auto diff = num1 - num2;
    bool success = diff.real () == 2 && diff.imag () == 3;
    if (success)
        std::cout << "Subtraction test succeeded :)\n";
    else
        std::cerr << "Subtraction test failed :(\n";
}

void test_multiplication () {
    // (+, +) * (+, +)
    MyComplex<int> num1 (3, 2);
    MyComplex<int> num2 (1, 7);
    auto product = num1 * num2;
    bool success = product.real () == -11 && product.imag () == 23;

    // (+, +) * (-, +)
    num1 = MyComplex (3, 2);
    num2 = MyComplex (-1, 7);
    product = num1 * num2;
    success = success && product.real () == -17 && product.imag () == 19;

    // (+, +) * (+, -)
    num1 = MyComplex (3, 2);
    num2 = MyComplex (1, -7);
    product = num1 * num2;
    success = success && product.real () == 17 && product.imag () == -19;

    // (+, +) * (-, -)
    num1 = MyComplex (3, 2);
    num2 = MyComplex (-1, -7);
    product = num1 * num2;
    success = success && product.real () == 11 && product.imag () == -23;

    // (+, -) * (+, +)
    num1 = MyComplex (3, -2);
    num2 = MyComplex (1, 7);
    product = num1 * num2;
    success = success && product.real () == 17 && product.imag () == 19;

    // (+, -) * (-, +)
    num1 = MyComplex (3, -2);
    num2 = MyComplex (-1, 7);
    product = num1 * num2;
    success = success && product.real () == 11 && product.imag () == 23;

    // (+, -) * (+, -)
    num1 = MyComplex (3, -2);
    num2 = MyComplex (1, -7);
    product = num1 * num2;
    success = success && product.real () == -11 && product.imag () == -23;

    // (+, -) * (-, -)
    num1 = MyComplex (3, -2);
    num2 = MyComplex (-1, -7);
    product = num1 * num2;
    success = success && product.real () == -17 && product.imag () == -19;

    // (-, +) * (+, +)
    num1 = MyComplex (-3, 2);
    num2 = MyComplex (1, 7);
    product = num1 * num2;
    success = success && product.real () == -17 && product.imag () == -19;

    // (-, +) * (-, +)
    num1 = MyComplex (-3, 2);
    num2 = MyComplex (-1, 7);
    product = num1 * num2;
    success = success && product.real () == -11 && product.imag () == -23;

    // (-, +) * (+, -)
    num1 = MyComplex (-3, 2);
    num2 = MyComplex (1, -7);
    product = num1 * num2;
    success = success && product.real () == 11 && product.imag () == 23;

    // (-, +) * (-, -)
    num1 = MyComplex (-3, 2);
    num2 = MyComplex (-1, -7);
    product = num1 * num2;
    success = success && product.real () == 17 && product.imag () == 19;

    // (-, -) * (+, +)
    num1 = MyComplex (-3, -2);
    num2 = MyComplex (1, 7);
    product = num1 * num2;
    success = success && product.real () == 11 && product.imag () == -23;

    // (-, -) * (-, +)
    num1 = MyComplex (-3, -2);
    num2 = MyComplex (-1, 7);
    product = num1 * num2;
    success = success && product.real () == 17 && product.imag () == -19;

    // (-, -) * (+, -)
    num1 = MyComplex (-3, -2);
    num2 = MyComplex (1, -7);
    product = num1 * num2;
    success = success && product.real () == -17 && product.imag () == 19;

    // (-, -) * (-, -)
    num1 = MyComplex (-3, -2);
    num2 = MyComplex (-1, -7);
    product = num1 * num2;
    success = success && product.real () == -11 && product.imag () == 23;

    if (success)
        std::cout << "Multiplication test succeeded :)\n";
    else
        std::cerr << "Multiplication test failed :(\n";
}

void test_conjugate () {
    MyComplex<double> num (2.5, 3.7);
    auto conj = num.conjugate ();
    bool success = conj.real () == num.real () && conj.imag () == -num.imag ();

    if (success)
        std::cout << "Conjugate test succeeded :)\n";
    else
        std::cerr << "Conjugate test failed :(\n";
}

void test_division () {
    MyComplex<double> num1 (2.5, 3.7);
    MyComplex<double> num2 (4.5, 5.5);
    MyComplex<double> remainder = num1 / num2;
    bool success = remainder.real () >= 0.6257 && remainder.imag () >= 0.0574;
    if (success)
        std::cout << "Division test succeeded :)\n";
    else
        std::cerr << "Division test failed :(\n";
}

void test_exp () {
    for (size_t i = 0; i < 100000; i++) {
#define RANDOM() range_rand (-1000, 1000)
        double nume = RANDOM () + 1.1;
        double denom = RANDOM () + 2.1;
        double real = nume / denom;

        nume = RANDOM () + 1.1;
        denom = RANDOM () + 2.1;
        double imag = nume / denom;

        MyComplex<double> mine (real, imag);
        std::complex<double> his (real, imag);

        auto num1 = mine.exp();
        auto num2 = exp(his);

        if (fabs(num2.real()) == INFINITY || fabs(num2.imag()) == INFINITY)
            continue;

        double tolerance = 1e-10;

        double real_diff = fabs(fabs(num1.real()) - fabs(num2.real()));
        double imag_diff = fabs(fabs(num1.imag()) - fabs(num2.imag()));

        bool success = real_diff < tolerance && imag_diff < tolerance;

        if (!success) {
            std::cerr << "Exp test failed :(\n";
            std::cerr << "Diff: (" << real_diff << ", " << imag_diff << ")\n";
            std::cerr << "Mine: (" << num1.real() << ", " << num1.imag() << ")\n";
            std::cerr << "His: " << num2 << "\n";
            return;
        }
#undef RANDOM
    }

    std::cout << "Exp test succeeded :)\n";
}

void complex_test () {
    test_initialization ();
    test_addition ();
    test_subtraction ();
    test_multiplication ();
    test_conjugate ();
    test_division ();
    test_exp ();
}