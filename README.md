# polynomial

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/7f6bcfc8326a4dde8ff0d56069e94ed3)](https://app.codacy.com/app/emsr/polynomial?utm_source=github.com&utm_medium=referral&utm_content=emsr/polynomial&utm_campaign=Badge_Grade_Dashboard)
[![Build Status](https://travis-ci.org/emsr/polynomial.svg?branch=master)](https://travis-ci.org/emsr/polynomial)

<em>Note: I originally envisioned providing this library as an extension to the gnu C++ standard library.
I've given up on that idea. The standard C++ library is difficult enough to keep up with as it is.
I've changed the namespace from __gnu_cxx to emsr (if nlohmann can do it so can I :-)). It stands for <b>Excellent Math System - Reborn</b>. Really!</em>

This project contains C++ polynomial classes and related algorithms.
The primary polynomial class is a dense polynomial array with all of your favorite
polynomial arithmetic operator overloads.
It is designed to work with C and Fortran algorithms that have a size, pointer interface.

I am working on a sparse polynomial that is a set of (possibly multivariate) monomials.

This library has implementations of several root finders including Jenkins-Traub (real and complex), Madsen-Reid and Bairstow, and Laguerre and quadratic factorization steppers.
