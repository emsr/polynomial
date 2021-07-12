/*
 g++ -g -o test_madsen_reid test_madsen_reid.cpp
 ./test_madsen_reid < test/input/complex_solver1.in
 ./test_madsen_reid < test/input/complex_solver2.in

 ./test_madsen_reid < test/input/test_solver1.in
 ./test_madsen_reid < test/input/test_solver2.in
 ./test_madsen_reid < test/input/test_solver3.in
 ./test_madsen_reid < test/input/test_solver4.in
 ./test_madsen_reid < test/input/test_solver5.in
 ./test_madsen_reid < test/input/test_solver6.in
 ./test_madsen_reid < test/input/test_solver7.in
 ./test_madsen_reid < test/input/test_solver8.in
*/

#include <complex>
#include <vector>
#include <iostream>

using Real = double;
using Cmplx = std::complex<Real>;

/**
 * Return the L1 sum of absolute values or Manhattan metric of a complex number.
 */
Real
norml1(Cmplx z)
{
    return std::abs(std::real(z)) + std::abs(std::imag(z));
}

/**
 * Return the L2 modulus of the complex number (this is std::norm).
 */
Real
norml2(Cmplx z)
{
    return std::norm(z);
}

/**
 * Evaluate polynomial at z, set fz, return squared modulus.
 */
Real
eval(Cmplx z, Cmplx& fz, int n1, const std::vector<Cmplx>& a)
{
    auto n = n1 - 1;
    auto p = a[1];
    for (int i = 1; i <= n; ++i)
    {
        p = p * z + a[i + 1];
    }

    fz = p;

    return norml2(p);
}

//
// Root search...
//
void
solve(std::vector<Cmplx>& a1, int m, std::vector<Cmplx>& root, std::vector<Cmplx>& a, int mp1)
{

    //Cmplx a1(m+1), root(m), a(mp1)
    Cmplx z0, f0z, z, dz, f1z, fz, w, fw, dzk;
    Real f0, ff, f, fa, fmin, f2;
    bool stage1, div2;

    Real r0, t, u, u0, u1, u2, r, r1;

    int i, n1;

    const Real BIG = 1.0e+70; // Overflow limit
    const Real SMALL = 1.0e-70; // Underflow limit.
    //data BASE/16;
    const Real BASE = 2;
    const Real EPS = 2.3e-16;

    const Real SSMALL = std::sqrt(SMALL);
    const Real ALOGB = std::log(BASE);

    int n = m;

    // Store original polynomial in a and in root.
    int j = m + 1;
    a[j] = a1[1];
    for (int i = 1; i <= m; ++i)
    {
        root[i] = a1[i];
        a[i] = a1[j];
        --j;
    }

    // Test for zeros at infinity.
    while (norml1(a[1]) <= 0 && n > 0)
    {
        for (int i = 1; i <= n; ++i)
        {
            a[i] = a[i + 1];
        }
        root[n] = BIG;
        --n;
    }
    if (n <= 0)
        return;

 _40:
    if (n <= 1) {
        goto _260;
    }

    n1 = n + 1;

    // Scale the coefficients.
    u1 = 0.0;
    u2 = BIG;
    for(int k = 1; k <= n1; ++k)
    {
        u = norml1(a[k]);
        if (u <= 0.0) {
            continue;
        }
        if (u > u1) {
            u1 = u;
        }
        if (u < u2) {
            u2 = u;
        }
    }
    u = std::sqrt(u1) * std::sqrt(u2);
    i = -std::log(u) / ALOGB; // ilogb?
    u = std::pow(BASE, Real(i));
    for (int k = 1; k <= n; ++k)
    {
        a[k] = u * a[k];
        a1[k] = a[k] * Real(n1 - k);
    }
    a[n1] = u * a[n1];

    // Test for zeros at (0, 0)
    z = Cmplx(0, 0);
    if (norml1(a[n1]) <= SSMALL)
    {
        goto _290;
    }
    z0 = Cmplx(0, 0);
    f0 = norml2(a[n1]);
    fmin = f0 * std::pow(Real(n) * 2 * EPS, 2); // 16 -> 2

    // z is te current point
    // f = |f(z)|^2
    // z0 is the last point
    // f0z = f'(z0)
    // f0 = |f(z0)|^2
    // r0 = 3|z - z0|
    // dz is the last tentative step if the last step was successful or is the required next step.
    ff = f0;
    u0 = f0;
    t = BIG;
    for (int k = 1; k <= n; ++k)
    {
        u = norml2(a[k]);
        if (u == 0.0) {
            continue;
        }
        u = log(u0 / u) / Real(2 *(n1 - k));
        if (u < t) {
            t = u;
        }
    }
    t = std::exp(t);
    f0z = a[n];
    z = Cmplx(1, 0);
    if (norml1(f0z) > Real{0})
    {
        z = -a[n1] / a[n];
    }
    u = 0.5 * t / norml1(z);
    z = u * z;
    dz = z;
    f =  eval(z, fz, n + 1, a);
    r0 = 0.5 * t;

    // Calculate tentative step.
_120:
    u = eval(z, f1z, n, a1);
    if (u == 0.0) {
        goto _140;
    }
    dz = -fz / f1z;
    f2 = norml2(f0z - f1z) / norml2(z0 - z);
    stage1 = (f * f2 / u > 0.25 * u) || (f != ff);
    r = norml1(dz);
    if (r <= 3.0 * r0) {
        goto _150;
    }

    dz = Cmplx(1.8 * std::real(dz) - 2.4 * std::imag(dz) * r0 / r,
               2.4 * std::real(dz) + 1.8 * std::imag(dz) * r0 / r);
    goto _150;

_140:
    dz = Cmplx(1.8 * std::real(dz) - 2.4 * std::imag(dz),
               2.4 * std::real(dz) + 1.8 * std::imag(dz));
    stage1 = true;

_150:
    f0z = f1z;

    // Find the next point in te iteration.
    // This is where iteration starts if the previous one was unsuccessful.
_160:
    z0 = z;
    f0 = f;
    dzk = dz;
    z = z0 + dz;
    // If either part of z is small replace by zero to avoid underflows.
    if (std::abs(std::real(z)) < EPS * std::abs(std::imag(z))) {
        z = Cmplx(0, std::imag(z));
    }
    if (std::abs(std::imag(z)) < EPS * std::abs(std::real(z))) {
        z = Cmplx(std::real(z), 0.0);
    }
    w = z;
    f = eval(z, fz, n + 1, a);
    ff = f;
    if (!stage1) {
        goto _240;
    }

    // Beginning of stage 1 search.
    j = 1;
    div2 = f >= f0;
_180:
    if (div2) {
        dz = 0.5 * dz;
        w = z0 + dz;
    }else{
        w += dz;
    }

    fa = eval(w, fw, n + 1, a);
    if (fa >= f) {
        goto _240;
    }
    f = fa;
    fz = fw;
    z = w;
    ++j;
    if (div2 && j == 3) {
        goto _220;
    }
    if (j <= n) {
        goto _180;
    }
    goto _240;
_220:
    dz = Cmplx(0.6 * std::real(dz) - 0.8 * std::imag(dz),
               0.8 * std::real(dz) + 0.6 * std::imag(dz));
    z = z0 + dz;
    f = eval(z, fz, n + 1, a);

    // End of stage 1 search.

_240:
    r0 = norml1(z0 - z);

    // Convergence test.
    if (f < f0) {
        goto _250;
    }
    z = z0;
_250:
    r1 = norml1(z);
    if (r0 < EPS * r1) {
        goto _270;
    }
    if (f < f0) {
        goto _120;
    }
    f = f0;
    if (f < fmin) {
        goto _270;
    }
    dz = Cmplx(-0.3 * std::real(dzk) + 0.4 * std::imag(dzk),
               -0.4 * std::real(dzk) - 0.3 * std::imag(dzk));
    stage1 = true;
    goto _160;

_260:
    // Deal with n == 1 case.
    z = -a[2] / a[1];
    goto _290;

    // Deflate, store root, restore coefficient of original polynomial.
_270:
    for (int k = 2; k <= n; ++k)
    {
        a[k] = a[k - 1] * z + a[k];
    }
_290:
    a1[n] = root[n];
    root[n] = z;
    --n;
    if (n < 1) {
        return;
    }else if (n == 1) {
        goto _260;
    }else{
        goto _40;
    }

    return;
}



int
main()
{
    int m;
    std::vector<Cmplx> a1(21), a(21), root(20);

    std::cout << "Enter degree: ";
    std::cin >> m;
    if (m <= 0)
        return 0;

    a1.resize(m + 1);
    a.resize(m + 1);
    root.resize(m);

    for (int k = 1; k <= m + 1; ++k)
    {
        std::cout << "Enter coefficient " << k - 1 << ": ";
        std::cin >> a1[m + 2 - k];
    }
    
    solve(a1, m, root, a, m+1);

    std::cout << '\n';
    for (int k = 1; k <= m; ++k)
    {
        std::cout << "Zero " << k - 1 << ": " << root[k] << '\n';
    }
}
