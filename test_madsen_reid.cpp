/*
 g++ -g -o test_madsen_reid test_madsen_reid.cpp
 ./test_madsen_reid < test/input/complex_solver1.in > out2
 ./test_madsen_reid < test/input/complex_solver2.in >> out2

 ./test_madsen_reid < test/input/test_solver1.in >> out2
 ./test_madsen_reid < test/input/test_solver2.in >> out2
 ./test_madsen_reid < test/input/test_solver3.in >> out2
 ./test_madsen_reid < test/input/test_solver4.in >> out2
 ./test_madsen_reid < test/input/test_solver5.in >> out2
 ./test_madsen_reid < test/input/test_solver6.in >> out2
 ./test_madsen_reid < test/input/test_solver7.in >> out2
 ./test_madsen_reid < test/input/test_solver8.in >> out2
*/

#include <limits>
#include <complex>
#include <vector>
#include <iostream>

using Real = double;
using Cmplx = std::complex<Real>;

/**
 * Return the L1 sum of absolute values or Manhattan metric of a complex number.
 */
Real
norm_l1(Cmplx z)
{
    return std::abs(std::real(z)) + std::abs(std::imag(z));
}

/**
 * Return the L2 modulus of the complex number (this is std::norm).
 */
Real
norm_l2(Cmplx z)
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

    return norm_l2(p);
}

/**
 * Store the root.
 */
void
push_root(std::vector<Cmplx>& a1, std::vector<Cmplx>& root, int& n, Cmplx z)
{
    a1[n] = root[n];
    root[n] = z;
    --n;
}

/**
 * Store the root.
 */
void
deflate(std::vector<Cmplx>& a, int n, Cmplx z)
{
    for (int k = 2; k <= n; ++k)
    {
        a[k] = a[k - 1] * z + a[k];
    }
}

//
// Root search...
//
void
solve(std::vector<Cmplx>& a1, int m, std::vector<Cmplx>& root, std::vector<Cmplx>& a, int mp1)
{
    Cmplx z0, f0z, z, dz, f1z, fz, w, fw, dzk;
    Real f0, ff, f, fa, fmin, f2;
    bool stage1, div2;

    Real r0, u0, r, r1;

    const auto theta = std::atan(Real{3} / Real{4});
    const auto phase = std::polar(Real{1}, theta);

    const Real DIGITS = std::numeric_limits<Real>::max_digits10;
    const Real BIG = std::numeric_limits<Real>::max(); // Overflow limit
    const Real SMALL = std::numeric_limits<Real>::min(); // Underflow limit.
    const Real BASE = std::numeric_limits<Real>::radix; // 16 -> 2
    const Real EPS = std::numeric_limits<Real>::epsilon();

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
    while (norm_l1(a[1]) <= 0 && n > 0)
    {
        for (int i = 1; i <= n; ++i)
        {
            a[i] = a[i + 1];
        }
        root[n] = BIG;
        --n;
    }

    while (n >= 1)
    {
        if (n == 1)
        {
            z = -a[2] / a[1];
            a1[n] = root[n];
            root[n] = z;
            return;
        }

        int n1 = n + 1;

        // Scale the coefficients.
        auto u1 = Real {0};
        auto u2 = BIG;
        for(int k = 1; k <= n1; ++k)
        {
            u = norm_l1(a[k]);
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
        auto u = std::sqrt(u1) * std::sqrt(u2);
        int i = -std::log(u) / ALOGB; // ilogb?
        u = std::pow(BASE, Real(i));
        for (int k = 1; k <= n; ++k)
        {
            a[k] = u * a[k];
            a1[k] = a[k] * Real(n1 - k);
        }
        a[n1] = u * a[n1];

        // Test for zeros at (0, 0)
        z = Cmplx(0, 0);
        if (norm_l1(a[n1]) <= SSMALL)
        {
            push_root(a1, root, n, z);
            continue;
        }
        z0 = Cmplx(0, 0);
        f0 = norm_l2(a[n1]);
        fmin = f0 * std::pow(Real(n) * DIGITS * EPS, 2);

        // z is the current point
        // f = |f(z)|^2
        // z0 is the last point
        // f0z = f'(z0)
        // f0 = |f(z0)|^2
        // r0 = 3|z - z0|
        // dz is the last tentative step if the last step was successful or is the required next step.
        ff = f0;
        u0 = f0;
        auto t = BIG;
        for (int k = 1; k <= n; ++k)
        {
            u = norm_l2(a[k]);
            if (u == Real{0}) {
                continue;
            }
            u = std::log(u0 / u) / Real(2 *(n1 - k));
            if (u < t) {
                t = u;
            }
        }
        t = std::exp(t);
        f0z = a[n];
        z = Cmplx(1, 0);
        if (norm_l1(f0z) > Real{0})
        {
            z = -a[n1] / a[n];
        }
        u = 0.5 * t / norm_l1(z);
        z = u * z;
        dz = z;
        f =  eval(z, fz, n + 1, a);
        r0 = 0.5 * t;

        // Calculate tentative step.
    _120:
        u = eval(z, f1z, n, a1);
        if (u == 0.0)
        {
            dz *= Real{3} * phase;
            stage1 = true;
        }
        else
        {
            dz = -fz / f1z;
            f2 = norm_l2(f0z - f1z) / norm_l2(z0 - z);
            stage1 = (f * f2 / u > 0.25 * u) || (f != ff);
            r = norm_l1(dz);
            if (r > Real{3} * r0)
            {
                dz *= Real{3} * phase * r0 / r;
            }
        }

        f0z = f1z;

        // Find the next point in the iteration.
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
        if (stage1)
        {
            j = 1;
            div2 = f >= f0;

            do
            {
                if (div2)
                {
                    dz = 0.5 * dz;
                    w = z0 + dz;
                }
                else
                {
                    w += dz;
                }

                fa = eval(w, fw, n + 1, a);
                if (fa >= f)
                {
                    break;
                }
                f = fa;
                fz = fw;
                z = w;
                ++j;
                if (div2 && j == 3)
                {
                    dz *= phase;
                    z = z0 + dz;
                    f = eval(z, fz, n + 1, a);
                    break;
                }
            }
            while (j <= n);
        }

        r0 = norm_l1(z0 - z);

        // Convergence test.
        if (f >= f0)
        {
            z = z0;
        }

        r1 = norm_l1(z);
        if (r0 < EPS * r1)
        {
            deflate(a, n, z);
            push_root(a1, root, n, z);
            continue;
        }

        if (f < f0)
        {
            goto _120;
        }

        f = f0;
        if (f < fmin)
        {
            deflate(a, n, z);
            push_root(a1, root, n, z);
            continue;
        }

        dz = -Real{0.5} * phase * dzk;
        stage1 = true;
        goto _160;
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
