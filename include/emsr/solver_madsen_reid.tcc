/**
 * @file solver_madsen_reid.tcc Class outline definitions for the Madsen-Reid solver.
 */

/**
 * @def  SOLVER_MADSEN_REID_TCC
 *
 * @brief  A guard for the _BairstowSolver class header.
 */
#ifndef SOLVER_MADSEN_REID_TCC
#define SOLVER_MADSEN_REID_TCC 1

  /**
   * Root search...
   */
  template<typename Real>
    void
    SolverMadsenReid<Real>::
    solve(std::vector<std::complex<Real>>& a1, int m,
          std::vector<std::complex<Real>>& root, std::vector<std::complex<Real>>& a)
    {
        using Cmplx = std::complex<Real>;

        Cmplx z0, f0z, z, dz, f1z, fz, w, fw, dzk;
        Real f0, ff, f, fa, fmin, f2;
        bool stage1, div2;

        Real r0, u0, r, r1;

        const Real SSMALL = std::sqrt(SMALL);
        const Real ALOGB = std::log(BASE);

        const auto THETA = std::atan(Real{3} / Real{4});
        const auto PHASE = std::polar(Real{1}, -THETA);

        int n = m;

        // Store original polynomial in a and in root.
        int j = m;
        a[j] = a1[0];
        for (int i = 0; i < m; ++i)
        {
            root[i] = a1[i];
            a[i] = a1[j];
            --j;
        }

        // Test for zeros at infinity.
        while (norm_l1(a[0]) <= 0 && n > 0)
        {
            for (int i = 0; i < n; ++i)
            {
                a[i] = a[i + 1];
            }
            root[n] = BIG;
            --n;
        }

        while (n >= 0)
        {
            if (n == 1)
            {
                z = -a[1] / a[0];
                a1[n-1] = root[n-1];
                root[n-1] = z;
                return;
            }

            // Scale the coefficients.
            auto u1 = Real {0};
            auto u2 = BIG;
            for(int k = 0; k <= n; ++k)
            {
                auto u = norm_l1(a[k]);
                if (u <= Real{0})
                    continue;
                if (u > u1)
                    u1 = u;
                if (u < u2)
                    u2 = u;
            }
            auto u = std::sqrt(u1) * std::sqrt(u2);
            int i = -std::log(u) / ALOGB; // ilogb?
            u = std::pow(BASE, Real(i));
            for (int k = 0; k < n; ++k)
            {
                a[k] = u * a[k];
                a1[k] = a[k] * Real(n - k);
            }
            a[n] = u * a[n];

            // Test for zeros at (0, 0)
            z = Cmplx(0, 0);
            if (norm_l1(a[n]) <= SSMALL)
            {
                push_root(a1, root, n, z);
                continue;
            }
            z0 = Cmplx(0, 0);
            f0 = norm_l2(a[n]);
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
            for (int k = 0; k < n; ++k)
            {
                u = norm_l2(a[k]);
                if (u == Real{0})
                    continue;
                u = std::log(u0 / u) / Real(2 *(n - k));
                if (u < t)
                    t = u;
            }
            t = std::exp(t);
            f0z = a[n - 1];
            z = Cmplx(1, 0);
            if (norm_l1(f0z) > Real{0})
            {
                z = -a[n] / a[n - 1];
            }
            u = 0.5 * t / norm_l1(z);
            z = u * z;
            dz = z;
            f =  eval(z, fz, n + 1, a);
            r0 = 0.5 * t;

        _120:
            // Calculate tentative step.
            u = eval(z, f1z, n, a1);
            if (u == Real{0})
            {
                dz *= Real{3} * PHASE;
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
                    dz *= Real{3} * PHASE * (r0 / r);
                }
            }

            f0z = f1z;

        _160:
            // Find the next point in the iteration.
            // This is where iteration starts if the previous one was unsuccessful.
            z0 = z;
            f0 = f;
            dzk = dz;
            z = z0 + dz;
            // If either part of z is small replace by zero to avoid underflows.
            if (std::abs(std::real(z)) < EPS * std::abs(std::imag(z)))
                z = Cmplx(Real{0}, std::imag(z));
            if (std::abs(std::imag(z)) < EPS * std::abs(std::real(z)))
                z = Cmplx(std::real(z), Real{0});
            w = z;
            f = eval(z, fz, n + 1, a);
            ff = f;
            if (stage1)
            {
                int j = 1;
                div2 = f >= f0;

                do
                {
                    if (div2)
                    {
                        dz *= Real{0.5L};
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
                        dz *= PHASE;
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

            dz = -Real{0.5} * PHASE * dzk;
            stage1 = true;
            goto _160;
        }

        return;
    }

#endif // SOLVER_MADSEN_REID_TCC
