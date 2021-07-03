/**
 * @def  SOLVER_JENKINS_TRAUB_COMPLEX_TCC
 *
 * @brief  A guard for the _JenkinsTraubSolver class header.
 */
#ifndef SOLVER_JENKINS_TRAUB_COMPLEX_TCC
#define SOLVER_JENKINS_TRAUB_COMPLEX_TCC 1

#include <vector>
#include <complex>
#include <limits>

namespace __gnu_cxx
{

template<typename Real>
class _JenkinsTraubSolver<std::complex<Real>>
{

public:

    using Cmplx = std::complex<Real>;

    _JenkinsTraubSolver(const std::vector<Cmplx>& op)
    : p(op)
    {
        if (p.size() == 0)
            throw std::domain_error("Polynomial degree must be at least 1.");

        // Algorithm fails if the leading-order coefficient is zero
        if (op[0] == ZERO)
            throw std::domain_error("Leading-order coefficient must be nonzero.");;

        degree = p.size() - 1;

        h.resize(degree + 1);
        qp.resize(degree + 1);
        qh.resize(degree + 1);
        sh.resize(degree + 1);
    }

    std::vector<Cmplx>
    solve()
    {
        std::vector<Cmplx> zero;
        solve(zero);
        return zero;
    }

private:

    Cmplx s;
    Cmplx t;
    Cmplx pv;

    Real add_rel_err, mul_rel_err, epsilon, infin;

    int degree;

    std::vector<Cmplx> p;
    std::vector<Cmplx> h;
    std::vector<Cmplx> qp;
    std::vector<Cmplx> qh;
    std::vector<Cmplx> sh;

    inline static constexpr Cmplx ZERO{0.0, 0.0};

    static constexpr auto s_sqrt2 = Real{1.4142'13562'37309'50488'01688'72420'96980'78569e+0L};
    static constexpr auto s_pi = Real{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L};
    static constexpr auto s_rotation = Real{94} * s_pi / Real{180};

    static void
    mcon(Real& epsilon, Real& infiny, Real& smalno, Real& base)
    {
        base = std::numeric_limits<Real>::radix;
        epsilon = std::numeric_limits<Real>::epsilon();
        infiny = std::numeric_limits<Real>::max();
        smalno = std::numeric_limits<Real>::min();
    }

    int
    solve(std::vector<Cmplx>& zero)
    {
        int cnt1, cnt2, i;
        bool converged;
        Cmplx z;

        Real smalno, base;
        mcon(epsilon, infin, smalno, base);
        add_rel_err = epsilon;
        mul_rel_err = 2 * s_sqrt2 * epsilon;
        Real xx = 1 / s_sqrt2;
        Real yy = -xx;
        Real cosr = std::cos(s_rotation);
        Real sinr = std::sin(s_rotation);

        zero.reserve(degree + 1);

        // Remove the zeros at the origin if any
        while (p[degree] == ZERO)
        {
            zero.push_back(ZERO);
            --degree;
        }

        // Get scales.
        for (i = 0; i <= degree; ++i)
        {
            //p[i] = op[i];
            sh[i] = std::abs(p[i]);
        }

        // Scale the polynomial
        auto bound = scale(degree, sh, epsilon, infin, smalno, base);
        if (bound != Real{1})
            for (i = 0; i <= degree; ++i)
                p[i] *= bound;

    search: 
        if (degree <= 1)
        {
            zero.push_back(-p[1] / p[0]);
            return degree;
        }

        // Calculate bound, a lower bound on the modulus of the zeros
        for (i = 0; i <= degree; ++i)
            sh[i] = std::abs(p[i]);

        bound = cauchy(degree, sh);

        // Outer loop to control 2 Major passes with different sequences of shifts
        for (cnt1 = 1; cnt1 <= 2; ++cnt1)
        {
            // First stage  calculation , no shift
            no_shift(5);

            // Inner loop to select a shift
            for (cnt2 = 1; cnt2 <= 9; ++cnt2)
            {
                // Shift is chosen with modulus bound and amplitude rotated by 94 degree from the previous shift
                auto xxx = cosr * xx - sinr * yy;
                yy = sinr * xx + cosr * yy;
                xx = xxx;
                s = bound * Cmplx(xx, yy);

                // Second stage calculation, fixed shift
                fixed_shift(10 * cnt2, z, converged);
                if (converged)
                {
                    // The second stage jumps directly to the third stage ieration.
                    // If successful the zero is stored and the polynomial deflated.
                    zero.push_back(z);
                    --degree;
                    for (i = 0; i <= degree; ++i)
                    {
                        p[i] = qp[i];
                    }
                    goto search;
                }
                // If the iteration is unsuccessful another shift is chosen
            }
            // If 9 shifts fail, the outer loop is repeated with another sequence of shifts
        }

        // The zerofinder has failed on two major passes
        // return empty handed with the number of roots found (less than the original degree)

        return degree;       
    }

    /**
     * Compute the derivative polynomial as the initial h
     * polynomial and computes num_no_shift_iters no-shift h polynomials.
     */
    void
    no_shift(int num_no_shift_iters)
    {
        auto n = degree;
        auto nm1 = n - 1;
        for (int i = 0; i < n; ++i)
        {
            Real xni = degree - i;
            h[i] = xni * p[i] / Real(n);
        }
        for (int jj = 1; jj <= num_no_shift_iters; ++jj)
        {
            if (std::abs(h[n - 1]) > epsilon * 10 * std::abs(p[n - 1]))
            {
                t = -p[degree] / h[n - 1];
                for (int i = 0; i < nm1; ++i)
                {
                    int j = degree - i - 1;
                    auto tt = h[j - 1];
                    h[j] = t * tt + p[j];
                }
                h[0] = p[0];
            }
            else
            {
                // If the constant term is essentially zero, shift H coefficients
                for (int i = 0; i < nm1; ++i)
                {
                    int j = degree - i - 1;
                    h[j] = h[j - 1];
                }
                h[0] = ZERO;
            }
        }
    }

    /**
     * Computes num_fixed_shift_iters fixed-shift h polynomials and tests for convergence.
     * Initiates a variable-shift iteration and returns with the approximate zero if successful.
     * @param[in] num_fixed_shift_iters  Limit of fixed shift steps
     * @param[out] z  - Approximate zero if converged is true
     * @param[out] converged  - Boolean indicating convergence of stage 3 iteration
     */
    void
    fixed_shift(int num_fixed_shift_iters, Cmplx& z, bool& converged)
    {
        auto n = degree;
        poly_eval(degree, s, p, qp, pv);
        bool test = true;
        bool pasd = false;

        // Calculate first T = -P(S)/H(S)
        bool h_is_tiny;
        calc_t(h_is_tiny);

        // Main loop for second stage
        for (int j = 1; j <= num_fixed_shift_iters; ++j)
        {
            auto ot = t;

            // Compute the next H Polynomial and new t
            next_h_poly(h_is_tiny);
            calc_t(h_is_tiny);
            z = s + t;

            // Test for convergence unless stage 3 has failed once or this
            // is the last H Polynomial
            if (!(h_is_tiny || !test || j == 12))
            {
                if (std::abs(t - ot) < 0.5 * std::abs(z))
                {
                    if (pasd)
                    {
                        // The weak convergence test has been passwed twice, start the third stage
                        // Iteration, after saving the current H polynomial and shift
                        for (int i = 0; i < n; ++i)
                        {
                            sh[i] = h[i];
                        }
                        auto svs = s;
                        variable_shift(10, z, converged);
                        if (converged)
                            return;

                        //The iteration failed to converge. Turn off testing and restore h,s,pv and T
                        test = 0;
                        for (int i = 0; i < n; ++i)
                        {
                            h[i] = sh[i];
                        }
                        s = svs;
                        poly_eval(degree, s, p, qp, pv);
                        calc_t(h_is_tiny);
                        continue;
                    }
                    pasd = true;
                }
                else
                    pasd = false;
            }
        }

        // Attempt an iteration with final H polynomial from second stage
        variable_shift(10, z, converged);
    }

    /**
     * Carries out the third stage iteration.
     *
     * @param[in] num_variable_shift_iters    Limit of steps in stage 3.
     * @param[inout] z  On entry contains the initial iterate, if the
     *                  iteration converges it contains the final iterate on exit.
     * @param[out] converged   True if iteration converges
     */
    void
    variable_shift(int num_variable_shift_iters, Cmplx &z, bool &converged)
    {
        bool h_is_tiny;
        Real omp, relstp;

        converged = false;
        bool b = false;
        s = z;

        // Main loop for stage three
        for (int i = 1; i <= num_variable_shift_iters; ++i)
        {
            // Evaluate P at S and test for convergence
            poly_eval(degree, s, p, qp, pv);
            auto mp = std::abs(pv);
            auto ms = std::abs(s);
            if (mp <= 20 * error_eval(degree, qp, ms, mp, add_rel_err, mul_rel_err))
            {
                // Polynomial value is smaller in value than a bound onthe error
                // in evaluationg P, terminate the ietartion
                converged = true;
                z = s;
                return;
            }
            if (i != 1)
            {
                if (!(b || mp < omp || relstp >= 0.05))
                {
                    // Iteration has stalled. Probably a cluster of zeros. Do 5 fixed 
                    // shift steps into the cluster to force one zero to dominate
                    Real tp = relstp;
                    b = 1;
                    if (relstp < epsilon)
                        tp = epsilon;
                    Real r1 = std::sqrt(tp);
                    Real r2 = s.real() * (1 + r1) - s.imag() * r1;
                    s.imag(s.real() * r1 + s.imag() * (1 + r1));
                    s.real(r2);
                    poly_eval(degree, s, p, qp, pv);
                    for (int j = 1; j <= 5; ++j)
                    {
                        calc_t(h_is_tiny);
                        next_h_poly(h_is_tiny);
                    }
                    omp = infin;
                    goto _20;
                }

                // Exit if polynomial value increase significantly
                if (mp * 0.1 > omp)
                    return;
            }

            omp = mp;

            // Calculate next iterate
      _20:  calc_t(h_is_tiny);
            next_h_poly(h_is_tiny);
            calc_t(h_is_tiny);
            if (!h_is_tiny)
            {
                relstp = std::abs(t) / std::abs(s);
                s += t;
            }
        }
    }

    /**
     * Compute  t = -p(s)/h(s).
     * @param[out]  h_is_tiny  True if h(s) is essentially zero.
     */
    void
    calc_t(bool& h_is_tiny)
    {
        auto n = degree;

        // evaluate h(s)
        Cmplx hv;
        poly_eval(n - 1, s, h, qh, hv);
        h_is_tiny = std::abs(hv) <= add_rel_err * 10 * std::abs(h[n - 1]) ? true : false;
        if (!h_is_tiny)
        {
            t = -pv / hv;
            return;
        }

        t = ZERO;
    }

    /**
     * Calculate the next shifted H polynomial.
     *
     * @param[in]  h_is_tiny  True if h(s) is essentially zero
     */
    void
    next_h_poly(bool h_is_tiny)
    {
        auto n = degree;
        if (!h_is_tiny)
        {
            for(int j = 1; j < n; ++j)
            {
                auto tt = qh[j - 1];
                h[j] = t * tt + qp[j];
            }
            h[0] = qp[0];
            return;
        }

        // If h[s] is zero replace h with qh
        for (int j = 1; j < n; ++j)
        {
            h[j] = qh[j - 1];
        }
        h[0] = ZERO;
    }

    /**
     * Evaluates a polynomial @c p at @c s by the Horner recurrence
     * placing the partial sums in q and the computed value in pv.
     */
    void
    poly_eval(int nn, const Cmplx& s, const std::vector<Cmplx>& p,
           std::vector<Cmplx>& q, Cmplx &pv)  
    {
        q[0] = p[0];
        pv = q[0];

        for (int i = 1; i <= nn; ++i)
        {
            pv = pv * s + p[i];
            q[i] = pv;
        }
    }

    /**
     * Bounds the error in evaluating the polynomial by the Horner recurrence.
     *
     * @param[in]  q  The partial sums
     * @param[in]  ms  Modulus of the point
     * @param[in]  mp  Modulus of polynomial value
     * @param[in]  are  Error bound on complex addition
     * @param[in]  mre  Error bound on complex multiplication
     */
    Real
    error_eval(int nn, const std::vector<Cmplx>& q,
          Real ms, Real mp, Real add_rel_err, Real mul_rel_err)
    {
        auto e = std::abs(q[0]) * mul_rel_err / (add_rel_err + mul_rel_err);
        for(int i = 0; i <= nn; ++i)
            e = e * ms + std::abs(q[i]);

        return e * (add_rel_err + mul_rel_err) - mp * mul_rel_err;
    }

    /**
     * Cauchy computes a lower bound on the moduli of the zeros of a polynomial.
     * @param[in]  pt  The modulus of the coefficients.
     */
    Real
    cauchy(int nn, std::vector<Cmplx>& pt)
    {
        pt[nn].real(-pt[nn].real());

        // Compute upper estimate bound
        auto n = nn;
        auto x = std::exp(std::log(-pt[nn].real()) - std::log(pt[0].real())) / n;
        if (pt[n - 1].real() != 0)
        {
            // Newton step at the origin is better, use it
            auto xm = -pt[nn].real() / pt[n - 1].real();
            if (xm < x)
                x = xm;
        }

        // Chop the interval (0,x) until f < 0
        while (true)
        {
            auto xm = x * 0.1;
            auto f = pt[0].real();
            for (int i = 1; i <= nn; ++i)
                f = f * xm + pt[i].real();
            if (f <= 0)
                break;
            x = xm;
        }
        auto dx = x;

        // Do Newton iteration until x converges to two decimal places
        while (std::fabs(dx / x) > 0.005)
        {
            pt[0].imag(pt[0].real());
            for (int i = 1; i <= nn; ++i)
                pt[i].imag(pt[i - 1].imag() * x + pt[i].real());
            auto f = pt[nn].imag();
            auto df = pt[0].imag();
            for (int i = 1; i < n; ++i)
                df = df * x + pt[i].imag();
            dx = f / df;
            x -= dx;
        }

        return x;
    }

    // Returns a scale factor to multiply the coefficients of the polynomial.
    // The scaling is done to avoid overflow and to avoid undetected underflow
    // interfering with the convergence criterion.  The factor is a power of the base.
    // pt - MODULUS OF COEFFICIENTS OF P
    // epsilon, INFIN, SMALNO, BASE - CONSTANTS DESCRIBING THE FLOATING POINT ARITHMETIC.
    Real
    scale(int nn, const std::vector<Cmplx>& pt,
          Real epsilon, Real infin, Real smalno, Real base)
    {
        int i, l;
        Real hi, lo, max, min, x, sc;
        Real fn_val;

        // Find largest and smallest moduli of coefficients
        hi = std::sqrt(infin);
        lo = smalno / epsilon;
        max = 0;
        min = infin;

        for (i = 0; i <= nn; ++i)
        {
            x = std::real(pt[i]);
            if (x > max)
                max = x;
            if (x != 0 && x < min)
                min = x;
        }

        // Scale only if there are very large or very small components
        fn_val = 1;
        if (min >= lo && max <= hi)
            return fn_val;
        x = lo / min;
        if (x <= 1)
            sc = 1 / (sqrt(max)* sqrt(min));
        else
        {
            sc = x;
            if (infin / sc > max)
                sc = 1;
        }
        l = int(std::log(sc) / std::log(base) + 0.5);
        fn_val = std::pow(base, l);
        return fn_val;
    }
};

} // namespace __gnu_cxx

#endif // SOLVER_JENKINS_TRAUB_COMPLEX_TCC
