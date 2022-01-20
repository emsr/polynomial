/**
 * @file solver_madsen_reid.h Class declaration for the Madsen-Reid solver.
 */

/**
 * @def  SOLVER_MADSEN_REID_H
 *
 * @brief  A guard for the SolverMadsenReid class header.
 */
#ifndef SOLVER_MADSEN_REID_H
#define SOLVER_MADSEN_REID_H 1

#include <limits>
#include <complex>
#include <vector>
#include <iostream>

/**
 * Return the L1 sum of absolute values or Manhattan metric of a complex number.
 */
template<typename Real>
  inline Real
  norm_l1(std::complex<Real> z)
  {
    return std::abs(std::real(z)) + std::abs(std::imag(z));
  }

/**
 * Return the L2 modulus of the complex number (this is std::norm).
 */
template<typename Real>
  inline Real
  norm_l2(std::complex<Real> z)
  {
    return std::norm(z);
  }

template<typename Real>
  class SolverMadsenReid
  {
  public:

    using Cmplx = std::complex<Real>;

    /**
     * Constructor.
     *
     * @param a_in Coefficient of the polynomial in "big-endian" - largest degree coefficient first - order.
     */
    SolverMadsenReid(const std::vector<Cmplx>& a_in)
    : poly(a_in),
      poly_work(a_in.size())
    {}

    /**
     * Solve the polynomial.
     */
    std::vector<Cmplx>
    solve()
    {
        std::vector<Cmplx> root;
        if (poly.size() <= 1)
            return root;

        int degree = poly.size() - 1;
        root.resize(degree);
        solve(poly, degree, root, poly_work);
        return root;
    }

  private:

    static constexpr Real DIGITS = std::numeric_limits<Real>::max_digits10;
    static constexpr Real BIG = std::numeric_limits<Real>::max(); // Overflow limit
    static constexpr Real SMALL = std::numeric_limits<Real>::min(); // Underflow limit.
    static constexpr Real BASE = std::numeric_limits<Real>::radix;
    static constexpr Real EPS = std::numeric_limits<Real>::epsilon();

    /// Big-endian polynomial.
    std::vector<Cmplx> poly;

    /// Big-endian working polynomial.
    std::vector<Cmplx> poly_work;

    /**
     * Evaluate polynomial at z, set fz, return squared modulus.
     *
     * @param  z     Argument of the polynomial.
     * @param  fz    Polynomial value at the given argument.
     * @param  size  Size (degree + 1) of the polynomial polynomial.
     * @param  a     Polynomial coefficients.
     * @return  Squared modulus of the function value.
     */
    Real
    eval(Cmplx z, Cmplx& fz, int size, const std::vector<Cmplx>& a)
    {
        auto deg = size - 1;
        auto p = a[0];
        for (int i = 0; i < deg; ++i)
        {
            p = p * z + a[i + 1];
        }

        fz = p;

        return norm_l2(p);
    }

    /**
     * Store the root.
     *
     * @param  a1  Working polynomial.
     * @param  root  Roots of the polynomial.
     * @param  z  New root of the polynomial.
     */
    void
    push_root(std::vector<Cmplx>& a1, std::vector<Cmplx>& root, int& n, Cmplx z)
    {
        a1[n - 1] = root[n - 1];
        root[n - 1] = z;
        --n;
    }

    /**
     * Deflate the polynomial.
     *
     * @param  a  Polynomial
     * @param  n  New degree.
     * @param  z  New root of the polynomial.
     */
    void
    deflate(std::vector<Cmplx>& a, int n, Cmplx z)
    {
        for (int k = 1; k < n; ++k)
        {
            a[k] = a[k - 1] * z + a[k];
        }
    }

    /**
     * Root search...
     *
     * @param  a     Input polynomial
     * @param  n     Degree.
     * @param  root  Roots of the polynomial.
     * @param  a     Working polynomial.
     */
    void solve(std::vector<Cmplx>& a1, int m, std::vector<Cmplx>& root, std::vector<Cmplx>& a);
  };

#include <emsr/solver_madsen_reid.tcc>

#endif // SOLVER_MADSEN_REID_H
