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
     */
    SolverMadsenReid(const std::vector<Cmplx>& a_in)
    : poly(a_in),
      poly_work(a_in.size())
    {}

    /**
     * Solve.
     */
    std::vector<Cmplx>
    solve()
    {
        std::vector<Cmplx> root;
        if (poly.size() <= 1)
            return root;

        int m = poly.size() - 2;
        root.resize(m + 1);
        solve(poly, m, root, poly_work);
        return root;
    }

  private:

    /// Big-endian polynomial.
    std::vector<Cmplx> poly;

    /// Big-endian working polynomial.
    std::vector<Cmplx> poly_work;

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
     * Deflate the polynomial.
     */
    void
    deflate(std::vector<Cmplx>& a, int n, Cmplx z)
    {
        for (int k = 2; k <= n; ++k)
        {
            a[k] = a[k - 1] * z + a[k];
        }
    }

    /**
     * Root search...
     */
    void solve(std::vector<Cmplx>& a1, int m, std::vector<Cmplx>& root, std::vector<Cmplx>& a);
  };

#include <ext/solver_madsen_reid.tcc>

#endif // SOLVER_MADSEN_REID_H
