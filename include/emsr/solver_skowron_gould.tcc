
#include <cmath>
#include <complex>
#include <numbers>
#include <vector>

template <typename Tp>
void cmplx_laguerre2newton(std::vector<std::complex<Tp>> &poly, int degree, std::complex<Tp> &root, int &iter, bool &success,
                           int starting_mode);
template <typename Tp> void sort_5_points_by_separation(std::vector<std::complex<Tp>> &points);
template <typename Tp>
void sort_5_points_by_separation_i(std::vector<int> &sorted_points, std::vector<std::complex<Tp>> &points);
template <typename Tp> void find_2_closest_from_5(int &i1, int &i2, Tp &d2min, std::vector<std::complex<Tp>> &points);
template <typename Tp>
void cmplx_laguerre(std::vector<std::complex<Tp>> &poly, int degree, std::complex<Tp> &root, int &iter, bool &success);
template <typename Tp>
void cmplx_newton_spec(std::vector<std::complex<Tp>> &poly, int degree, std::complex<Tp> &root, int &iter, bool &success);
template <typename Tp>
void solve_quadratic_eq(std::complex<Tp> &x0, std::complex<Tp> &x1, std::vector<std::complex<Tp>> &poly);
template <typename Tp>
void solve_cubic_eq(std::complex<Tp> &x0, std::complex<Tp> &x1, std::complex<Tp> &x2, std::vector<std::complex<Tp>> &poly);
template <typename Tp>
void divide_poly_1(std::vector<std::complex<Tp>> &polyout, std::complex<Tp> &remainder, std::complex<Tp> &p,
                   std::vector<std::complex<Tp>> &polyin, int degree);

/**
 * @param roots  Array which will hold all roots that had been found.
 *               If the flag @c use_roots_as_starting_points is set to
 *               true, then instead of point(0, 0) we use value from
 *               this array as starting point for cmplx_laguerre
 * @param poly  An array of polynomial cooefs, length = degree + 1,
 *               P(x) = poly[0] x ^ 0 + poly[1] x ^ 1 + poly[2] x ^ 2 + ...
 * @param degree  Degree of the polynomial and size of roots array
 * @param polish_roots_after  After all roots have been found by dividing
 *               original polynomial by each root found,
 *               you can opt in to polish all roots using full polynomial
 * @param  use_roots_as_starting_points  usually we start Laguerre's method from point(0, 0),
 *               but you can decide to use the values of 'roots' array as starting point
 *               for each new root that is searched for. This is useful if you have
 *               very rough idea where some of the roots can be.
 */
template <typename Tp>
void
cmplx_roots_gen(std::vector<std::complex<Tp>> &roots, std::vector<std::complex<Tp>> &poly, int degree, bool polish_roots_after,
                bool use_roots_as_starting_points)
{
    std::vector<std::complex<Tp>> poly2(MAXM);
    static int i, j, n, iter;
    bool success;
    std::complex<Tp> coef, prev;

    if (!use_roots_as_starting_points)
    {
        for (int jj = 0; jj < degree; ++jj)
        {
            roots[jj] = std::complex<Tp>(0, 0);
        }
    }

    for (j = 0; j <= degree; ++j)
        poly2[j] = poly[j];

    // Don't do Laguerre's for small degree polynomials
    if (degree <= 1)
    {
        if (degree == 1)
            roots[0] = -poly[0] / poly[1];
        return;
    }

    for (n = degree; n >= 3; --n)
    {
        cmplx_laguerre2newton(poly2, n, roots[n - 1], iter, success, 2);
        if (!success)
        {
            printf("Not a success!\n");
            roots[n - 1] = std::complex<Tp>(0, 0);
            cmplx_laguerre(poly2, n, roots[n - 1], iter, success);
        }

        // Divide by root
        // for (i = 0; i <= degree; ++i)
        //  printf("Before Division...\npoly coef%d = (%f, %f)\n", i, poly[i].real(), poly[i].imag());
        coef = poly2[n];
        for (i = n - 1; i >= 0; --i)
        {
            prev = poly2[i];
            poly2[i] = coef;
            coef = prev + roots[n - 1] * coef;
        }
        // for (i = 0; i <= degree; ++i)
        //   printf("After Division...\npoly coef%d = (%f, %f)\n", i, poly[i].real(), poly[i].imag());
    }

    // Find the to last 2 roots
    cmplx_laguerre2newton(poly2, 2, roots[1], iter, success, 2);
    if (!success)
    {
        printf("Not a success!\n");
        solve_quadratic_eq(roots[1], roots[0], poly2);
    }
    else
    {
        roots[0] = -(roots[1] + poly2[1] / poly2[2]); // Viete's Formula for the last root
    }
    if (polish_roots_after)
    {
        for (n = 0; n < degree; ++n)
        {
            cmplx_newton_spec(poly, degree, roots[n], iter, success); // Polish roots with full polynomial
            if (!success)
            {
                printf("not a success\n");
            }
        }
    }

    return;
}

/**
 * Sort array of five points
 * Most isolated point will become the first point in the array
 * The closest points will be the last two points in the array
 *
 * Algorithm works well for all dimensions. We put n = 5 as
 * a hardcoded value just for optimization purposes
 */
template <typename Tp>
void
sort_5_points_by_separation(std::vector<std::complex<Tp>> &points)
{
    const int n = 5; // works for different n as well, but is faster for n as constant (optimization)
    int sorted_points[n];

    std::complex<Tp> savepoints[n];

    sort_5_points_by_separation_i(sorted_points, points);
    for (int i = 0; i < n; ++i)
    {
        savepoints[i] = points[i];
    }
    for (int i = 1; i <= n; ++i)
    {
        points[i - 1] = savepoints[sorted_points[i - 1] - 1];
    }
    return;
}

/**
 * Return index array that sorts array of five points
 * Index of the most isolated point will appear on the firts place
 * of the output array.
 * The indices of the closest 2 points will be at the last two
 * places in the 'sorted_points' array
 *
 * Algorithm works well for all dimensions. We put n=5 as
 * a hardcoded value just for optimization purposes.
 */
template <typename Tp>
void
sort_5_points_by_separation_i(std::vector<int> &sorted_points, std::vector<std::complex<Tp>> &points)
{
    const int n = 5;

    Tp d1, d2, d;
    Tp distances2[n][n];
    int ind2, put;
    Tp neigh1st[n], neigh2nd[n];
    std::complex<Tp> p;

    for (int kj = 1; kj <= n; ++kj)
    {
        for (int ki = 1; ki <= kj - 1; ++ki)
        {
            p = points[ki - 1] - points[kj - 1];
            d = std::norm(p);
            distances2[ki - 1][kj - 1] = d;
            distances2[kj - 1][ki - 1] = d;
        }
    }

    // find neighbours
    for (int kj = 1; kj <= n; ++kj)
    {
        for (int ki = 1; ki <= kj - 1; ++ki)
        {
            d = distances2[kj - 1][ki - 1];
            if (d < neigh2nd[kj - 1])
            {
                if (d < neigh1st[kj - 1])
                {
                    neigh2nd[kj - 1] = neigh1st[kj - 1];
                    neigh1st[kj - 1] = d;
                }
                else
                {
                    neigh2nd[kj - 1] = d;
                }
            }
        }
    }

    // initialize sorted_points
    for (int ki = 1; ki <= n; ++ki)
    {
        sorted_points[ki - 1] = ki;
    }

    // sort the rest 1..n-2
    for (int kj = 2; kj <= n; ++kj)
    {
        d1 = neigh1st[kj - 1];
        d2 = neigh2nd[kj - 1];
        put = 1;
        for (int ki = kj - 1; ki >= 1; --ki)
        {
            ind2 = sorted_points[ki - 1];
            d = neigh1st[ind2 - 1];
            if (d >= d1)
            {
                if (d == d1)
                {
                    if (neigh2nd[ind2 - 1] > d2)
                    {
                        put = ki + 1;
                        break;
                    }
                }
                else
                {
                    put = ki + 1;
                    break;
                }
            }
            sorted_points[ki] = sorted_points[ki - 1];
        }
        sorted_points[put - 1] = kj;
    }
    return;
}

/**
 * Returns indices of the two cloest points out of array of 5
 * d2min is the square of minimal distance
 */
template <typename Tp>
void
find_2_closest_from_5(int &i1, int &i2, Tp &d2min, std::vector<std::complex<Tp>> &points)
{
    int n = 5;
    // Tp distances2[n][n];

    Tp d2min1, d2;
    int i, j;
    std::complex<Tp> p;
    d2min1 = Tp{ 1 };
    for (j = 1; j <= n; ++j)
    {
        for (i = 1; i <= j - 1; ++i)
        {
            p = points[i - 1] - points[j - 1];
            d2 = real(conj(p) * p);
            // distances2[i][j] = d2;
            // distances2[j][i] = d2;
            if (d2 <= d2min1)
            {
                i1 = i;
                i2 = j;
                d2min1 = d2;
            }
        }
    }
    d2min = d2min1;
}

/**
 * Subroutine finds one root of a complex polynomial using Laguerre's method.
 * In every loop it calculates simplified Adams' stopping criterion
 * for the value of the polynomial.
 *
 * Uses 'root' value as a starting point(!!!!!)
 * Remember to initialize 'root' to some initial guess or to
 * point(0, 0) if you have no prior knowledge.
 *
 * @param  poly An array of polynomial cooefs
 *
 * length = degree + 1, poly(1) is constant
 *     1              2                3
 * poly(1) x ^ 0 + poly(2) x ^ 1 + poly(3) x ^ 2 + ...
 *
 * @param  degree  Degree of the polynomial
 *
 * root - input: guess for the value of a root
 * output : a root of the polynomial
 * iter - number of iterations performed(the number of polynomial
 * evaluations and stopping criterion evaluation)
 *
 * success - is false if routine reaches maximum number of iterations
 *
 * For a summary of the method go to :
 * http://en.wikipedia.org/wiki/Laguerre's_method
 *
 */
template <typename Tp>
void
cmplx_laguerre(std::vector<std::complex<Tp>> &poly, int degree, std::complex<Tp> &root, int &iter, bool &success)
{
    static int FRAC_JUMP_EVERY = 10;
    constexpr int FRAC_JUMP_LEN = 10;
    constexpr Tp FRAC_JUMPS[FRAC_JUMP_LEN]{
        0.64109297, 0.91577881, 0.25921289, 0.50487203, 0.08177045, 0.13653241, 0.306162, 0.37794326, 0.04618805, 0.75132137
    }; // some random numbers
    constexpr Tp PI = std::numbers::pi_v<Tp>;

    Tp faq;                                                // jump length
    Tp FRAC_ERR = 20 * std::numeric_limits<Tp>::epsilon(); // 2e-15;
    std::complex<Tp> p, dp, d2p_half;                      // value of polynomial, 1st derivative, and 2nd derivative
    static int i, j, k;
    bool good_to_go;
    std::complex<Tp> denom, denom_sqrt, dx, newroot;
    Tp ek, absroot, abs2p;
    std::complex<Tp> fac_newton, fac_extra, F_half, c_one_nth;
    Tp one_nth, n_1_nth, two_n_div_n_1;
    std::complex<Tp> c_one(1, 0);
    std::complex<Tp> zero(0, 0);
    Tp stopping_crit2;

    // EXTREME FAILSAFE! not usually needed but kept here just to be on the safe side. Takes care of first coefficient being 0
    if (true)
    {
        if (degree < 0)
        {
            printf("Error: cmplx_laguerre: degree<0");
            return;
        }
        if (poly[degree] == std::complex<Tp>(0, 0))
        {
            if (degree == 0)
                return;
            cmplx_laguerre(poly, degree - 1, root, iter, success);
        }
        if (degree <= 1)
        {
            if (degree == 0)
            {
                success = false; // we just checked if poly[0] is zero and it isnt
                printf("Warning: cmplx_laguerre: degree = 0 and poly[0] does not equal zero, no roots");
                return;
            }
            else
            {
                root = -poly[0] / poly[1];
                return;
            }
        }
    } // End of EXTREME failsafe

    good_to_go = false;
    one_nth = Tp{ 1 } / degree;
    n_1_nth = (degree - Tp{ 1 }) * one_nth;
    two_n_div_n_1 = Tp{ 2 } / n_1_nth;
    c_one_nth = std::complex<Tp>(one_nth, Tp{ 0 });
    for (i = 1; i <= MAXIT; ++i)
    {
        ek = std::abs(poly[degree]); // Preparing stopping criterion
        absroot = std::abs(root);
        // Calculate the values of polynomial and its first and second derivatives
        p = poly[degree];
        dp = zero;
        d2p_half = zero;
        for (k = degree - 1; k >= 0; --k)
        {
            d2p_half = dp + d2p_half * root;
            dp = p + dp * root;
            p = poly[k] + p * root; // b_k
                                    // Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
                                    // Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
                                    // ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
                                    // Eq 8.
            ek = absroot * ek + std::abs(p);
        }
        iter += 1;

        abs2p = real(conj(p) * p);
        if (abs2p == 0)
            return;

        stopping_crit2 = pow(FRAC_ERR * ek, Tp{ 2 });
        if (abs2p < stopping_crit2)
        {
            //(simplified a little Eq. 10 of Adams 1967)
            // do additional iteration if we are less than 10x from stopping criterion
            if (abs2p < Tp{ 0.01 } * stopping_crit2)
                return; // we are at a good place!
            else
                good_to_go = true;
        }
        else
        {
            good_to_go = false;
        }

        faq = Tp{ 1 };
        denom = zero;
        if (dp != zero)
        {
            fac_newton = p / dp;
            fac_extra = d2p_half / dp;
            F_half = fac_newton * fac_extra;
            denom_sqrt = std::sqrt(c_one - two_n_div_n_1 * F_half);

            // NEXT LINE PROBABLY CAN BE COMMENTED OUT. Check if compiler outputs positive real
            if (real(denom_sqrt) >= Tp{ 0 })
                denom = c_one_nth + n_1_nth * denom_sqrt;
            else
                denom = c_one_nth - n_1_nth * denom_sqrt;
        }

        if (std::abs(denom) == Tp{ 0 })
            dx = (absroot + Tp{ 1 }) * std::exp(std::complex<Tp>(Tp{ 0 }, FRAC_JUMPS[i % FRAC_JUMP_LEN] * 2 * PI));
        else
            dx = fac_newton / denom;

        newroot = root - dx;
        if (newroot == root)
            return; // nothing changes so return
        if (good_to_go)
        {
            root = newroot;
            return;
        }
        if (i % FRAC_JUMP_EVERY == 0)
        { // decide whether to do a jump of modified length (to break cycles)
            faq = FRAC_JUMPS[(i / FRAC_JUMP_EVERY - 1) % FRAC_JUMP_LEN];
            newroot = root - faq * dx; // do jump of semi-random length
        }
        root = newroot;
    }
    success = false; // too many iterations here
    return;
}

/**
 * Subroutine finds one root of a complex polynomial
 * Newton's method. It calculates simplified Adams' stopping
 * criterion for the value of the polynomial once per 10 iterations (!),
 * after initial iteration. This is done to speed up calculations
 * when polishing roots that are known preety well, and stopping
 *  criterion does significantly change in their neighborhood.
 *
 * Uses 'root' value as a starting point (!!!!!)
 * Remember to initialize 'root' to some initial guess.
 * Do not initilize 'root' to point (0,0) if the polynomial
 * coefficients are strictly real, because it will make going
 * to imaginary roots impossible.
 *
 * @param  poly An array of polynomial cooefs
 *              length = degree+1, poly(1) is constant:
 *              P(x) = poly[0] x^0 + poly[1] x^1 + poly[2] x^2 + ...
 * @param  degree  Degree of the polynomial
 * @param  root  Input: guess for the value of a root
 *               Output: a root of the polynomial
 * @param  iter  Number of iterations performed (the number of polynomial evaluations)
 * @param  success  False if routine reaches maximum number of iterations
 *
 * For a summary of the method go to:
 * http://en.wikipedia.org/wiki/Newton's_method
 */
template <typename Tp>
void
cmplx_newton_spec(std::vector<std::complex<Tp>> &poly, int degree, std::complex<Tp> &root, int &iter, bool &success)
{
    constexpr Tp PI = std::numbers::pi_v<Tp>;
    int FRAC_JUMP_EVERY = 10;
    constexpr int FRAC_JUMP_LEN = 10;
    Tp FRAC_JUMPS[FRAC_JUMP_LEN] = { 0.64109297, 0.91577881, 0.25921289, 0.50487203, 0.08177045,
                                     0.13653241, 0.306162,   0.37794326, 0.04618805, 0.75132137 }; // some random numbers
    Tp faq;                                                                                        // jump length
    Tp FRAC_ERR = 20 * std::numeric_limits<Tp>::epsilon();                                         // 2e-15;
    std::complex<Tp> p;                                                                            // value of polynomial
    std::complex<Tp> dp;                                                                           // value of 1st derivative
    int i, k;
    bool good_to_go;
    std::complex<Tp> dx, newroot;
    Tp ek, absroot, abs2p;
    std::complex<Tp> zero(0, 0);
    Tp stopping_crit2;

    iter = 0;
    success = true;

    // the next if block is an EXTREME failsafe, not usually needed, and thus turned off in this version
    if (true)
    { // change false to true if you would like to use caustion about haveing first coefficient == 0
        if (degree < 0)
        {
            printf("Error: cmplx_newton_spec: degree<0");
            return;
        }
        if (poly[degree] == zero)
        {
            if (degree == 0)
                return;
            cmplx_newton_spec(poly, degree, root, iter, success);
            return;
        }
        if (degree <= 1)
        {
            if (degree == 0)
            {
                success = false;
                printf("Warning: cmplx_newton_spec: degree=0 and poly[0]!=0, no roots");
                return;
            }
            else
            {
                root = -poly[0] / poly[1];
                return;
            }
        }
    }
    // end EXTREME Failsafe
    good_to_go = false;

    stopping_crit2 = Tp{ 0 }; // value not important, will be initialized anyway on the first loop
    for (i = 1; i <= MAXIT; ++i)
    {
        faq = Tp{ 1 };
        // prepare stoping criterion
        // calculate value of polynomial and its first two derivatives
        p = poly[degree];
        dp = zero;
        if (i % 10 == 1)
        { // calculate stopping criterion every tenth iteration
            ek = std::abs(poly[degree]);
            absroot = std::abs(root);
            for (k = degree - 1; k >= 0; --k)
            {
                dp = p + dp * root;
                p = poly[k] + p * root; // b_k
                                        // Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
                                        // Communications of ACM, Volume 10 Issue 10, Oct. 1967, p. 655
                                        // ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
                                        // Eq. 8
                ek = absroot * ek + std::abs(p);
            }
            stopping_crit2 = pow(FRAC_ERR * ek, 2);
        }
        else
        { // calculate just the value and derivative
            for (k = degree - 1; k >= 0; --k)
            { // Horner Scheme, see for eg. Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
                dp = p + dp * root;
                p = poly[k] + p * root;
            }
        }

        iter = iter + 1;

        abs2p = real(conj(p) * p);
        if (abs2p == Tp{ 0 })
            return;

        if (abs2p < stopping_crit2)
        { // simplified a little Eq. 10 of Adams 1967
            if (dp == zero)
                return; // if we have problem with zero, but we are close to the root, just accept
                        // do additional iteration if we are less than 10x from stopping criterion
            if (abs2p < 0.01 * stopping_crit2)
                return; // return immediatley because we are at very good place
            else
                good_to_go = true; // do one iteration more
        }

        else
        {
            good_to_go = false; // reset if we are outside the zone of the root
        }

        if (dp == zero)
            // problem with zero
            dx = (std::abs(root) + Tp{ 1 })
                 * std::exp(std::complex<Tp>(Tp{ 0 }, FRAC_JUMPS[i % FRAC_JUMP_LEN] * Tp{ 2 } * PI));
        else
            dx = p / dp; // Newton method, see http://en.wikipedia.org/wiki/Newton's_method

        newroot = root - dx;
        if (newroot == root)
            return; // nothing changes -> return
        if (good_to_go)
        { // this was jump already after stopping criterion was met
            root = newroot;
            return;
        }
        if (i % FRAC_JUMP_EVERY == 0)
        { // decide whether to do a jump of modified length (to break cycles)
            faq = FRAC_JUMPS[(i / FRAC_JUMP_EVERY - 1) % FRAC_JUMP_LEN];
            newroot = root - faq * dx;
        }
        root = newroot;
    }
    success = false;
    return;
    // too many iterations here
}

/**
 * Subroutine finds one root of a complex polynomial using
 * Laguerre's method, Second-order General method and Newton's
 * method - depending on the value of function F, which is a
 * combination of second derivative, first derivative and
 * value of polynomial [F=-(p"*p)/(p'p')].
 *
 * Subroutine has 3 modes of operation. It starts with mode=2
 * which is the Laguerre's method, and continues until F
 * becames F<0.50, at which point, it switches to mode=1,
 * i.e., SG method (see paper). While in the first two
 * modes, routine calculates stopping criterion once per every
 * iteration. Switch to the last mode, Newton's method, (mode=0)
 * happens when becomes F<0.05. In this mode, routine calculates
 * stopping criterion only once, at the beginning, under an
 * assumption that we are already very close to the root.
 * If there are more than 10 iterations in Newton's mode,
 * it means that in fact we were far from the root, and
 * routine goes back to Laguerre's method (mode=2).
 *
 * Uses 'root' value as a starting point (!!!!!)
 * Remember to initialize 'root' to some initial guess or to
 * point (0,0) if you have no prior knowledge.
 *
 * @param  poly - is an array of polynomial cooefs
 *                P(x) = poly[0] x^0 + poly[1] x^1 + poly[2] x^2
 * @param  degree  Degree of the polynomial
 * @param  root  Input: guess for the value of a root
 *               Output: a root of the polynomial
 * @param  iter  Number of iterations performed (the number of polynomial
 *               evaluations and stopping criterion evaluation)
 * @param  success  False if routine reaches maximum number of iterations
 * @param  starting_mode  This should be by default = 2. However if you
 *                   choose to start with SG method put 1 instead.
 *                   Zero will cause the routine to
 *                   start with Newton for first 10 iterations, and
 *                   then go back to mode 2.
 *
 * For a summary of the method see the paper: Skowron & Gould (2012)
 */
template <typename Tp>
void
cmplx_laguerre2newton(std::vector<std::complex<Tp>> &poly, int degree, std::complex<Tp> &root, int &iter, bool &success,
                      int starting_mode)
{
    constexpr Tp PI = std::numbers::pi_v<Tp>;
    int FRAC_JUMP_EVERY = 10;
    const int FRAC_JUMP_LEN = 10;
    Tp FRAC_JUMPS[FRAC_JUMP_LEN] = { 0.64109297, 0.91577881, 0.25921289, 0.50487203, 0.08177045,
                                     0.13653241, 0.306162,   0.37794326, 0.04618805, 0.75132137 }; // some random numbers

    Tp faq;                                                // jump length
    Tp FRAC_ERR = 20 * std::numeric_limits<Tp>::epsilon(); // 2e-15;

    std::complex<Tp> p;        // value of polynomial
    std::complex<Tp> dp;       // value of 1st derivative
    std::complex<Tp> d2p_half; // value of 2nd derivative
    int i, j, k;
    bool good_to_go;
    // std::complex<Tp> G, H, G2;
    std::complex<Tp> denom, denom_sqrt, dx, newroot;
    Tp ek, absroot, abs2p, abs2_F_half;
    std::complex<Tp> fac_netwon, fac_extra, F_half, c_one_nth;
    Tp one_nth, n_1_nth, two_n_div_n_1;
    int mode;
    std::complex<Tp> c_one(1, 0);
    std::complex<Tp> zero(0, 0);
    Tp stopping_crit2;

    iter = 0;
    success = true;
    stopping_crit2 = 0; // value not important, will be initialized anyway on the first loop

    // next if block is an EXTREME failsafe, not usually needed, and thus turned off in this version.
    if (false)
    { // change false to true if you would like to use caution about having first coefficent == 0
        if (degree < 0)
        {
            printf("Error: cmplx_laguerre2newton: degree < 0");
            return;
        }
        if (poly[degree] == zero)
        {
            if (degree == 0)
                return;
            cmplx_laguerre2newton(poly, degree, root, iter, success, starting_mode);
            return;
        }
        if (degree <= 1)
        {
            if (degree == 0)
            { //// we know from previous check that poly[0] not equal zero
                success = false;
                printf("Warning: cmplx_laguerre2newton: degree = 0 and poly[0] = 0, no roots");
                return;
            }
            else
            {
                root = -poly[0] / poly[1];
                return;
            }
        }
    }
    // end EXTREME failsafe

    j = 1;
    good_to_go = false;

    mode = starting_mode; // mode = 2 full laguerre, mode = 1 SG, mode = 0 newton

    while (true)
    {
        // infinite loop, just to be able to come back from newton, if more than 10 iteration there

        // == mode 2 ==

        if (mode >= 2)
        {                                     // Laguerre's method
            one_nth = Tp{ 1 } / degree;       ///
            n_1_nth = (degree - 1) * one_nth; ////
            two_n_div_n_1 = Tp{ 2 } / n_1_nth;
            c_one_nth = std::complex<Tp>(one_nth, Tp{ 0 });

            for (i = 1; i <= MAXIT; ++i)
            {
                faq = Tp{ 1 };

                // prepare stoping criterion
                ek = std::abs(poly[degree]);
                absroot = std::abs(root);
                // calculate value of polynomial and its first two derivative
                p = poly[degree];
                dp = zero;
                d2p_half = zero;
                for (k = degree; k >= 1; --k)
                { // Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
                    d2p_half = dp + d2p_half * root;
                    dp = p + dp * root;
                    p = poly[k - 1] + p * root; // b_k
                                                // Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
                                                // Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
                                                // ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
                                                // Eq 8.
                    ek = absroot * ek + std::abs(p);
                }
                abs2p = std::norm(p);
                iter = iter + 1;
                if (abs2p == 0)
                    return;

                stopping_crit2 = std::pow(FRAC_ERR * ek, Tp{ 2 });
                if (abs2p < stopping_crit2)
                { //(simplified a little Eq. 10 of Adams 1967)
                  // do additional iteration if we are less than 10x from stopping criterion
                    if (abs2p < Tp{ 0.01 } * stopping_crit2)
                        return; // ten times better than stopping criterion
                                // return immediately, because we are at very good place
                    else
                        good_to_go = true; // do one iteration more
                }
                else
                {
                    good_to_go = false; // reset if we are outside the zone of the root
                }

                denom = zero;
                if (dp != zero)
                {
                    fac_netwon = p / dp;
                    fac_extra = d2p_half / dp;
                    F_half = fac_netwon * fac_extra;

                    abs2_F_half = std::norm(F_half);
                    if (abs2_F_half <= Tp{ 0.0625 })
                    { // F<0.50, F/2<0.25
                        // go to SG method
                        if (abs2_F_half <= Tp{ 0.000625 })
                        {             // F<0.05, F/2<0.02
                            mode = 0; // go to Newton's
                        }
                        else
                        {
                            mode = 1; // go to SG
                        }
                    }

                    denom_sqrt = std::sqrt(c_one - two_n_div_n_1 * F_half);

                    // NEXT LINE PROBABLY CAN BE COMMENTED OUT
                    if (std::real(denom_sqrt) > Tp{ 0 })
                    {
                        // real part of a square root is positive for probably all compilers. You can ù
                        // test this on your compiler and if so, you can omit this check
                        denom = c_one_nth + n_1_nth * denom_sqrt;
                    }
                    else
                    {
                        denom = c_one_nth - n_1_nth * denom_sqrt;
                    }
                }
                if (denom == zero)
                { // test if demoninators are > 0.0 not to divide by zero
                    dx = (std::abs(root) + Tp{ 1 })
                         + std::exp(
                             std::complex<Tp>(Tp{ 0 }, FRAC_JUMPS[i % FRAC_JUMP_LEN] * Tp{ 2 } * PI)); // make some random jump
                }
                else
                {
                    dx = fac_netwon / denom;
                }
                newroot = root - dx;
                if (newroot == root)
                    return; // nothing changes -> return
                if (good_to_go)
                { // this was jump already after stopping criterion was met
                    root = newroot;
                    return;
                }
                if (mode != 2)
                {
                    root = newroot;
                    j = i + 1; // remember iteration index
                    break;     // go to Newton's or SG
                }
                if (i % FRAC_JUMP_EVERY == 0)
                { // decide whether to do a jump of modified length (to break cycles)
                    faq = FRAC_JUMPS[((i / FRAC_JUMP_EVERY - 1) % FRAC_JUMP_LEN)];
                    newroot = root - faq * dx; // do jump of some semi-random length (0 < faq < 1)
                }
                root = newroot;
            } // do mode 2

            if (i >= MAXIT)
            {
                success = false;
                return;
            }
        }

        // == mode 1 ==

        if (mode == 1)
        { // SECOND-ORDER GENERAL METHOD (SG)
            for (i = j; i <= MAXIT; ++i)
            {
                faq = Tp{ 1 };
                // calculate value of polynomial and its first two derivatives
                p = poly[degree];
                dp = zero;
                d2p_half = zero;
                if ((i - j) % 10 == 0)
                {
                    // prepare stopping criterion
                    ek = std::abs(poly[degree]);
                    absroot = std::abs(root);
                    for (k = degree; k >= 1; --k)
                    { // Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
                        d2p_half = dp + d2p_half * root;
                        dp = p + dp * root;
                        p = poly[k - 1]
                            + p * root; // b_k
                                        // Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
                                        // Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
                                        // ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
                                        // Eq 8.
                        ek = absroot * ek + std::abs(p);
                    }
                    stopping_crit2 = pow(FRAC_ERR * ek, 2);
                }
                else
                {
                    for (k = degree; k >= 1; --k)
                    { // Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
                        d2p_half = dp + d2p_half * root;
                        dp = p + dp * root;
                        p = poly[k - 1] + p * root; // b_k
                    }
                }
                abs2p = std::norm(p);
                ++iter;
                if (abs2p == Tp{ 0 })
                    return;

                if (abs2p < stopping_crit2)
                { //(simplified a little Eq. 10 of Adams 1967)
                    if (dp == zero)
                        return;
                    // do additional iteration if we are less than 10x from stopping criterion
                    if (abs2p < Tp{ 0.01 } * stopping_crit2)
                        return; // ten times better than stopping criterion
                    else
                        good_to_go = true; // do one iteration more
                }
                else
                {
                    good_to_go = false; // reset if we are outside the zone of the root
                }
                if (dp == zero)
                { // test if denominators are > 0.0 not to divide by zero
                    dx = (std::abs(root) + Tp{ 1 })
                         * std::exp(
                             std::complex<Tp>(Tp{ 0 }, FRAC_JUMPS[i % FRAC_JUMP_LEN] * Tp{ 2 } * PI)); // make some random jump
                }
                else
                {
                    fac_netwon = p / dp;
                    fac_extra = d2p_half / dp;
                    F_half = fac_netwon * fac_extra;

                    abs2_F_half = real(conj(F_half) * F_half);
                    if (abs2_F_half <= Tp{ 0.000625 })
                    {             // F<0.05, F/2<0.025
                        mode = 0; // set Newton's, go there after jump
                    }
                    dx = fac_netwon * (c_one + F_half); // SG
                }
                newroot = root - dx;
                if (newroot == root)
                    return; // nothing changes -> return
                if (good_to_go)
                {
                    root = newroot; // this was jump already after stopping criterion was met
                    return;
                }
                if (mode != 1)
                {
                    root = newroot;
                    j = i + 1; // remember iteration number
                    break;     // go to Newton's
                }
                if ((i % FRAC_JUMP_EVERY) == 0)
                { // decide whether to do a jump of modified length (to break cycles)
                    faq = FRAC_JUMPS[(i / FRAC_JUMP_EVERY - 1) % FRAC_JUMP_LEN];
                    newroot = root - faq * dx; // do jump of some semi random lenth (0 < faq < 1)
                }
                root = newroot;
            }
            if (i >= MAXIT)
            {
                success = false;
                return;
            }
        }

        // == mode 0 ==

        if (mode == 0)
        { // Newton's Method

            for (i = j; i <= j + 10; ++i)
            { // Do only 10 iterations the most then go back to Laguerre
                faq = Tp{ 1 };

                // calc polynomial and first two derivatives
                p = poly[degree];
                dp = zero;
                if (i == j)
                { // Calculating stopping criterion only at the beginning
                    ek = std::abs(poly[degree]);
                    absroot = std::abs(root);
                    for (k = degree; k >= 1; --k)
                    {
                        dp = p + dp * root;
                        p = poly[k - 1] + p * root;
                        ek = absroot * ek + std::abs(p);
                    }
                    stopping_crit2 = pow(FRAC_ERR * ek, Tp{ 2 });
                }
                else
                {
                    for (k = degree; k >= 1; --k)
                    {
                        dp = p + dp * root;
                        p = poly[k - 1] + p * root;
                    }
                }
                abs2p = std::norm(p);
                iter = iter + 1;
                if (abs2p == Tp{ 0 })
                    return;

                if (abs2p < stopping_crit2)
                {
                    if (dp == zero)
                        return;
                    // do additional iteration if we are less than 10x from stopping criterion
                    if (abs2p < Tp{ 0.01 } * stopping_crit2)
                        return; // return immediately since we are at a good place
                    else
                        good_to_go = true; // do one more iteration
                }
                else
                {
                    good_to_go = false;
                }

                if (dp == zero)
                {
                    dx = (std::abs(root) + Tp{ 1 })
                         * std::exp(
                             std::complex<Tp>(Tp{ 0 }, Tp{ 2 } * PI * FRAC_JUMPS[i % FRAC_JUMP_LEN])); // make a random jump
                }
                else
                {
                    dx = p / dp;
                }

                newroot = root - dx;
                if (newroot == root)
                    return;
                if (good_to_go)
                {
                    root = newroot;
                    return;
                }
                root = newroot;
            }
            if (iter >= MAXIT)
            {
                // too many iterations
                success = false;
                return;
            }
            mode = 2; // go back to Laguerre's. Happens when could not converge with 10 steps of Newton
        }

    } /// end of infinite loop
}

template <typename Tp>
void
solve_quadratic_eq(std::complex<Tp> &x0, std::complex<Tp> &x1, std::vector<std::complex<Tp>> &poly)
{
    std::complex<Tp> a, b, c, b2, delta;
    a = poly[2];
    b = poly[1];
    c = poly[0];
    b2 = b * b;
    delta = std::sqrt(b2 - Tp{ 4 } * a * c);
    if (std::real(std::conj(b) * delta) >= 0)
    {
        x0 = -Tp{ 0.5 } * (b + delta);
    }
    else
    {
        x0 = -Tp{ 0.5 } * (b - delta);
    }
    if (x0 == std::complex<Tp>(Tp{ 0 }, Tp{ 0 }))
    {
        x1 = std::complex<Tp>(Tp{ 0 }, Tp{ 0 });
    }
    else
    { // Viete's formula
        x1 = c / x0;
        x0 = x0 / a;
    }
    return;
    // The low-degree solvers are real coeffs only!!!! That's why we can't have nice thinga!
/*
    auto roots = emsr::quadratic<Tp>(poly);
    x0 = emsr::to_complex(roots[0]);
    x1 = emsr::to_complex(roots[1]);
*/
    return;
}

/**
 * Cubic equation solver for complex polynomial (degree=3)
 * http://en.wikipedia.org/wiki/Cubic_function   Lagrange's method
 * @param  poly  An array of polynomial cooefs, length = degree+1, poly[0] is constant
 *               P(x) = poly[0] x^0 + poly[1] x^1 + poly[2] x^2 + poly[3] x^3
 */
template <typename Tp>
void
solve_cubic_eq(std::complex<Tp> &x0, std::complex<Tp> &x1, std::complex<Tp> &x2, std::vector<std::complex<Tp>> &poly)
{
/*
    constexpr Tp SQRT3_2 = std::numbers::sqrt3_v<Tp> / Tp{ 2 };
    std::complex<Tp> zeta = std::complex<Tp>(-Tp{ 0.5 }, SQRT3_2);
    std::complex<Tp> zeta2 = std::complex<Tp>(-Tp{ 0.5 }, -SQRT3_2);
    constexpr Tp third = Tp{ 1 } / Tp{ 3 };
    std::complex<Tp> s0, s1, s2;
    std::complex<Tp> E1; // x0+x1+x2
    std::complex<Tp> E2; // x0*x1+x1*x2+x2*x0
    std::complex<Tp> E3; // x0*x1*x2
    std::complex<Tp> A, B, a_1, E12, delta, A2;

    std::complex<Tp> val, x;
    a_1 = 1 / poly[3];
    E1 = -poly[2] * a_1;
    E2 = poly[1] * a_1;
    E3 = -poly[0] * a_1;

    s0 = E1;
    E12 = E1 * E1;
    A = Tp{ 2 } * E1 * E12 - Tp{ 9 } * E1 * E2 + Tp{ 27 } * E3;
    B = E12 - Tp{ 3 } * E2;
    // quadratic equation z^2 - A * z + B^3 where roots are equal to s1^3 and s2^3
    A2 = A * A;
    delta = std::sqrt(A2 - Tp{ 4 } * (B * B * B));
    if (std::real(conj(A) * delta) >= Tp{ 0 })
    { // scalar product to decide the sign yielding bigger magnitude
        s1 = std::cbrt(Tp{ 0.5 } * (A + delta));
    }
    else
    {
        s1 = std::cbrt(Tp{ 0.5 } * (A - delta));
    }
    if (std::real(s1) == Tp{ 0 } && std::imag(s1) == Tp{ 0 })
    {
        s2 = std::complex<Tp>(0, 0);
    }
    else
    {
        s2 = B / s1;
    }

    x0 = third * (s0 + s1 + s2);
    x1 = third * (s0 + s1 * zeta2 + s2 * zeta);
    x2 = third * (s0 + s1 * zeta + s2 * zeta2);
*/
    auto roots = emsr::cubic(poly);
    x0 = roots[0];
    x1 = roots[1];
    x2 = roots[2];
    return;
}

/**
 * Subroutine will divide polynomial 'polyin' by(x - p)
 * results will be returned in polynomial 'polyout' of degree - 1
 * The remainder of the division will be returned in 'remainder'
 *
 * You can provide same array as 'polyin' and 'polyout' - this
 * routine will work fine, though it will not set to zero the
 * unused, highest coefficient in the output array.You just have
 * remember the proper degree of a polynomial.
 * @param  poly  An array of polynomial cooefs, length = degree + 1, poly(1) is constant
 *               P(x) = poly(1) x^0 + poly(2) x^1 + poly(3) x^2 + ...
 */
template <typename Tp>
void
divide_poly_1(std::vector<std::complex<Tp>> &polyout, std::complex<Tp> &remainder, std::complex<Tp> &p,
              std::vector<std::complex<Tp>> &polyin, int degree)
{
    std::complex<Tp> coef;
    coef = polyin[degree];
    for (int k = 1; k <= degree; ++k)
    {
        polyout[k - 1] = polyin[k - 1];
    }
    for (int i = degree; i >= 1; --i)
    {
        auto prev = polyout[i - 1];
        polyout[i - 1] = coef;
        coef = prev + p * coef;
    }
    remainder = coef;
    return;
}
