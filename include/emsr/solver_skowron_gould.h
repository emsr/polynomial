
constexpr int MR = 8;
constexpr int MT = 10;
constexpr int MAXIT = MT * MR;
constexpr int MAXM = 30;

template <typename Tp>
void cmplx_roots_gen(std::vector<std::complex<Tp>> &roots, std::vector<std::complex<Tp>> &poly, int degree,
                     bool polish_roots_after, bool use_roots_as_starting_points);

#include "solver_skowron_gould.tcc"
