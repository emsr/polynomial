#ifndef SOLVER_JENKINS_TRAUB_H
#define SOLVER_JENKINS_TRAUB_H 1

#include <ext/math_const.h> // For __const_pi, etc.
#include "solution.h" // For solution_t

namespace __gnu_cxx
{

/**
 * A solver for real-coefficient polynomials due to Jenkins and Traub.
 */
template<typename _Real>
  class _JenkinsTraubSolver
  {
  public:

    _JenkinsTraubSolver(const std::vector<_Real>& op);
    _JenkinsTraubSolver(std::vector<_Real>&& op);

    std::vector<solution_t<_Real>> solve();

  private:

    enum NormalizationType
    {
      none,
      divide_by_c,
      divide_by_d,
      near_h_root
    };

    void quadratic(_Real a, _Real b, _Real c,
		   solution_t<_Real> &z_small, solution_t<_Real> &z_large);
    int fxshfr(int l2);
    int iter_quadratic(_Real uu, _Real vv);
    int iter_real(_Real sss, int& iflag);
    NormalizationType init_next_h_poly();
    void next_h_poly(NormalizationType type);
    std::pair<_Real, _Real> quadratic_coefficients(NormalizationType type);
    void remquo_quadratic(int n, _Real u, _Real v,
			  std::vector<_Real>& poly, std::vector<_Real>& quot,
			  _Real& a, _Real& b);

    static constexpr auto _S_eps = std::numeric_limits<_Real>::epsilon();
    static constexpr auto _S_base = _Real{std::numeric_limits<_Real>::radix};
    static constexpr auto _S_tiny = _S_eps * _S_eps * _S_eps; // 1.0e-50; //std::numeric_limits<_Real>::min();
    static constexpr auto _S_huge = _Real{1} / _S_tiny; // 1.0e+50; //std::numeric_limits<_Real>::max();
    static constexpr auto _S_low = _S_tiny / _S_eps;
    static constexpr auto _S_pi = __gnu_cxx::__const_pi<_Real>();
    static constexpr auto _S_rotation = _Real{94} * _S_pi / _Real{180};

    int max_iter_quadratic = 20;
    _Real min_log_deriv = _Real{0.005L};
    int max_iter_real = 10;
    // Epsilon parameter.
    _Real __are = _S_eps;
    // Epsilon parameter.
    _Real __mre = _S_eps;

    std::vector<_Real> _P;
    std::vector<_Real> _P_quot;
    std::vector<_Real> _H, _H_quot, _H_save;
    _Real __sr, __si;
    _Real __u, __v;
    _Real __a;
    _Real __b;
    _Real __c;
    _Real __d;
    _Real __e;
    _Real __f;
    _Real __g;
    _Real __h;
    _Real __a1;
    _Real __a2;
    _Real __a3;
    _Real __a7;
    solution_t<_Real> __z_small;
    solution_t<_Real> __z_large;
    int __order;
    bool __zerok;
    int __num_iters = 0;
  };

} // namespace __gnu_cxx

#include "solver_jenkins_traub.tcc"

#endif // SOLVER_JENKINS_TRAUB_H
