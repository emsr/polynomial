
#ifndef _EXT_HORNER_H
#define _EXT_HORNER_H 1

#pragma GCC system_header

#include <type_traits>

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

template<typename _ArgT, typename _Coef0>
  constexpr std::conditional_t<std::is_integral_v<_ArgT>, double, _ArgT>
  horner(_ArgT __x, _Coef0 __c0)
  {
    using __arg_t = std::conditional_t<std::is_integral_v<_ArgT>, double, _ArgT>;
    return __arg_t{__c0};
  }

template<typename _ArgT, typename _Coef0, typename... _Coef>
  constexpr std::conditional_t<std::is_integral_v<_ArgT>, double, _ArgT>
  horner(_ArgT __x, _Coef0 __c0, _Coef... __c)
  {
    using __arg_t = std::conditional_t<std::is_integral_v<_ArgT>, double, _ArgT>;
    return __arg_t{__c0} + __x * horner(__x, __c...);
  }


template<typename _ArgT, typename _Coef0>
  constexpr std::conditional_t<std::is_integral_v<_ArgT>, double, _ArgT>
  horner_big_end(_ArgT, _Coef0 __c0)
  {
    using __arg_t = std::conditional_t<std::is_integral_v<_ArgT>, double, _ArgT>;
    return __arg_t{__c0};
  }

template<typename _ArgT, typename _Coef1, typename _Coef0>
  constexpr std::conditional_t<std::is_integral_v<_ArgT>, double, _ArgT>
  horner_big_end(_ArgT __x, _Coef1 __c1, _Coef0 __c0)
  {
    using __arg_t = std::conditional_t<std::is_integral_v<_ArgT>, double, _ArgT>;
    return horner_big_end(__x, __x * __arg_t{__c1} + __arg_t{__c0});
  }

template<typename _ArgT, typename _CoefN, typename _CoefNm1, typename... _Coef>
  constexpr std::conditional_t<std::is_integral_v<_ArgT>, double, _ArgT>
  horner_big_end(_ArgT __x, _CoefN __cn, _CoefNm1 __cnm1, _Coef... __c)
  {
    using __arg_t = std::conditional_t<std::is_integral_v<_ArgT>, double, _ArgT>;
    return horner_big_end(__x, __x * __arg_t{__cn} + __arg_t{__cnm1}, __c...);
  }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#endif // _EXT_HORNER_H
