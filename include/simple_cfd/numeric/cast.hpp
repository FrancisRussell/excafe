#ifndef SIMPLE_CFD_CAPTURE_NUMERIC_CAST_HPP
#define SIMPLE_CFD_CAPTURE_NUMERIC_CAST_HPP

#include <cln/cln.h>
#include <boost/numeric/conversion/cast.hpp>

namespace cfd
{

namespace detail
{

template<typename source_type, typename result_type>
struct RawConverter
{
  static result_type low_level_convert(source_type s) { return static_cast<result_type>(s); }
};

template<>
struct RawConverter<cln::cl_F, float>
{
  static float low_level_convert(const cln::cl_F& s) 
  { 
    return cln::float_approx(s);
  }
};

template<>
struct RawConverter<cln::cl_F, double>
{
  static double low_level_convert(const cln::cl_F& s) 
  { 
    return cln::double_approx(s);
  }
};

template<>
struct RawConverter<cln::cl_R, double>
{
  static double low_level_convert(const cln::cl_R& s) 
  { 
    return cln::double_approx(s);
  }
};

template<>
struct RawConverter<cln::cl_R, float>
{
  static double low_level_convert(const cln::cl_R& s) 
  { 
    return cln::float_approx(s);
  }
};


}


template<typename Target, typename Source> inline
typename boost::numeric::converter<Target,Source>::result_type
numeric_cast ( Source arg )
{
  using namespace boost::numeric;
  typedef conversion_traits<Target, Source> Traits;

  return boost::numeric::converter<
    Target,
    Source,
    Traits,
    def_overflow_handler,
    Trunc< typename Traits::source_type >,
    detail::RawConverter<typename Traits::source_type, typename Traits::result_type>,
    UseInternalRangeChecker
    >::convert(arg);
}


}

#endif
