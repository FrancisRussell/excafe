#ifndef SIMPLE_CFD_CAPTURE_NUMERIC_CAST_HPP
#define SIMPLE_CFD_CAPTURE_NUMERIC_CAST_HPP

#include <simple_cfd/mp/integer.hpp>
#include <simple_cfd/mp/rational.hpp>
#include <simple_cfd/mp/float.hpp>
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
struct RawConverter<mp::Float, float>
{
  static float low_level_convert(const mp::Float& s) 
  { 
    return s.toFloat();
  }
};

template<>
struct RawConverter<mp::Float, double>
{
  static double low_level_convert(const mp::Float& s) 
  { 
    return s.toDouble();
  }
};

template<>
struct RawConverter<mp::Integer, float>
{
  static float low_level_convert(const mp::Integer& s) 
  { 
    return s.toFloat();
  }
};

template<>
struct RawConverter<mp::Integer, double>
{
  static double low_level_convert(const mp::Integer& s) 
  { 
    return s.toDouble();
  }
};

template<>
struct RawConverter<mp::Rational, float>
{
  static float low_level_convert(const mp::Rational& s) 
  { 
    return s.toFloat();
  }
};

template<>
struct RawConverter<mp::Rational, double>
{
  static double low_level_convert(const mp::Rational& s) 
  { 
    return s.toDouble();
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
