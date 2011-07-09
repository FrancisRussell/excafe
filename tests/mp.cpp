#include <cassert>
#include <string>
#include <set>
#include <simple_cfd/mp/integer.hpp>
#include <simple_cfd/mp/rational.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

#define BOOST_TEST_MODULE "BigNum"
#include <boost/test/unit_test.hpp>

using namespace cfd::mp;

BOOST_AUTO_TEST_CASE(BasicSanity)
{
  std::set<int> integers;
  integers.insert(0);
  integers.insert(1);
  integers.insert(7);
  integers.insert(76);
  integers.insert(5850);
  integers.insert(53288);
  
  BOOST_FOREACH(const int intValue, integers)
  {
    const Integer integerValue(intValue);

    BOOST_CHECK_EQUAL(integerValue, integerValue);
    BOOST_CHECK_EQUAL(-integerValue, -integerValue);

    BOOST_CHECK_EQUAL(integerValue, intValue);
    BOOST_CHECK_EQUAL(-integerValue, -intValue);
    
    BOOST_CHECK_EQUAL(-integerValue*-integerValue, integerValue*integerValue);
    BOOST_CHECK_EQUAL(integerValue*1, integerValue);
    BOOST_CHECK_EQUAL(integerValue*-1, -integerValue);
    BOOST_CHECK_EQUAL(-(-integerValue), integerValue);

    if (intValue != 0)
      BOOST_CHECK_NE(-intValue, integerValue);

    const Integer less = integerValue-1;
    BOOST_CHECK_LT(less, intValue);
    BOOST_CHECK_LT(less, integerValue);

    BOOST_CHECK_GT(-less, -intValue);
    BOOST_CHECK_GT(-less, -integerValue);
  }
}

BOOST_AUTO_TEST_CASE(Powers)
{
  const Integer large = pow(Integer(898311), 5);
  Integer raised(1);

  for(int i=0; i<10; ++i)
  {
    BOOST_CHECK_EQUAL(raised, pow(large, i));
    raised *= large;
  }
}

BOOST_AUTO_TEST_CASE(IntegerSquareRoot)
{
  BOOST_CHECK_EQUAL(isqrt(Integer(48)), 6);
  BOOST_CHECK_EQUAL(isqrt(Integer(49)), 7);
  BOOST_CHECK_EQUAL(isqrt(Integer(50)), 7);

  BOOST_CHECK_EQUAL(isqrt(Integer(100)), 10);
  BOOST_CHECK_EQUAL(isqrt(Integer(143)), 11);

  BOOST_CHECK_EQUAL(isqrt(pow(Integer(105253), 6)), pow(Integer(105253), 3));
  BOOST_CHECK_EQUAL(isqrt(pow(Integer(462711), 10)), pow(Integer(462711), 5));
  BOOST_CHECK_EQUAL(isqrt(pow(Integer(23451), 14)), pow(Integer(23451), 7));
}

BOOST_AUTO_TEST_CASE(MulDivMod)
{
  const Integer small = pow(Integer(1848), 5);
  const Integer large = pow(Integer(140849), 5);

  BOOST_REQUIRE_LT(small, large);
  BOOST_CHECK_EQUAL((small*large)/small, large);
  BOOST_CHECK_EQUAL((small*large)/large, small);
  BOOST_CHECK_EQUAL((large%small) + (large/small)*small, large);
}

BOOST_AUTO_TEST_CASE(AddSub)
{
  const Integer small = pow(Integer(1848), 5);
  const Integer large = pow(Integer(140849), 5);

  BOOST_REQUIRE_LT(small, large);
  BOOST_CHECK_EQUAL(small-small, 0);
  BOOST_CHECK_EQUAL(large-large, 0);
  BOOST_CHECK_EQUAL((small-large)+large, small);
  BOOST_CHECK_EQUAL((large-small)+small, large);
}

BOOST_AUTO_TEST_CASE(IncDec)
{
  const Integer value = pow(Integer(-18484), 5);
  Integer inc(value);
  Integer dec(value);

  ++inc;
  --dec;

  BOOST_CHECK_EQUAL(value+1, inc);
  BOOST_CHECK_EQUAL(value-1, dec);

  Integer incZero(0);
  Integer decZero(0);

  ++incZero;
  --decZero;
  BOOST_CHECK_EQUAL(incZero, 1);
  BOOST_CHECK_EQUAL(decZero, -1);
}

BOOST_AUTO_TEST_CASE(Shifts)
{
  const Integer value = pow(Integer(46272), 7);
  BOOST_CHECK_EQUAL(value << 53, value * pow(Integer(2), 53));
  BOOST_CHECK_EQUAL(value >> 41, value / pow(Integer(2), 41));
}

BOOST_AUTO_TEST_CASE(InputOutput)
{
  const Integer positive = pow(Integer(47383), 19);
  const Integer negative = pow(Integer(-97748), 17);

  BOOST_REQUIRE_GT(positive, 0);
  BOOST_REQUIRE_LT(negative, 0);

  const std::string positiveString = boost::lexical_cast<std::string>(positive);
  const std::string negativeString = boost::lexical_cast<std::string>(negative);

  BOOST_CHECK_EQUAL(positive, Integer(positiveString.c_str()));
  BOOST_CHECK_EQUAL(negative, Integer(negativeString.c_str()));
}
