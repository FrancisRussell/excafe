#include <cassert>
#include <set>
#include <simple_cfd/mp/integer.hpp>
#include <simple_cfd/mp/rational.hpp>
#include <boost/foreach.hpp>

#define BOOST_TEST_MODULE "BigNum"
#include <boost/test/unit_test.hpp>

using namespace cfd::mp;

BOOST_AUTO_TEST_CASE( BasicSanity )
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
  }
}

BOOST_AUTO_TEST_CASE( IntegerSquareRoot )
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


BOOST_AUTO_TEST_CASE( mul_div_mod )
{
  const Integer small = pow(Integer(1848), 5);
  const Integer large = pow(Integer(140849), 5);

  BOOST_REQUIRE_LT(small, large);
  BOOST_CHECK_EQUAL((small*large)/small, large);
  BOOST_CHECK_EQUAL((small*large)/large, small);
  BOOST_CHECK_EQUAL((large%small) + (large/small)*small, large);
}

BOOST_AUTO_TEST_CASE( add_sub )
{
  const Integer small = pow(Integer(1848), 5);
  const Integer large = pow(Integer(140849), 5);

  BOOST_REQUIRE_LT(small, large);
  BOOST_CHECK_EQUAL(small-small, 0);
  BOOST_CHECK_EQUAL(large-large, 0);
  BOOST_CHECK_EQUAL((small-large)+large, small);
  BOOST_CHECK_EQUAL((large-small)+small, large);
}

BOOST_AUTO_TEST_CASE( increments )
{
  const Integer value = pow(Integer(-18484), 5);
  Integer inc(value);
  Integer dec(value);

  ++inc;
  --dec;

  BOOST_CHECK_EQUAL(value+1, inc);
  BOOST_CHECK_EQUAL(value-1, dec);
}
