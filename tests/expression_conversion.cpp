#include <iostream>
#include <string>
#include <simple_cfd/numeric/convert_expression.hpp>
#include <simple_cfd/numeric/polynomial_fraction.hpp>
#include <simple_cfd/numeric/ginac_expression.hpp>
#include <simple_cfd/numeric/excafe_expression.hpp>


template<typename T>
struct TypeNames
{
};

template<typename V>
struct TypeNames< cfd::PolynomialFraction<V> >
{
  static std::string name() { return "PolynomialFraction"; }
};

template<typename V>
struct TypeNames< cfd::GinacExpression<V> >
{
  static std::string name() { return "GinacExpression"; }
};

template<typename V>
struct TypeNames< cfd::ExcafeExpression<V> >
{
  static std::string name() { return "ExcafeExpression"; }
};

template<typename from, typename to>
void printConversion(const from& f)
{
  std::cout << f << " -> " << cfd::detail::convert_expression<to>(f) << std::endl;
}

template<typename from, typename to>
void printConversions()
{
  const std::string k("k");
  const std::string x("x");
  const std::string y("y");

  std::cout << "Converting from " << TypeNames<from>::name() << " to " << TypeNames<to>::name() << ":" << std::endl;
  printConversion<from, to>(from(7.5));
  printConversion<from, to>(from(x) + from(y));
  printConversion<from, to>(from(7, k, 5));
  printConversion<from, to>(from(4, k, 2) + from(2, x, 3));
  printConversion<from, to>(from(2, k, 2) / from(6, x, 3));
  printConversion<from, to>(from(4, k, 2) * from(1, x, 2));
  std::cout << std::endl;
}

int main(int, char**)
{
  using namespace cfd;

  printConversions< PolynomialFraction<std::string>, PolynomialFraction<std::string> >();
  printConversions< PolynomialFraction<std::string>, GinacExpression<std::string> >();
  printConversions< PolynomialFraction<std::string>, ExcafeExpression<std::string> >();

  printConversions< GinacExpression<std::string>, PolynomialFraction<std::string> >();
  printConversions< GinacExpression<std::string>, GinacExpression<std::string> >();
  printConversions< GinacExpression<std::string>, ExcafeExpression<std::string> >();

  printConversions< ExcafeExpression<std::string>, PolynomialFraction<std::string> >();
  printConversions< ExcafeExpression<std::string>, GinacExpression<std::string> >();
  printConversions< ExcafeExpression<std::string>, ExcafeExpression<std::string> >();
}
