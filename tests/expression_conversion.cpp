#include <iostream>
#include <string>
#include <simple_cfd/numeric/convert_expression.hpp>
#include <simple_cfd/numeric/polynomial_fraction.hpp>
#include <simple_cfd/numeric/ginac_expression.hpp>

template<typename from, typename to>
void test(const from& f)
{
  std::cout << f << " -> " << cfd::detail::convert_expression<to>(f) << std::endl;
}

template<typename from, typename to>
void test()
{
  test<from, to>(from(7.5));
  test<from, to>(from("x") + from("y"));
  test<from, to>(from(7, "k", 5));
  test<from, to>(from(4, "k", 2) + from(2, "x", 3));
  test<from, to>(from(2, "k", 2) / from(6, "x", 3));
  test<from, to>(from(4, "k", 2) * from(1, "x", 2));
}


int main(int, char**)
{
  using namespace cfd;

  std::cout << "Converting from PolynomialFraction to GinacExpression:" << std::endl; 
  test< PolynomialFraction<std::string>, GinacExpression<std::string> >();

  std::cout << std::endl;

  std::cout << "Converting from GinacExpression to PolynomialFraction:" << std::endl; 
  test< GinacExpression<std::string>, PolynomialFraction<std::string> >();
}
