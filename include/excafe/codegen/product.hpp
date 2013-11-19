#ifndef EXCAFE_CODEGEN_PRODUCT_HPP
#define EXCAFE_CODEGEN_PRODUCT_HPP

#include "codegen_fwd.hpp"
#include <excafe/util/lazy_copy.hpp>
#include <excafe/mp/float.hpp>
#include <ostream>
#include <map>
#include <string>

namespace excafe
{

namespace codegen
{

class Product
{
private:
  typedef mp::Float numeric_t;
  typedef std::map<std::string, int> exp_map_t;
  typedef util::LazyCopy<exp_map_t> lazy_exp_map_t;

  numeric_t coefficient;
  lazy_exp_map_t exponents;

  void write(std::ostream& out, bool sse) const;
  std::string negate(const std::string& exp, bool sse) const;
  std::string constructPositiveProduct(const exp_map_t& exps, bool sse) const;

public:
  Product() : coefficient(1.0)
  {
  }

  Product(const std::string& v) : coefficient(1.0)
  {
    ++(*exponents)[v];
  }

  Product(const numeric_t& v) : coefficient(v)
  {
  }

  Product pow(const int exponent) const;
  Product& operator*=(const Product& p);
  void write(std::ostream& out) const;
  void writeSSE(std::ostream& out) const;
};

}

}


#endif
