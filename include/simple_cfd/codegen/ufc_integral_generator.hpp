#ifndef SIMPLE_CFD_CODEGEN_UFC_INTEGRAL_GENERATOR_HPP
#define SIMPLE_CFD_CODEGEN_UFC_INTEGRAL_GENERATOR_HPP

#include "codegen_fwd.hpp"
#include <simple_cfd/exception.hpp>
#include <simple_cfd/util/lazy_copy.hpp>
#include <simple_cfd/capture/assembly/scalar_placeholder.hpp>
#include <simple_cfd/capture/assembly/position_component.hpp>
#include <simple_cfd/capture/assembly/cell_vertex_component.hpp>
#include <simple_cfd/capture/assembly/scalar_access.hpp>
#include <simple_cfd/capture/assembly/basis_coefficient.hpp>
#include <simple_cfd/cse/factorised_expression_visitor.hpp>
#include <simple_cfd/cse/polynomial_index.hpp>
#include <ostream>
#include <stack>
#include <vector>
#include <map>
#include <string>
#include <boost/variant.hpp>
#include <cln/cln.h>

namespace cfd
{

namespace codegen
{

namespace detail
{

class Product
{
private:
  typedef cln::cl_R numeric_t;
  typedef std::map<std::string, int> exp_map_t;
  typedef util::LazyCopy<exp_map_t> lazy_exp_map_t;

  numeric_t coefficient;
  lazy_exp_map_t exponents;

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
};

}

class UFCIntegralGenerator : public cse::FactorisedExpressionVisitor<cfd::detail::ScalarPlaceholder>
{
public:
  typedef std::map<boost::variant<Field::expr_ptr, Scalar::expr_ptr>, unsigned> coefficient_index_map_t;

private:
  friend class detail::ScalarPlaceholderNamer;
  typedef util::LazyCopy< std::vector<detail::Product> > sum_t;

  static const std::string factorisedTermPrefix;
  static const std::string resultName;
  static const std::string coefficientsName;
  static const std::string coordinatesName;
  static long nextClassID;

  const long classID;
  const coefficient_index_map_t coefficientIndices;
  std::ostream& out;
  std::stack<sum_t> stack;
  std::map<cse::PolynomialIndex, std::string> factorisedTermNames;

  void pushProduct(const detail::Product& p);
  detail::Product popProduct();
  std::string getName(const cse::PolynomialIndex& i);
  std::size_t getFieldIndex(const Field::expr_ptr& f) const;
  std::size_t getScalarIndex(const Scalar::expr_ptr& s) const;

public:
  UFCIntegralGenerator(std::ostream& _out, const coefficient_index_map_t& _coefficientIndices);
  void writeSum();
  void outputPrefix();
  void outputPostfix();
  void visitConstant(const value_t& s);
  void visitVariable(const variable_t& var);
  void visitExponent(const int exponent);
  void postSummation(const std::size_t nops);
  void postProduct(const std::size_t nops);
  void visitOriginalTerm(const unsigned index);
  void visitFactorisedTerm(const cse::PolynomialIndex& index);
  void postOriginalTerm(const unsigned index);
  void postFactorisedTerm(const cse::PolynomialIndex& index);
  std::string getClassName() const;
};

namespace detail
{

class ScalarPlaceholderNamer : public boost::static_visitor<std::string>
{
private:
  const UFCIntegralGenerator& generator;

public:
  ScalarPlaceholderNamer(const UFCIntegralGenerator& _generator);
  result_type operator()(const cfd::detail::PositionComponent& c) const;
  result_type operator()(const cfd::detail::CellVertexComponent& c) const;
  result_type operator()(const cfd::detail::ScalarAccess& s) const;
  result_type operator()(const cfd::detail::BasisCoefficient& c) const;
  result_type operator()(const boost::blank&) const;
};

}

}

}

#endif
