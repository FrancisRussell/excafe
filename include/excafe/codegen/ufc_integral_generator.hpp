#ifndef EXCAFE_CODEGEN_UFC_INTEGRAL_GENERATOR_HPP
#define EXCAFE_CODEGEN_UFC_INTEGRAL_GENERATOR_HPP

#include "codegen_fwd.hpp"
#include "product.hpp"
#include <excafe/exception.hpp>
#include <excafe/util/lazy_copy.hpp>
#include <excafe/capture/assembly/scalar_placeholder.hpp>
#include <excafe/capture/assembly/position_component.hpp>
#include <excafe/capture/assembly/cell_vertex_component.hpp>
#include <excafe/capture/assembly/scalar_access.hpp>
#include <excafe/capture/assembly/basis_coefficient.hpp>
#include <excafe/capture/assembly/generic_symbol.hpp>
#include <excafe/cse/factorised_expression_visitor.hpp>
#include <excafe/cse/polynomial_index.hpp>
#include <excafe/mp/float.hpp>
#include <ostream>
#include <stack>
#include <vector>
#include <map>
#include <string>
#include <boost/variant.hpp>

namespace excafe
{

namespace codegen
{

class UFCIntegralGenerator : public cse::FactorisedExpressionVisitor<excafe::detail::ScalarPlaceholder>
{
public:
  typedef std::map<boost::variant<Field::expr_ptr, Scalar::expr_ptr>, unsigned> coefficient_index_map_t;

private:
  friend class detail::ScalarPlaceholderNamer;
  typedef util::LazyCopy< std::vector<Product> > sum_t;

  static const std::string factorisedTermPrefix;
  static const std::string resultName;
  static const std::string coefficientsName;
  static const std::string coordinatesName;
  static long nextClassID;

  const long classID;
  const coefficient_index_map_t coefficientIndices;
  const int geometricDimension;
  std::ostream& out;
  std::stack<sum_t> stack;
  std::map<cse::PolynomialIndex, std::string> factorisedTermNames;

  void pushProduct(const Product& p);
  Product popProduct();
  std::string getName(const cse::PolynomialIndex& i);
  std::size_t getFieldIndex(const Field::expr_ptr& f) const;
  std::size_t getScalarIndex(const Scalar::expr_ptr& s) const;
  static void writeSum(std::ostream& out, const sum_t& sum);

public:
  UFCIntegralGenerator(std::ostream& _out, const coefficient_index_map_t& _coefficientIndices, int geometricDimension);
  void outputPrefix();
  void outputPostfix();
  void visitConstant(const float_t& s);
  void visitConstant(const integer_t& s);
  void visitVariable(const variable_t& var);
  void visitExponent(const int exponent);
  void visitAbsoluteValue();
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
  result_type operator()(const excafe::detail::PositionComponent& c) const;
  result_type operator()(const excafe::detail::CellVertexComponent& c) const;
  result_type operator()(const excafe::detail::ScalarAccess& s) const;
  result_type operator()(const excafe::detail::BasisCoefficient& c) const;
  result_type operator()(const excafe::detail::GenericSymbol& s) const;
  result_type operator()(const boost::blank&) const;
};

}

}

}

#endif
