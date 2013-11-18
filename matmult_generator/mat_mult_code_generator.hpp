#ifndef MAT_MULT_CODE_GENERATOR_HPP
#define MAT_MULT_CODE_GENERATOR_HPP

#include <string>
#include <iosfwd>
#include <vector>
#include <stack>
#include <excafe/cse/factorised_expression_visitor.hpp>
#include <excafe/cse/polynomial_index.hpp>
#include <excafe/codegen/product.hpp>
#include <excafe/util/lazy_copy.hpp>
#include "vector_entry.hpp"

class MatMultCodeGenerator : public excafe::cse::FactorisedExpressionVisitor<VectorEntry>
{
private:
  typedef excafe::util::LazyCopy< std::vector<excafe::codegen::Product> > sum_t;

  static const std::string factorisedTermPrefix;
  static const std::string resultName;
  static const std::string inputName;
  static const std::string resultNameLocal;
  static const std::string inputNameLocal;
  static const std::string numVectorsName;

  std::ostream& out;
  std::stack<sum_t> stack;
  std::map<excafe::cse::PolynomialIndex, std::string> factorisedTermNames;

  const size_t inputVectorLength;
  const size_t outputVectorLength;

  void pushProduct(const excafe::codegen::Product& p);
  excafe::codegen::Product popProduct();
  std::string getName(const excafe::cse::PolynomialIndex& i);
  static void writeSum(std::ostream& out, const sum_t& sum);

public:
  MatMultCodeGenerator(std::ostream& _out, size_t inputVectorLength, size_t outputVectorLength);
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
  void visitFactorisedTerm(const excafe::cse::PolynomialIndex& index);
  void postOriginalTerm(const unsigned index);
  void postFactorisedTerm(const excafe::cse::PolynomialIndex& index);
};

#endif
