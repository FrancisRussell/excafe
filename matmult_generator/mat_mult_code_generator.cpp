#include <excafe/cse/polynomial_index.hpp>
#include <excafe/util/lazy_copy.hpp>
#include <excafe/exception.hpp>
#include "mat_mult_code_generator.hpp"
#include <ostream>
#include <sstream>
#include <string>
#include <map>
#include <cassert>

using namespace excafe;

const std::string MatMultCodeGenerator::factorisedTermPrefix("subterm");
const std::string MatMultCodeGenerator::resultNameLocal("b_local");
const std::string MatMultCodeGenerator::inputNameLocal("x_local");
const std::string MatMultCodeGenerator::numVectorsName("count");

const std::string MatMultCodeGenerator::resultName("b");
const std::string MatMultCodeGenerator::inputName("x");

void MatMultCodeGenerator::pushProduct(const codegen::Product& p)
{
  sum_t sum;
  sum->push_back(p);
  stack.push(sum);
}

codegen::Product MatMultCodeGenerator::popProduct()
{
  assert(!stack.empty());

  const sum_t sum = stack.top();
  stack.pop();

  if (sum->size() != 1)
  {
    CFD_EXCEPTION("Tried to pop a product, got a sum of products.");
  }
  else
  {
    return sum->front();
  }
}

std::string MatMultCodeGenerator::getName(const cse::PolynomialIndex& i)
{
  const std::map<cse::PolynomialIndex, std::string>::const_iterator iter = factorisedTermNames.find(i);

  if (iter != factorisedTermNames.end())
  {
    return iter->second;
  }
  else
  {
    std::ostringstream stream;
    stream << factorisedTermPrefix << "_" << i.getIndex();
    const std::string name = stream.str();

    factorisedTermNames.insert(std::make_pair(i, name));
    return name;
  }
}

void MatMultCodeGenerator::writeSum(std::ostream& out, const sum_t& sum)
{
  for(std::vector<codegen::Product>::const_iterator prodIter = sum->begin(); prodIter != sum->end(); ++prodIter)
  {
    if (prodIter != sum->begin())
      out << " + ";

    prodIter->write(out);
  }
}

MatMultCodeGenerator::MatMultCodeGenerator(
  std::ostream& _out, const size_t _inputVectorLength, const size_t _outputVectorLength) :
  out(_out), inputVectorLength(_inputVectorLength), outputVectorLength(_outputVectorLength)
{
}

void MatMultCodeGenerator::visitOriginalTerm(const unsigned index)
{
  std::ostringstream stream;
  stream << resultNameLocal << "[" << index << " * " << numVectorsName << "]";
  pushProduct(stream.str());
}

void MatMultCodeGenerator::visitFactorisedTerm(const cse::PolynomialIndex& index)
{
  pushProduct(getName(index));
}

void MatMultCodeGenerator::postOriginalTerm(const unsigned index)
{
  out << "    " << resultNameLocal << "[" << index << " * " << numVectorsName << "] = ";

  assert(stack.size() == 1);
  const sum_t sum = stack.top(); stack.pop();
  writeSum(out, sum);

  out << ";" << std::endl;
}

void MatMultCodeGenerator::postFactorisedTerm(const cse::PolynomialIndex& index)
{
  out << "    const double " << getName(index) << " = ";

  assert(stack.size() == 1);
  const sum_t sum = stack.top(); stack.pop();
  writeSum(out, sum);

  out << ";" << std::endl;
}

void MatMultCodeGenerator::visitConstant(const float_t& s)
{
  pushProduct(s);
}

void MatMultCodeGenerator::visitConstant(const integer_t& s)
{
  pushProduct(float_t(s));
}

void MatMultCodeGenerator::visitVariable(const variable_t& var)
{
  std::ostringstream nameStream;
  nameStream << inputNameLocal << "[" << var.index << " * " << numVectorsName << "]";
  pushProduct(nameStream.str());
}

void MatMultCodeGenerator::visitExponent(const int exponent)
{
  const codegen::Product p = popProduct();
  pushProduct(p.pow(exponent));
}

void MatMultCodeGenerator::visitAbsoluteValue()
{
  assert(!stack.empty());

  const sum_t sum = stack.top(); stack.pop();
  std::ostringstream stream;
  stream << "std::abs(";
  writeSum(stream, sum);
  stream << ")";

  pushProduct(stream.str());
}

void MatMultCodeGenerator::postSummation(const std::size_t nops)
{
  assert(stack.size() >= nops);

  sum_t sum;
  for(std::size_t i=0; i<nops; ++i)
  {
    const sum_t subTerm = stack.top();
    stack.pop();
    sum->insert(sum->end(), subTerm->begin(), subTerm->end());
  }

  stack.push(sum);
}

void MatMultCodeGenerator::postProduct(const std::size_t nops)
{
  assert(stack.size() >= nops);

  codegen::Product product;
  for(std::size_t i=0; i<nops; ++i)
  {
    const codegen::Product subTerm = popProduct();
    product *= subTerm;
  }

  pushProduct(product);
}

void MatMultCodeGenerator::outputPrefix()
{
  const std::string loopIndex = "index";

  out << "#include <stddef.h>" << std::endl;
  out << std::endl;
  out << "void frgmv(const double *" << inputName << ", double *" << resultName << ", const size_t " << numVectorsName << ")" << std::endl;
  out << "{" << std::endl;
  out << "  size_t " << loopIndex << " = 0;" << std::endl;
  out << "  for (; " << loopIndex << " < count; ++" << loopIndex << ")" << std::endl;
  out << "  {" << std::endl;
  out << "    const double *" << inputNameLocal << " = " << inputName << " + index;" << std::endl;
  out << "    double *" << resultNameLocal << " = " << resultName << " + index;" << std::endl;
}

void MatMultCodeGenerator::outputPostfix()
{
  out << "  }" << std::endl;
  out << "}" << std::endl;
}
