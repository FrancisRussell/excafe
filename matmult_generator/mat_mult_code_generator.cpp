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
const std::string MatMultCodeGenerator::loopIndex("index");


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
  if (generateSSE)
  {
    for(size_t count = 1; count < sum->size(); ++count)
      out << "_mm_add_pd(";

    for(std::vector<codegen::Product>::const_iterator prodIter = sum->begin(); prodIter != sum->end(); ++prodIter)
    {
      if (prodIter != sum->begin())
        out << ", ";

      prodIter->writeSSE(out);

      if (prodIter != sum->begin())
        out << ")";
    }
  }
  else
  {
    for(std::vector<codegen::Product>::const_iterator prodIter = sum->begin(); prodIter != sum->end(); ++prodIter)
    {
      if (prodIter != sum->begin())
        out << " + ";

      prodIter->write(out);
    }
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

  if (generateSSE)
    pushProduct(std::string("_mm_load_pd(&") + stream.str() + ")");
  else
    pushProduct(stream.str());
}

void MatMultCodeGenerator::visitFactorisedTerm(const cse::PolynomialIndex& index)
{
  pushProduct(getName(index));
}

void MatMultCodeGenerator::postOriginalTerm(const unsigned index)
{
  if (generateSSE)
    out << "      _mm_store_pd(&" << resultNameLocal << "[" << index << " * " << numVectorsName << "], ";
  else
    out << "      " << resultNameLocal << "[" << index << " * " << numVectorsName << "] = ";

  assert(stack.size() == 1);
  const sum_t sum = stack.top(); stack.pop();
  writeSum(out, sum);

  if (generateSSE)
    out << ");" << std::endl;
  else
    out << ";" << std::endl;
}

void MatMultCodeGenerator::postFactorisedTerm(const cse::PolynomialIndex& index)
{
  if (generateSSE)
  {
    out << "      const __m128d " << getName(index) << " = ";
  }
  else
  {
    out << "      const double " << getName(index) << " = ";
  }

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

  if (generateSSE)
    nameStream << "_mm_load_pd(&";

  nameStream << inputNameLocal << "[" << var.index << " * " << numVectorsName << "]";

  if (generateSSE)
    nameStream << ")";

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
  out << "#include <stddef.h>" << std::endl;
  out << "#include <stdint.h>" << std::endl;
  out << "#include <stdio.h>" << std::endl;
  out << "#ifdef __SSE2__" << std::endl;
  out << "#include <emmintrin.h>" << std::endl;
  out << "#endif" << std::endl;
  out << std::endl;
  out << "#ifdef __cplusplus" << std::endl;
  out << "extern \"C\"" << std::endl;
  out << "{" << std::endl;
  out << "void frgmv(const double*, double*, size_t);" << std::endl;
  out << "}" << std::endl;
  out << "#endif" << std::endl;
  out << std::endl;
  out << "void frgmv(const double *" << inputName << ", double *" << resultName << ", const size_t " << numVectorsName << ")" << std::endl;
  out << "{" << std::endl;
  out << "  size_t " << loopIndex << " = 0;" << std::endl;
  out << "#ifdef __SSE2__" << std::endl;
  out << "  if ((((uintptr_t) " << inputName << " | (uintptr_t) " << resultName << ") & 15) == 0 && " << numVectorsName << " % 2 == 0)" << std::endl;
  out << "  {" << std::endl;
  out << "    for (; " << loopIndex << " + 1 < " << numVectorsName << "; " << loopIndex << " += 2)" << std::endl;
  out << "    {" << std::endl;
  out << "      const double *" << inputNameLocal << " = " << inputName << " + index;" << std::endl;
  out << "      double *" << resultNameLocal << " = " << resultName << " + index;" << std::endl;
}

void MatMultCodeGenerator::outputInfix()
{
  out << "    }" << std::endl;
  out << "  }" << std::endl;
  out << "#endif" << std::endl;
  out << "  if (" << loopIndex << " < " << numVectorsName << ")" << std::endl;
  out << "    fprintf(stderr, \"Warning: using non-SSE2 implementation\\n\");" << std::endl;
  out << std::endl;
  out << "  for (; " << loopIndex << " < " << numVectorsName << "; ++" << loopIndex << ")" << std::endl;
  out << "  {" << std::endl;
  out << "      const double *" << inputNameLocal << " = " << inputName << " + index;" << std::endl;
  out << "      double *" << resultNameLocal << " = " << resultName << " + index;" << std::endl;

}

void MatMultCodeGenerator::outputPostfix()
{
  out << "  }" << std::endl;
  out << "}" << std::endl;
}
