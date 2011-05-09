#ifndef SIMPLE_CFD_CODEGEN_UFC_KERNEL_GENERATOR_HPP
#define SIMPLE_CFD_CODEGEN_UFC_KERNEL_GENERATOR_HPP

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
#include <simple_cfd/numeric/cast.hpp>
#include <ostream>
#include <stack>
#include <vector>
#include <map>
#include <string>
#include <cassert>
#include <sstream>
#include <utility>
#include <boost/foreach.hpp>
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

  Product pow(const int exponent) const
  {
    Product result;
    result.coefficient = cln::expt(coefficient, exponent);

    BOOST_FOREACH(const exp_map_t::value_type& exp, *exponents)
      (*result.exponents)[exp.first] = exp.second*exponent;

    return result;
  }

  Product& operator*=(const Product& p)
  {
    coefficient *= p.coefficient;
    BOOST_FOREACH(const exp_map_t::value_type& exp, *p.exponents)
      (*exponents)[exp.first] += exp.second;

    return *this;
  }

  void write(std::ostream& out) const
  {
    exp_map_t numerators, denominators;

    BOOST_FOREACH(const exp_map_t::value_type& exp, *exponents)
    {
      if (exp.second > 0)
        numerators.insert(exp);
      else
        denominators.insert(exp_map_t::value_type(exp.first, -exp.second));
    }

    bool firstNumerator = true;
    const double coefficientAsDouble = cfd::numeric_cast<double>(coefficient);
    if (coefficientAsDouble != 1.0 || numerators.empty())
    {
      firstNumerator = false;
      out << coefficientAsDouble;
    }
    
    for(exp_map_t::const_iterator numIter = numerators.begin(); 
        numIter != numerators.end();
        ++numIter)
    {
      for(int i=0; i<numIter->second; ++i)
      {
        if (!firstNumerator)
          out << "*";
        else
          firstNumerator = false;

        out << numIter->first;
      }
    }

    if (!denominators.empty())
    {
      out << "/(";
      bool firstDenominator = true;
      
      for(exp_map_t::const_iterator denomIter = denominators.begin(); 
          denomIter != denominators.end();
          ++denomIter)
      {
        for(int i=0; i<denomIter->second; ++i)
        {
          if (!firstDenominator)
            out << "*";
          else
            firstDenominator = false;

          out << denomIter->first;
        }
      }

      out << ")";
    }
  }
};

}

class UFCKernelGenerator : public cse::FactorisedExpressionVisitor<cfd::detail::ScalarPlaceholder>
{
private:
  friend class detail::ScalarPlaceholderNamer;
  typedef util::LazyCopy< std::vector<detail::Product> > sum_t;

  static const std::string factorisedTermPrefix;
  static const std::string resultName;
  static const std::string coefficientsName;
  static const std::string cellName;

  const std::map<Field::expr_ptr, unsigned> fieldIndices;
  const std::map<Scalar::expr_ptr, unsigned> scalarIndices;
  std::ostream& out;
  std::stack<sum_t> stack;
  std::map<cse::PolynomialIndex, std::string> factorisedTermNames;

  void pushProduct(const detail::Product& p)
  {
    sum_t sum;
    sum->push_back(p);
    stack.push(sum);
  }

  detail::Product popProduct()
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

  std::string getName(const cse::PolynomialIndex& i)
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

  std::size_t getFieldIndex(const Field::expr_ptr& f) const
  {
    const std::map<Field::expr_ptr, unsigned>::const_iterator iter = fieldIndices.find(f);
    if (iter != fieldIndices.end())
    {
      return iter->second;
    }
    else
    {
      CFD_EXCEPTION("Unable to locate field index for UFC code generation.");
    }
  }

  std::size_t getScalarIndex(const Scalar::expr_ptr& s) const
  {
    const std::map<Scalar::expr_ptr, unsigned>::const_iterator iter = scalarIndices.find(s);
    if (iter != scalarIndices.end())
    {
      return iter->second;
    }
    else
    {
      CFD_EXCEPTION("Unable to locate scalar index for UFC code generation.");
    }
  }

public:
  UFCKernelGenerator(std::ostream& _out,
                     const std::map<Field::expr_ptr, unsigned>& _fieldIndices,
                     const std::map<Scalar::expr_ptr, unsigned>& _scalarIndices) : 
    fieldIndices(_fieldIndices), scalarIndices(_scalarIndices), out(_out)
  {
  }

  void writeSum()
  {
    assert(stack.size() == 1);
    const sum_t sum = stack.top(); stack.pop();
    
    for(std::vector<detail::Product>::const_iterator prodIter = sum->begin(); prodIter != sum->end(); ++prodIter)
    {
      if (prodIter != sum->begin())
        out << " + ";

      prodIter->write(out);
    }
  }

  void outputPrefix()
  {
    out << "void tabulate_tensor(double* const " << resultName;
    out << ", const double * const * " << coefficientsName;
    out << ", const ufc::cell& " << cellName << ")" << std::endl;
    out << "{" << std::endl;
  }

  void outputPostfix()
  {
    out << "}" << std::endl;
  }

  void visitConstant(const value_t& s)
  {
    pushProduct(s);
  }

  void visitVariable(const variable_t& var);

  void visitExponent(const int exponent)
  {
    const detail::Product p = popProduct();
    pushProduct(p.pow(exponent));
  }

  void postSummation(const std::size_t nops)
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

  void postProduct(const std::size_t nops)
  {
    assert(stack.size() >= nops);

    detail::Product product;
    for(std::size_t i=0; i<nops; ++i)
    {
      const detail::Product subTerm = popProduct();
      product *= subTerm;
    }

    pushProduct(product);
  }

  void visitOriginalTerm(const unsigned index)
  {
    std::ostringstream stream;
    stream << resultName << "[" << index << "]";
    pushProduct(stream.str());
  }

  void visitFactorisedTerm(const cse::PolynomialIndex& index)
  {
    pushProduct(getName(index));
  }

  void postOriginalTerm(const unsigned index)
  {
    out << "  " << resultName << "[" << index << "] = ";
    writeSum();
    out << ";" << std::endl;
  }

  void postFactorisedTerm(const cse::PolynomialIndex& index)
  {
    out << "  const double " << getName(index) << " = ";
    writeSum();
    out << ";" << std::endl;
  }
};

namespace detail
{

class ScalarPlaceholderNamer : public boost::static_visitor<std::string>
{
private:
  const UFCKernelGenerator& generator;

public:
  ScalarPlaceholderNamer(const UFCKernelGenerator& _generator) : generator(_generator)
  {
  }

  result_type operator()(const cfd::detail::PositionComponent& c) const
  {
    CFD_EXCEPTION("PositionComponent detected in cell integral.");
  }

  result_type operator()(const cfd::detail::CellVertexComponent& c) const
  {
    std::ostringstream stream;
    stream << generator.cellName << ".coordinates[" << c.getVertexID() << "][" << c.getComponent() << "]";
    return stream.str();
  }

  result_type operator()(const cfd::detail::ScalarAccess& s) const
  {
    std::ostringstream stream;
    stream << generator.coefficientsName << "[?]";
    return stream.str();
  }

  result_type operator()(const cfd::detail::BasisCoefficient& c) const
  {
    std::ostringstream stream;
    stream << generator.coefficientsName << "[?]";
    return stream.str();
  }

  result_type operator()(const boost::blank&) const
  {
    CFD_EXCEPTION("boost::blank found in ScalarPlaceholder. This should never happen.");
  }
};

}

}

}

#endif
