#include <excafe/codegen/codegen_fwd.hpp>
#include <excafe/codegen/ufc_integral_generator.hpp>
#include <excafe/codegen/product.hpp>
#include <excafe/exception.hpp>
#include <excafe/util/lazy_copy.hpp>
#include <excafe/capture/assembly/scalar_placeholder.hpp>
#include <excafe/capture/assembly/position_component.hpp>
#include <excafe/capture/assembly/cell_vertex_component.hpp>
#include <excafe/capture/assembly/generic_symbol.hpp>
#include <excafe/capture/assembly/scalar_access.hpp>
#include <excafe/capture/assembly/basis_coefficient.hpp>
#include <excafe/mp/float.hpp>
#include <excafe/cse/factorised_expression_visitor.hpp>
#include <excafe/cse/polynomial_index.hpp>
#include <excafe/numeric/cast.hpp>
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

namespace excafe
{

namespace codegen
{

namespace detail
{

ScalarPlaceholderNamer::ScalarPlaceholderNamer(const UFCIntegralGenerator& _generator) : generator(_generator)
{
}

ScalarPlaceholderNamer::result_type ScalarPlaceholderNamer::operator()(const excafe::detail::PositionComponent& c) const
{
  CFD_EXCEPTION("PositionComponent detected in cell integral.");
}

ScalarPlaceholderNamer::result_type ScalarPlaceholderNamer::operator()(const excafe::detail::CellVertexComponent& c) const
{
  std::ostringstream stream;
  stream << generator.coordinatesName;
  stream << "[" << c.getVertexID() << " * " << generator.geometricDimension << " + " << c.getComponent() << "]";
  return stream.str();
}

ScalarPlaceholderNamer::result_type ScalarPlaceholderNamer::operator()(const excafe::detail::ScalarAccess& s) const
{
  std::ostringstream stream;
  stream << generator.coefficientsName << "[" << generator.getScalarIndex(s.getExpr()) << "][0]";
  return stream.str();
}

ScalarPlaceholderNamer::result_type ScalarPlaceholderNamer::operator()(const excafe::detail::BasisCoefficient& c) const
{
  std::ostringstream stream;
  stream << generator.coefficientsName << "[" << generator.getFieldIndex(c.getField()) << "][" << c.getIndex() << "]";
  return stream.str();
}

ScalarPlaceholderNamer::result_type ScalarPlaceholderNamer::operator()(const excafe::detail::GenericSymbol& s) const
{
  CFD_EXCEPTION("GenericSymbol found in ScalarPlaceholder. This should never happen.");
}


ScalarPlaceholderNamer::result_type ScalarPlaceholderNamer::operator()(const boost::blank&) const
{
  CFD_EXCEPTION("boost::blank found in ScalarPlaceholder. This should never happen.");
}

}

const std::string UFCIntegralGenerator::factorisedTermPrefix = "var";
const std::string UFCIntegralGenerator::resultName           = "A";
const std::string UFCIntegralGenerator::coefficientsName     = "w";
const std::string UFCIntegralGenerator::coordinatesName      = "x";

long UFCIntegralGenerator::nextClassID = 0;

void UFCIntegralGenerator::pushProduct(const Product& p)
{
  sum_t sum;
  sum->push_back(p);
  stack.push(sum);
}

Product UFCIntegralGenerator::popProduct()
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

std::string UFCIntegralGenerator::getName(const cse::PolynomialIndex& i)
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

std::size_t UFCIntegralGenerator::getFieldIndex(const Field::expr_ptr& f) const
{
  const coefficient_index_map_t::const_iterator iter = coefficientIndices.find(f);
  if (iter != coefficientIndices.end())
  {
    return iter->second;
  }
  else
  {
    CFD_EXCEPTION("Unable to locate field index for UFC code generation.");
  }
}

std::size_t UFCIntegralGenerator::getScalarIndex(const Scalar::expr_ptr& s) const
{
  const coefficient_index_map_t::const_iterator iter = coefficientIndices.find(s);
  if (iter != coefficientIndices.end())
  {
    return iter->second;
  }
  else
  {
    CFD_EXCEPTION("Unable to locate scalar index for UFC code generation.");
  }
}

std::string UFCIntegralGenerator::getClassName() const
{
  std::ostringstream stream;
  stream << "ExcafeCellIntegral_" << classID;
  return stream.str();
}

UFCIntegralGenerator::UFCIntegralGenerator(std::ostream& _out,
  const coefficient_index_map_t& _coefficientIndices,
  const int _geometricDimension) :
  classID(nextClassID++),
  coefficientIndices(_coefficientIndices),
  out(_out),
  geometricDimension(_geometricDimension)
{
}

void UFCIntegralGenerator::writeSum(std::ostream& out, const sum_t& sum)
{
  for(std::vector<Product>::const_iterator prodIter = sum->begin(); prodIter != sum->end(); ++prodIter)
  {
    if (prodIter != sum->begin())
      out << " + ";

    prodIter->write(out);
  }
}

void UFCIntegralGenerator::outputPrefix()
{
  out << "class " << getClassName() << " : public ufc::cell_integral\n";
  out << "{\n";
  out << "public:\n";
  out << "  void tabulate_tensor(double *const " << resultName;
  out << ", const double *const *" << coefficientsName;
  out << ", const double *vertex_coordinates";
  out << ", const int cell_orientation) const\n";
  out << "  {\n";
  out << "    const double *" << coordinatesName << " = vertex_coordinates;\n\n";
}

void UFCIntegralGenerator::outputPostfix()
{
  out << "  }\n\n";
  out << "};\n";
}

void UFCIntegralGenerator::visitConstant(const integer_t& s)
{
  pushProduct(float_t(s));
}

void UFCIntegralGenerator::visitConstant(const float_t& s)
{
  pushProduct(s);
}

void UFCIntegralGenerator::visitVariable(const variable_t& var)
{
  detail::ScalarPlaceholderNamer namer(*this);
  pushProduct(var.apply(namer));
}

void UFCIntegralGenerator::visitExponent(const int exponent)
{
  const Product p = popProduct();
  pushProduct(p.pow(exponent));
}

void UFCIntegralGenerator::postSummation(const std::size_t nops)
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

void UFCIntegralGenerator::postProduct(const std::size_t nops)
{
  assert(stack.size() >= nops);

  Product product;
  for(std::size_t i=0; i<nops; ++i)
  {
    const Product subTerm = popProduct();
    product *= subTerm;
  }

  pushProduct(product);
}

void UFCIntegralGenerator::visitOriginalTerm(const unsigned index)
{
  std::ostringstream stream;
  stream << resultName << "[" << index << "]";
  pushProduct(stream.str());
}

void UFCIntegralGenerator::visitFactorisedTerm(const cse::PolynomialIndex& index)
{
  pushProduct(getName(index));
}

void UFCIntegralGenerator::visitAbsoluteValue()
{
  assert(!stack.empty());

  const sum_t sum = stack.top(); stack.pop();
  std::ostringstream stream;
  stream << "std::abs(";
  writeSum(stream, sum);
  stream << ")";

  pushProduct(stream.str());
}

void UFCIntegralGenerator::postOriginalTerm(const unsigned index)
{
  out << "    " << resultName << "[" << index << "] = ";

  assert(stack.size() == 1);
  const sum_t sum = stack.top(); stack.pop();
  writeSum(out, sum);

  out << ";" << std::endl;
}

void UFCIntegralGenerator::postFactorisedTerm(const cse::PolynomialIndex& index)
{
  out << "    const double " << getName(index) << " = ";

  assert(stack.size() == 1);
  const sum_t sum = stack.top(); stack.pop();
  writeSum(out, sum);

  out << ";" << std::endl;
}

}

}
