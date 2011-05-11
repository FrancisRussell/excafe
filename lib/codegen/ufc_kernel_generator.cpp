#include <string>
#include <simple_cfd/codegen/ufc_kernel_generator.hpp>

namespace cfd
{

namespace codegen
{

const std::string UFCKernelGenerator::factorisedTermPrefix = "var";
const std::string UFCKernelGenerator::resultName           = "A";
const std::string UFCKernelGenerator::coefficientsName     = "w";
const std::string UFCKernelGenerator::coordinatesName      = "x";

long UFCKernelGenerator::nextClassID = 0;

void UFCKernelGenerator::visitVariable(const variable_t& var)
{
  detail::ScalarPlaceholderNamer namer(*this);
  pushProduct(var.apply(namer));
}

}

}
