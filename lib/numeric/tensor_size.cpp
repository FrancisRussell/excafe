#include <excafe/numeric/tensor_size.hpp>
#include <ostream>

namespace std
{

std::ostream& operator<<(std::ostream& o, const excafe::TensorSize& s)
{
  s.write(o);
  return o;
}

}
