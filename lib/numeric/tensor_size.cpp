#include <simple_cfd/numeric/tensor_size.hpp>
#include <ostream>

namespace std
{

std::ostream& operator<<(std::ostream& o, const cfd::TensorSize& s)
{
  s.write(o);
  return o;
}

}
