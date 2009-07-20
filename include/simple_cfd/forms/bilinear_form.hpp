#ifndef SIMPLE_CFD_FORMS_BILINEAR_FORM_HPP
#define SIMPLE_CFD_FORMS_BILINEAR_FORM_HPP

#include <cstddef>

namespace cfd
{

template<std::size_t D>
class BilinearForm
{
private:
  static const std::size_t dimension = D;
  boost::shared_ptr< Field<dimension> > field;

public:
};

}

#endif
