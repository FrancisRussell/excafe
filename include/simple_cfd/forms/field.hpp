#ifndef SIMPLE_CFD_FORMS_FIELD_HPP
#define SIMPLE_CFD_FORMS_FIELD_HPP

#include <cstddef>
#include <boost/shared_ptr.hpp>
#include "field_visitor.hpp"

namespace cfd
{

namespace forms
{

class Field
{
public:
  typedef boost::shared_ptr<Field> reference_t;

  virtual std::size_t getRank() const = 0;
  virtual void accept(FieldVisitor& visitor) = 0;
};

}

}

#endif
