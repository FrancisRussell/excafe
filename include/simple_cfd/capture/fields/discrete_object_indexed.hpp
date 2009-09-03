#ifndef SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_OBJECT_INDEXED_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_OBJECT_INDEXED_HPP

namespace cfd
{

namespace detail
{

template<typename discrete_object_tag>
class DiscreteObjectIndexed : public DiscreteTraits<discrete_object_tag>::expr_t
{
private:
 typedef IndexableValue<discrete_object_tag> parent_t; 
 typedef typename parent_t::value_ptr parent_ptr;

public:
  DiscreteObjectIndexed(const parent_ptr& _parent, const TemporalIndexExpr indexExpr)
  {
  }

  virtual void accept(DiscreteExprVisitor& visitor)
  {
    visitor.visit(*this);
  }
};

}

}

#endif
