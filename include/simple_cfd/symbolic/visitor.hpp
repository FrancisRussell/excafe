#ifndef SIMPLE_CFD_SYMBOLIC_VISITOR_HPP
#define SIMPLE_CFD_SYMBOLIC_VISITOR_HPP

namespace cfd
{

namespace symbolic
{


// Base class for acyclic visitor pattern
class Visitor {
public:
  virtual void accept(const Basic& b) = 0;
};


}

}

#endif
