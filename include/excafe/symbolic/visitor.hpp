#ifndef EXCAFE_SYMBOLIC_VISITOR_HPP
#define EXCAFE_SYMBOLIC_VISITOR_HPP

namespace excafe
{

namespace symbolic
{


// Base class for acyclic visitor pattern
class Visitor {
public:
  virtual void visit(const Basic& b) = 0;
};


}

}

#endif
