#ifndef EXCAFE_SYMBOLIC_FLAGS_HPP
#define EXCAFE_SYMBOLIC_FLAGS_HPP

namespace excafe
{

namespace symbolic
{

struct Flags
{
  static const unsigned DO_NOT_SIMPLIFY            = 0x1;
  static const unsigned DO_NOT_COLLECT             = 0x2;
  static const unsigned DO_NOT_SIMPLIFY_SUBST_MAP  = 0x4;
};

}

}

#endif
