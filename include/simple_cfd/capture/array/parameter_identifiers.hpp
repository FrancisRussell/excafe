#ifndef SIMPLE_CFD_CAPTURE_ARRAY_PARAMETER_IDENTIFIERS_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_PARAMETER_IDENTIFIERS_HPP

namespace cfd
{

namespace detail
{

class ArrayIndexID
{
private:
  std::size_t offset;

public:
  ArrayIndexID(const std::size_t _offset) : offset(_offset)
  {
  }

  std::size_t getOffset() const
  {
    return offset;
  }
  
  bool operator<(const ArrayIndexID& i) const
  {
    return offset < i.offset;
  }
};

class TensorIndexID
{
private:
  std::size_t offset;

public:
  TensorIndexID(const std::size_t _offset) : offset(_offset)
  {
  }

  std::size_t getOffset() const
  {
    return offset;
  }

  bool operator<(const TensorIndexID& i) const
  {
    return offset < i.offset;
  }
};

}

}

#endif
