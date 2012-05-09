#ifndef EXCAFE_LOCAL_ASSEMBLY_MATRIX_HPP
#define EXCAFE_LOCAL_ASSEMBLY_MATRIX_HPP

#include <cassert>
#include <cstddef>
#include <set>
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/operators.hpp>
#include <ostream>

namespace excafe
{

namespace detail
{

template<std::size_t D, typename T>
class LocalAssemblyMatrix : boost::multiplicative<LocalAssemblyMatrix<D,T>, T>
{
public:
  static const std::size_t dimension = D;
  typedef T value_type;
  typedef const FiniteElement<dimension> element_t;
  typedef typename std::vector<value_type>::iterator iterator;
  typedef typename std::vector<value_type>::const_iterator const_iterator;

private:
  std::set<const element_t*> testElements;
  std::set<const element_t*> trialElements;
  std::vector<value_type> values;

  static std::size_t getSpaceDimension(const std::set<const element_t*>& elements)
  {
    std::size_t size = 0;
    BOOST_FOREACH(const element_t* element, elements)
    {
      size += element->spaceDimension();
    }
    return size;
  }

  std::size_t getOffset(const std::set<const element_t*>& elements, const element_t* element, const std::size_t index) const
  {
    assert(index < element->spaceDimension());

    std::size_t offset = 0;
    BOOST_FOREACH(const element_t* currentElement, elements)
    {
      if (element == currentElement)
      {
        offset += index;
        return offset;
      }
      else
      {
        offset += currentElement->spaceDimension();
      }
    }

    CFD_EXCEPTION("Failed to locate requested element in local assembly matrix");
  }

  std::vector< Dof<dimension> > getDofs(const std::set<const element_t*>& elements, const std::size_t cid) const
  {
    std::vector< Dof<dimension> > dofs;

    BOOST_FOREACH(const element_t* element, elements)
    {
      const std::size_t spaceDimension = element->spaceDimension();
      for (std::size_t index=0; index<spaceDimension; ++index)
      {
        dofs.push_back(Dof<dimension>(element, cid, index));
      }
    }
    return dofs;
  }

  void write(std::ostream& o, const element_t* element) const
  {
    o << "Finite Elemement: " << "address=" << element << ", rank=" << element->getRank();
    o << ", dimension=" << element->getDimension() << ", space_dimension=" << element->spaceDimension();
    o << std::endl;
  }

  void write(std::ostream& o, const std::set<const element_t*>& elements) const
  {
    BOOST_FOREACH(const element_t* element, elements)
    {
      write(o, element);
    }
  }

public:
  LocalAssemblyMatrix(const std::set<const element_t*>& _testElements, const std::set<const element_t*> _trialElements) : 
    testElements(_testElements), trialElements(_trialElements),
    values(getTrialSize()*getTestSize())
  {
  }
  
  void clear()
  {
    std::fill(begin(), end(), value_type());
  }

  const value_type* data() const
  {
    return &values[0];
  }

  value_type* data()
  {
    return &values[0];
  }

  iterator begin()
  {
    return values.begin();
  }

  const_iterator begin() const
  {
    return values.begin();
  }

  iterator end()
  {
    return values.end();
  }

  const_iterator end() const
  {
    return values.end();
  }

  std::size_t getTrialSize() const
  {
    return getSpaceDimension(trialElements);
  }

  std::size_t getTestSize() const
  {
    return getSpaceDimension(testElements);
  }

  std::size_t getTrialOffset(const FiniteElement<dimension>& trialElement, const std::size_t index) const
  {
    return getOffset(trialElements, &trialElement, index);
  }

  std::size_t getTestOffset(const FiniteElement<dimension>& testElement, const std::size_t index) const
  {
    return getOffset(testElements, &testElement, index);
  }

  std::vector< Dof<dimension> > getTrialDofs(const std::size_t cid) const
  {
    return getDofs(trialElements, cid);
  }

  std::vector< Dof<dimension> > getTestDofs(const std::size_t cid) const
  {
    return getDofs(testElements, cid);
  }

  LocalAssemblyMatrix& operator*=(const value_type& v)
  {
    BOOST_FOREACH(value_type& val, values)
    {
      val *= v;
    }
    return *this;
  }

  LocalAssemblyMatrix& operator/=(const value_type& v)
  {
    BOOST_FOREACH(value_type& val, values)
    {
      val /= v;
    }
    return *this;
  }

  value_type& operator()(const std::size_t test, const std::size_t trial)
  {
    const std::size_t testSize = getTestSize();
    const std::size_t trialSize = getTrialSize();
    assert(test < testSize);
    assert(trial < trialSize);

    return values[test * trialSize + trial];
  }

  const value_type& operator()(const std::size_t test, const std::size_t trial) const
  {
    const std::size_t testSize = getTestSize();
    const std::size_t trialSize = getTrialSize();
    assert(test < testSize);
    assert(trial < trialSize);

    return values[test * trialSize + trial];
  }

  template<typename unary_function>
  LocalAssemblyMatrix<dimension, typename unary_function::result_type> transform(const unary_function& f) const
  {
    LocalAssemblyMatrix<dimension, typename unary_function::result_type> result(trialElements, testElements);
    std::transform(begin(), end(), result.begin(), f);
    return result;
  }

  void write(std::ostream& o) const
  {
    o << "Local assembly matrix: size(trials)=" << getTrialSize() << ", size(test)=" << getTestSize() << std::endl;
    o << "Trial functions:" << std::endl;
    write(o, trialElements);
    o << "Test functions:" << std::endl;
    write(o, testElements);
    o << std::endl;

    for(std::size_t test=0; test<getTestSize(); ++test)
    {
      for(std::size_t trial=0; trial<getTrialSize(); ++trial)
      {
        o << "(" << test << ", " << trial << ") = " << (*this)(test, trial) << std::endl;
      }
    }
    o << std::endl;
  }
};

}

}

namespace std
{

template<std::size_t D, typename T>
std::ostream& operator<<(std::ostream& o, const excafe::detail::LocalAssemblyMatrix<D,T>& m)
{
  m.write(o);
  return o;
}

}

#endif
