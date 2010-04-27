#ifndef SIMPLE_CFD_LOCAL_ASSEMBLY_MATRIX_HPP
#define SIMPLE_CFD_LOCAL_ASSEMBLY_MATRIX_HPP

#include <cassert>
#include <cstddef>
#include <set>
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/operators.hpp>

namespace cfd
{

namespace detail
{

template<std::size_t D, typename T>
class LocalAssemblyMatrix : boost::multipliable<LocalAssemblyMatrix<D,T>, T,
                            boost::dividable<LocalAssemblyMatrix<D,T>, T
                            > >
{
public:
  static const std::size_t dimension = D;
  typedef T value_type;
  typedef const FiniteElement<dimension> element_t;
  typedef typename std::vector<value_type>::iterator iterator;
  typedef typename std::vector<value_type>::const_iterator const_iterator;

private:
  std::set<const element_t*> trialElements;
  std::set<const element_t*> testElements;
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
        offset += element->spaceDimension();
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

public:
  LocalAssemblyMatrix(const std::set<const element_t*>& _trialElements, const std::set<const element_t*> _testElements) : 
    trialElements(_trialElements), testElements(_testElements),
    values(getTrialSize()*getTestSize())
  {
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
};

}

}

#endif
