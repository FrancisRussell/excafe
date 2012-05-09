#ifndef EXCAFE_NUMERIC_VALUE_MAP_HPP
#define EXCAFE_NUMERIC_VALUE_MAP_HPP

#include <map>
#include <utility>

namespace excafe
{

namespace detail
{

template<typename V, typename S>
class ValueMap
{
public:
  typedef V variable_type;
  typedef S scalar_type;
  typedef std::map<variable_type, variable_type> var_subst_map;
  typedef std::map<variable_type, scalar_type>   scalar_subst_map;

private:
  var_subst_map    variableSubsts;
  scalar_subst_map scalarSubsts;

public:
  void bind(const variable_type& var, const double s)
  {
    const typename scalar_subst_map::iterator varIter = scalarSubsts.find(var);

    if (varIter != scalarSubsts.end())
      varIter->second = s;
    else
      scalarSubsts.insert(std::make_pair(var, s));

    variableSubsts.erase(var);
  }

  void bind(const variable_type& var, const variable_type& newVar)
  {
    const typename var_subst_map::iterator varIter = variableSubsts.find(var);

    if (varIter != variableSubsts.end())
      varIter->second = newVar;
    else
      variableSubsts.insert(std::make_pair(var, newVar));

    scalarSubsts.erase(var);
  }

  const scalar_subst_map& getScalarSubstitutions() const
  {
    return scalarSubsts;
  }

  const var_subst_map& getVariableSubstitutions() const
  {
    return variableSubsts;
  }
};

}

}

#endif
