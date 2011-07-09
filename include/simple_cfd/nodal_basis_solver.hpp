#ifndef SIMPLE_CFD_NODAL_BASIS_SOLVER_HPP
#define SIMPLE_CFD_NODAL_BASIS_SOLVER_HPP

#include "mp/rational.hpp"
#include "vertex.hpp"
#include "exception.hpp"
#include "numeric/excafe_expression.hpp"
#include "capture/assembly/scalar_placeholder.hpp"
#include "capture/assembly/position_placeholder.hpp"
#include <vector>
#include <cassert>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/optional.hpp>

namespace cfd
{

template<std::size_t D>
class NodalBasisSolver
{
public:
  static const std::size_t dimension = D;

private:
  typedef boost::numeric::ublas::matrix<mp::Rational> matrix_t;

  std::vector< vertex<dimension, mp::Rational> > points;
  std::vector< ExcafeExpression<detail::ScalarPlaceholder> > primeBases;
  std::size_t spaceDimension;

  static matrix_t invert(const matrix_t& matrix)
  {
    using namespace boost::numeric::ublas;

    matrix_t input(matrix);
    permutation_matrix<std::size_t> permutationMatrix(input.size1());

    const int result = lu_factorize(input, permutationMatrix);
    if (result != 0)
      CFD_EXCEPTION("LU factorisation failure.");

    matrix_t inverse;
    inverse = identity_matrix<matrix_t::value_type>(input.size1());
    lu_substitute(input, permutationMatrix, inverse);

    return inverse;
  }

public:
  NodalBasisSolver(const std::vector< vertex<dimension, mp::Rational> >& _points,
                   const std::vector< ExcafeExpression<detail::ScalarPlaceholder> >& _primeBases) :
    points(_points), primeBases(_primeBases), spaceDimension(points.size())
  {
    assert(primeBases.size() == spaceDimension);
    assert(points.size() == spaceDimension);
  }

  std::vector< ExcafeExpression<detail::ScalarPlaceholder> > getBases()
  {
    using detail::ScalarPlaceholder;
    using mp::Rational;

    matrix_t evaluatedBases(spaceDimension, spaceDimension);

    for(std::size_t pointID=0; pointID < spaceDimension; ++pointID)
    {
      const vertex<dimension, Rational> point = points[pointID];

      ExcafeExpression<ScalarPlaceholder>::value_map componentValues;
      detail::PositionPlaceholder position;

      for(std::size_t d=0; d<dimension; ++d)
        componentValues.bind(position[d], point[d]);

      for(std::size_t basisID=0; basisID < spaceDimension; ++basisID)
      {
        const ExcafeExpression<ScalarPlaceholder> basis = primeBases[basisID];
        const boost::optional<Rational> basisValue = basis.evaluateRational(componentValues);

        if (basisValue)
          evaluatedBases(basisID, pointID) = *basisValue;
        else
          CFD_EXCEPTION("Couldn't evaluate basis function as a rational value.");
      }
    }

    const matrix_t coefficients = invert(evaluatedBases);
    std::vector< ExcafeExpression<ScalarPlaceholder> > bases;

    for(std::size_t basisID = 0; basisID < spaceDimension; ++basisID)
    {
      ExcafeExpression<ScalarPlaceholder> basis;

      for(std::size_t primeID=0; primeID < spaceDimension; ++primeID)
        basis += coefficients(basisID, primeID) * primeBases[primeID];

      bases.push_back(basis);
    }

    return bases;
  }
};

}
#endif
