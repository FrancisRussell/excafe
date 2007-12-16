#ifndef SIMPLE_CFD_STOKES_SYSTEM
#define SIMPLE_CFD_STOKES_SYSTEM

#include <iostream>
#include <map>
#include "mesh.hpp"
#include "simple_cfd_fwd.hpp"
#include "dof_map_builder.hpp"
#include "dof_map.hpp"
#include "lagrange_triangle_linear.hpp"
#include "lagrange_triangle_quadratic.hpp"
#include <boost/tuple/tuple.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace cfd
{

template<typename C>
class stokes_system
{
public:
  typedef C cell_type;

private:
  typedef typename cell_type::vertex_type vertex_type;
  const mesh<cell_type>& m;
  lagrange_triangle_linear pressure;
  lagrange_triangle_quadratic velocity_x;
  lagrange_triangle_quadratic velocity_y;
  dof_map<cell_type> dofMap;
  boost::numeric::ublas::compressed_matrix<double> stiffness_matrix;
  boost::numeric::ublas::vector<double> load_vector;

public:
  stokes_system(const mesh<cell_type>& _m) : m(_m), pressure(m), velocity_x(m), velocity_y(m)
  {  
    dof_map_builder<cell_type> mapBuilder(m);
    mapBuilder.addFiniteElement(pressure);
    mapBuilder.addFiniteElement(velocity_x);
    mapBuilder.addFiniteElement(velocity_y);
    mapBuilder.handleCells(m.getCells());
    dofMap = mapBuilder.getDofMap();
    
    std::cout << "Size of dof map: " << dofMap.getMappingSize() << std::endl;
    std::cout << "Degrees of freedom: " << dofMap.getDegreesOfFreedomCount() << std::endl;
    std::cout << "Degrees of freedom on boundary: " << dofMap.getBoundaryDegreesOfFreedomCount() << std::endl;
  }

  void assemble()
  {
    const std::map<cell_id, cell_type> cells(m.getCells());
    const unsigned dofs = dofMap.getDegreesOfFreedomCount();
    stiffness_matrix.resize(dofs, dofs, false);
    load_vector.resize(dofs, false);

    // Iterate over cells
    for(typename std::map<cell_id, cell_type>::const_iterator cellIter(cells.begin()); cellIter != cells.end(); ++cellIter)
    {
      // Iterate over quadrature points
      const std::map<vertex_type, double> quadrature(cellIter->second.getQuadrature(m.getGeometry()));
      for(typename std::map<vertex_type, double>::const_iterator quadIter(quadrature.begin()); quadIter != quadrature.end(); ++quadIter)
      {
        // Iterate over quadratic test functions
        for(unsigned test=0; test<velocity_x.space_dimension(); ++test)
        {
          const unsigned global_test_x = dofMap.getGlobalIndex(boost::make_tuple(&velocity_x, cellIter->first, test));
          const unsigned global_test_y = dofMap.getGlobalIndex(boost::make_tuple(&velocity_y, cellIter->first, test));
          evaluated_basis evaluated_test = velocity_x.evaluate_basis(cellIter->first, test, quadIter->first);

          // Iterate over quadratic trial functions
          for(unsigned trial=0; trial<velocity_x.space_dimension(); ++trial)
          {
            // assemble velocity related part of momentum equation
            const unsigned global_trial_x = dofMap.getGlobalIndex(boost::make_tuple(&velocity_x, cellIter->first, trial));
            const unsigned global_trial_y = dofMap.getGlobalIndex(boost::make_tuple(&velocity_y, cellIter->first, trial));

            evaluated_basis evaluated_trial = velocity_x.evaluate_basis(cellIter->first, trial, quadIter->first);

            stiffness_matrix(global_test_x, global_trial_x) += quadIter->second * (evaluated_test.dx*evaluated_trial.dx + evaluated_test.dy*evaluated_trial.dy);
            stiffness_matrix(global_test_y, global_trial_y) += quadIter->second * (evaluated_test.dx*evaluated_trial.dx + evaluated_test.dy*evaluated_trial.dy);
          }
          
          // Iterate over linear trial functions
          for(unsigned trial=0; trial<pressure.space_dimension(); ++trial)
          {
            // assemble pressure related part of momentum equation
            const unsigned global_trial = dofMap.getGlobalIndex(boost::make_tuple(&pressure, cellIter->first, trial));
            evaluated_basis evaluated_trial = pressure.evaluate_basis(cellIter->first, trial, quadIter->first);

            stiffness_matrix(global_test_x, global_trial) += quadIter->second * evaluated_test.value*evaluated_trial.dx;
            stiffness_matrix(global_test_y, global_trial) += quadIter->second * evaluated_test.value*evaluated_trial.dy;
          }
        }
          
        // TODO: work out why we're using linear test functions,
        // even though we're dealing with velocity in this formula.
        
        // Iterate over linear test functions 
        for(unsigned test=0; test<pressure.space_dimension(); ++test)
        {
          const unsigned global_test = dofMap.getGlobalIndex(boost::make_tuple(&pressure, cellIter->first, test));
          evaluated_basis evaluated_test = pressure.evaluate_basis(cellIter->first, test, quadIter->first);

          //Iterate over quadratic trial functions
          for(unsigned trial=0; trial<velocity_x.space_dimension(); ++trial)
          {
            // assemble continuity equation
            const unsigned global_trial_x = dofMap.getGlobalIndex(boost::make_tuple(&velocity_x, cellIter->first, trial));
            const unsigned global_trial_y = dofMap.getGlobalIndex(boost::make_tuple(&velocity_y, cellIter->first, trial));
            evaluated_basis evaluated_trial = velocity_x.evaluate_basis(cellIter->first, trial, quadIter->first);

            stiffness_matrix(global_test, global_trial_x) += quadIter->second * evaluated_test.value * evaluated_trial.dx;
            stiffness_matrix(global_test, global_trial_y) += quadIter->second * evaluated_test.value * evaluated_trial.dy;
          }
        }
      }
    }   
  }

  void applyBoundaryConditions()
  {
    const std::map<unsigned, std::set< boost::tuple<cell_id, unsigned> > > pressure_dofs(dofMap.getBoundaryDegreesOfFreedom(&pressure));
    const std::map<unsigned, std::set< boost::tuple<cell_id, unsigned> > > velocity_x_dofs(dofMap.getBoundaryDegreesOfFreedom(&velocity_x));
    const std::map<unsigned, std::set< boost::tuple<cell_id, unsigned> > > velocity_y_dofs(dofMap.getBoundaryDegreesOfFreedom(&velocity_y));
  }

  void print() const
  {
    std::cout << stiffness_matrix << std::endl;
  }
};

}

#endif
