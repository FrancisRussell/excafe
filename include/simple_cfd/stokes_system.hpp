#ifndef SIMPLE_CFD_STOKES_SYSTEM
#define SIMPLE_CFD_STOKES_SYSTEM

#include "mesh.hpp"
#include "simple_cfd_fwd.hpp"
#include "dof_map_builder.hpp"
#include "dof_map.hpp"
#include "lagrange_triangle_linear.hpp"
#include "lagrange_triangle_quadratic.hpp"
#include <iostream>
#include <map>
#include <algorithm>
#include <ostream>
#include <fstream>
#include <utility>
#include <boost/tuple/tuple.hpp>
#include <mtl/mtl.h>
#include <itl/interface/mtl.h>
#include <itl/krylov/bicgstab.h>
#include <itl/preconditioner/ilu.h>

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
  typedef mtl::matrix< double, mtl::rectangle<>, mtl::compressed<> >::type matrix_type;
  typedef mtl::dense1D<double> vector_type;
  const unsigned dofs;
  matrix_type stiffness_matrix;
  vector_type unknown_vector;
  vector_type load_vector;

  enum Location
  {
    TOP_EDGE,
    BOTTOM_EDGE,
    LEFT_EDGE,
    RIGHT_EDGE,
    BODY
  };

  static dof_map<cell_type> buildDofMap(const mesh<cell_type>& m,
                                        const lagrange_triangle_linear& pressure, 
                                        const lagrange_triangle_quadratic& velocity_x, 
                                        const lagrange_triangle_quadratic& velocity_y)
  {
    dof_map_builder<cell_type> mapBuilder(m);
    mapBuilder.addFiniteElement(pressure);
    mapBuilder.addFiniteElement(velocity_x);
    mapBuilder.addFiniteElement(velocity_y);
    mapBuilder.handleCells(m.getCells());
    return mapBuilder.getDofMap();
  }

  static void zeroRow(matrix_type& m, const unsigned row)
  {
    std::fill(mtl::rows(m)[row].begin(), mtl::rows(m)[row].end(), 0.0);
  }


public:
  stokes_system(const mesh<cell_type>& _m) : m(_m), pressure(m), velocity_x(m), velocity_y(m), 
                                             dofMap(buildDofMap(m, pressure, velocity_x, velocity_y)),
                                             dofs(dofMap.getDegreesOfFreedomCount()), stiffness_matrix(dofs, dofs),
                                             unknown_vector(dofs), load_vector(dofs)
  {  
    std::cout << "Size of dof map: " << dofMap.getMappingSize() << std::endl;
    std::cout << "Degrees of freedom: " << dofMap.getDegreesOfFreedomCount() << std::endl;
    std::cout << "Degrees of freedom on boundary: " << dofMap.getBoundaryDegreesOfFreedomCount() << std::endl;
  }

  void assemble()
  {
    const std::map<cell_id, cell_type> cells(m.getCells());

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


  Location getLocation(const vertex_type& v)
  {
    // TODO: make me into a function that doesn't depend on the specifics of the mesh generation

    Location location = BODY;
    if (v[0] == 0.0)
      location = LEFT_EDGE;

    if (v[0] == 1.0)
      location = RIGHT_EDGE;

    if (v[1] == 0.0)
      location = BOTTOM_EDGE;

    if (v[1] == 1.0)
      location = TOP_EDGE;

    return location;
  }

  void applyBoundaryConditions()
  {
    typedef std::map<unsigned, std::set< boost::tuple<cell_id, unsigned> > > dof_map;
    const dof_map pressure_dofs(dofMap.getBoundaryDegreesOfFreedom(&pressure));
    const dof_map velocity_x_dofs(dofMap.getBoundaryDegreesOfFreedom(&velocity_x));
    const dof_map velocity_y_dofs(dofMap.getBoundaryDegreesOfFreedom(&velocity_y));

    // Assign a value for the pressure degrees around the edges
    for(dof_map::const_iterator pressureIter(pressure_dofs.begin()); pressureIter != pressure_dofs.end(); ++pressureIter)
    {
      const boost::tuple<cell_id, unsigned> dofInfo(*pressureIter->second.begin());
      const vertex_type position(pressure.getDofCoordinate(boost::get<0>(dofInfo), boost::get<1>(dofInfo)));
      const Location location = getLocation(position);

      if (location == BOTTOM_EDGE || location == TOP_EDGE)
      {
        const unsigned pressureGlobalDof(pressureIter->first);
        zeroRow(stiffness_matrix, pressureGlobalDof);
        stiffness_matrix(pressureGlobalDof, pressureGlobalDof) = 1.0;
        load_vector[pressureGlobalDof] = 0.0;

        //To help convergence
        unknown_vector[pressureGlobalDof] = load_vector[pressureGlobalDof];
      }
    }

    // Assign values for the x velocity degrees of freedom (x^2 + y^2 around the entire boundary)
    for(dof_map::const_iterator velocity_x_iter(velocity_x_dofs.begin()); velocity_x_iter!=velocity_x_dofs.end(); ++velocity_x_iter)
    {
      assert(!velocity_x_iter->second.empty());  // If this failed, it would mean a degree of freedom tied to no cell
      const boost::tuple<cell_id, unsigned> dofInfo(*velocity_x_iter->second.begin());
      

      const vertex_type position(velocity_x.getDofCoordinate(boost::get<0>(dofInfo), boost::get<1>(dofInfo)));
      // Check this really is an edge cell
      assert(getLocation(position) != BODY);

      const unsigned velocity_x_globalDof(velocity_x_iter->first);
      zeroRow(stiffness_matrix, velocity_x_globalDof);
      stiffness_matrix(velocity_x_globalDof, velocity_x_globalDof) = 1.0;
      load_vector[velocity_x_globalDof] = position[0] * position[0] + position[1] * position[1]; // x^2 + y^2

      // To help convergence
      unknown_vector[velocity_x_globalDof] = load_vector[velocity_x_globalDof];
    }

    // Assign values for the y velocity degrees of freedom (zero along top and bottom edges) 
    for(dof_map::const_iterator velocity_y_iter(velocity_y_dofs.begin()); velocity_y_iter!=velocity_y_dofs.end(); ++velocity_y_iter)
    {
      assert(!velocity_y_iter->second.empty());  // If this failed, it would mean a degree of freedom tied to no cell
      const boost::tuple<cell_id, unsigned> dofInfo(*velocity_y_iter->second.begin());
      const vertex_type position(velocity_y.getDofCoordinate(boost::get<0>(dofInfo), boost::get<1>(dofInfo)));
      const Location location = getLocation(position);

      // Check this really is an edge cell
      assert(location != BODY);

      if (location == LEFT_EDGE || location == RIGHT_EDGE)
      {
        const unsigned velocity_y_globalDof(velocity_y_iter->first);
        zeroRow(stiffness_matrix, velocity_y_globalDof);
        stiffness_matrix(velocity_y_globalDof, velocity_y_globalDof) = 1.0;
        load_vector[velocity_y_globalDof] = 0.0;
        
        // To help convergence
        unknown_vector[velocity_y_globalDof] = load_vector[velocity_y_globalDof];
      }
    }
  }

  void solve()
  {
    const int max_iter = 2048;
    itl::ILU<matrix_type> precond;
    itl::noisy_iteration<double> iter(load_vector, max_iter, 1e-6);
    itl::bicgstab(stiffness_matrix, unknown_vector, load_vector, precond(), iter);
  }

  void print() const
  {
    std::cout << "Stiffness Matrix: " << std::endl;
    print_all_matrix(stiffness_matrix);
    std::cout << std::endl << "Load Vector: " << std::endl;
    print_vector(load_vector);
    std::cout << std::endl;
  }

  void outputToFile(const std::string& filename)
  {
    std::ofstream outFile(filename.c_str());
    render(10, 10, outFile);
    outFile.close();
  }

  std::pair<double, double> getVelocityVector(const vertex_type& vertex) const
  {
    const std::map<cell_id, cell_type> cells(m.getCells());
    const mesh_geometry<mesh<cell_type>::dimension> geometry(m.getGeometry());
    typename std::map<cell_id, cell_type>::const_iterator cellIter(cells.begin());

    while(cellIter!=cells.end())
    {
      if (cellIter->second.contains(geometry, vertex))
      {
        double xVelocity(0.0), yVelocity(0.0);

        for(unsigned dof=0; dof<velocity_x.space_dimension(); ++dof)
        {
          xVelocity += velocity_x.evaluate_basis(cellIter->first, dof, vertex).value * unknown_vector[dofMap.getGlobalIndex(boost::make_tuple(&velocity_x, cellIter->first, dof))];
          yVelocity += velocity_y.evaluate_basis(cellIter->first, dof, vertex).value * unknown_vector[dofMap.getGlobalIndex(boost::make_tuple(&velocity_y, cellIter->first, dof))];
        }
        return std::make_pair(xVelocity, yVelocity);
      }
      else
      {
        ++cellIter;
      }
    }
    return std::make_pair(0.0, 0.0);
  }

  void render(const unsigned xPoints, const unsigned yPoints, std::ostream& out)
  {
    const double xSpacing = 1.0/(xPoints-1);
    const double ySpacing = 1.0/(yPoints-1);

    out << "# vtk DataFile Version 2.0" << std::endl;
    out << "Simple Stokes Solver" << std::endl;
    out << "ASCII" << std::endl;
    out << "DATASET STRUCTURED_POINTS" << std::endl;
    out << "DIMENSIONS " << xPoints << " " << yPoints << " " << 1 << std::endl;
    out << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << std::endl;
    out << "SPACING " << xSpacing << " " << ySpacing << " " << 0.0 << std::endl;
    out << "POINT_DATA " << xPoints * yPoints << std::endl;
    out << "VECTORS velocity_field DOUBLE" << std::endl;

    for(unsigned yPoint=0; yPoint<yPoints; ++yPoint)
    {
      for(unsigned xPoint=0; xPoint<xPoints; ++xPoint)
      {
        const double x = static_cast<double>(xPoint)/(xPoints-1);
        const double y = static_cast<double>(yPoint)/(yPoints-1);

        const std::pair<double, double> velocity(getVelocityVector(vertex_type(x, y)));
        out << velocity.first << " " << velocity.second << " " << 0.0 << std::endl;
      }
    }
  }
};

}

#endif
