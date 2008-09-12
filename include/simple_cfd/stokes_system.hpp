#ifndef SIMPLE_CFD_STOKES_SYSTEM
#define SIMPLE_CFD_STOKES_SYSTEM

#include "mesh.hpp"
#include "simple_cfd_fwd.hpp"
#include "dof_map_builder.hpp"
#include "dof_map.hpp"
#include "lagrange_triangle_linear.hpp"
#include "lagrange_triangle_quadratic.hpp"
#include "numeric/solver.hpp"
#include "numeric/tensor.hpp"
#include "fe_matrix.hpp"
#include "fe_vector.hpp"
#include "fe_binary_function.hpp"
#include <iostream>
#include <map>
#include <algorithm>
#include <ostream>
#include <fstream>
#include <utility>
#include <boost/tuple/tuple.hpp>

namespace cfd
{

template<typename C>
class GradTrialInnerGradTest : public FEBinaryFunction<C>
{
private:
  typedef typename FEBinaryFunction<C>::cell_type cell_type;
  typedef typename FEBinaryFunction<C>::vertex_type vertex_type;
  typedef typename FEBinaryFunction<C>::finite_element_t finite_element_t;

  const lagrange_triangle_quadratic<1>* const trial;
  const lagrange_triangle_quadratic<1>* const test;

public:
  GradTrialInnerGradTest(const lagrange_triangle_quadratic<1>* const _trial, const lagrange_triangle_quadratic<1>* const _test) : trial(_trial), test(_test)
  {
  }

  virtual const finite_element_t* getTestFunction() const
  {
    return test;
  }

  virtual const finite_element_t* getTrialFunction() const 
  {
    return trial;
  }

  virtual double evaluate(const cell_type& cell, const std::size_t testDof, const std::size_t trialDof, const vertex_type& location) const
  {
    static const unsigned dimension = cell_type::dimension;

    Tensor<dimension, 2, double> trial_gradient = trial->evaluate_gradient(cell, trialDof, location);
    Tensor<dimension, 2, double> test_gradient = test->evaluate_gradient(cell, testDof, location);
    const double result = (test_gradient.colon_product(trial_gradient)).toScalar();
    return result;
  }
};

template<typename TrialType, typename TestType>
class NegTrialInnerDivTest : public FEBinaryFunction<typename TrialType::cell_type>
{
private:
  typedef typename TrialType::cell_type cell_type;
  typedef typename FEBinaryFunction<cell_type>::vertex_type vertex_type;
  typedef typename FEBinaryFunction<cell_type>::finite_element_t finite_element_t;

  const TrialType* const trial;
  const TestType* const test;

public:
  NegTrialInnerDivTest(const TrialType* const _trial, const TestType* const _test) : trial(_trial), test(_test)
  {
  }

  virtual const finite_element_t* getTestFunction() const
  {
    return test;
  }

  virtual const finite_element_t* getTrialFunction() const 
  {
    return trial;
  }

  virtual double evaluate(const cell_type& cell, const std::size_t testDof, const std::size_t trialDof, const vertex_type& location) const
  {
    typename TrialType::value_type trial_value = trial->evaluate_tensor(cell, trialDof, location);
    typename TestType::divergence_type test_divergence = test->evaluate_divergence(cell, testDof, location);
    const double result = - (trial_value * test_divergence).toScalar();
    return result;
  }
};

template<typename TrialType, typename TestType>
class DivTrialInnerTest : public FEBinaryFunction<typename TrialType::cell_type>
{
private:
  typedef typename TrialType::cell_type cell_type;
  typedef typename FEBinaryFunction<cell_type>::vertex_type vertex_type;
  typedef typename FEBinaryFunction<cell_type>::finite_element_t finite_element_t;

  const TrialType* const trial;
  const TestType* const test;

public:
  DivTrialInnerTest(const TrialType* const _trial, const TestType* const _test) : trial(_trial), test(_test)
  {
  }

  virtual const finite_element_t* getTestFunction() const
  {
    return test;
  }

  virtual const finite_element_t* getTrialFunction() const 
  {
    return trial;
  }

  virtual double evaluate(const cell_type& cell, const std::size_t testDof, const std::size_t trialDof, const vertex_type& location) const
  {
    typename TrialType::divergence_type trial_divergence = trial->evaluate_divergence(cell, trialDof, location);
    typename TestType::value_type test_value = test->evaluate_tensor(cell, testDof, location);
    const double result = (test_value * trial_divergence).toScalar();
    return result;
  }
};


template<typename C>
class stokes_system
{
public:
  typedef C cell_type;

private:
  static const unsigned dimension = cell_type::dimension;
  typedef typename cell_type::vertex_type vertex_type;
  typedef lagrange_triangle_linear<0> pressure_basis_t;
  typedef lagrange_triangle_quadratic<1> velocity_basis_t;

  const mesh<cell_type>& m;
  pressure_basis_t pressure;
  velocity_basis_t velocity;
  dof_map<cell_type> dofMap;

  FEMatrix<cell_type> stiffness_matrix;
  FEVector<cell_type> unknown_vector;
  FEVector<cell_type> load_vector;

  GradTrialInnerGradTest<cell_type> convective_term;
  NegTrialInnerDivTest<pressure_basis_t, velocity_basis_t> pressure_term;
  DivTrialInnerTest<velocity_basis_t, pressure_basis_t> continuity_term;

  FEVector<cell_type> velocity_vector;
  FEVector<cell_type> pressure_vector;

  enum Location
  {
    TOP_EDGE,
    BOTTOM_EDGE,
    LEFT_EDGE,
    RIGHT_EDGE,
    BODY
  };

  static dof_map<cell_type> buildDofMap(const mesh<cell_type>& m,
                                        const pressure_basis_t& pressure, 
                                        const velocity_basis_t& velocity)
  {
    dof_map_builder<cell_type> mapBuilder(m);
    mapBuilder.addFiniteElement(pressure);
    mapBuilder.addFiniteElement(velocity);
    mapBuilder.handleCells(m.getCells());
    return mapBuilder.getDofMap();
  }

public:
  stokes_system(const mesh<cell_type>& _m) : m(_m), pressure(m), velocity(m), 
                                             dofMap(buildDofMap(m, pressure, velocity)),
                                             stiffness_matrix(dofMap, dofMap), 
                                             unknown_vector(dofMap), load_vector(dofMap),
                                             convective_term(&velocity, &velocity),
                                             pressure_term(&pressure, &velocity),
                                             continuity_term(&velocity, &pressure),
                                             velocity_vector(dofMap.extractDofs(&velocity)),
                                             pressure_vector(dofMap.extractDofs(&pressure))
  {  
    std::cout << "Size of dof map: " << dofMap.getMappingSize() << std::endl;
    std::cout << "Degrees of freedom: " << dofMap.getDegreesOfFreedomCount() << std::endl;
    std::cout << "Degrees of freedom on boundary: " << dofMap.getBoundaryDegreesOfFreedomCount() << std::endl;
  }

  void assemble()
  {
    const std::map<cell_id, cell_type> cells(m.getCells());

    // assemble velocity related part of momentum equation
    std::cout << "Assembling convective term..." << std::endl;
    stiffness_matrix.addTerm(m, convective_term);
    std::cout << "Assembling pressure term..." << std::endl;
    stiffness_matrix.addTerm(m, pressure_term);
    std::cout << "Assembling continuity term..." << std::endl;
    stiffness_matrix.addTerm(m, continuity_term);
    std::cout << "Assembled all terms." << std::endl;
/*
    // Iterate over cells
    for(typename std::map<cell_id, cell_type>::const_iterator cellIter(cells.begin()); cellIter != cells.end(); ++cellIter)
    {
      // Iterate over quadrature points
      const std::map<vertex_type, double> quadrature(cellIter->second.getQuadrature(m.getGeometry()));
      for(typename std::map<vertex_type, double>::const_iterator quadIter(quadrature.begin()); quadIter != quadrature.end(); ++quadIter)
      {
        // Iterate over quadratic test functions
        for(unsigned test=0; test<velocity.space_dimension(); ++test)
        {
          const typename dof_map<cell_type>::dof_t global_test = boost::make_tuple(&velocity, cellIter->first, test);
          Tensor<dimension, 0, double> test_velocity_divergence = velocity.evaluate_divergence(cellIter->second, test, quadIter->first);
          Tensor<dimension, 2, double> test_velocity_gradient = velocity.evaluate_gradient(cellIter->second, test, quadIter->first);

          // Iterate over quadratic trial functions
          for(unsigned trial=0; trial<velocity.space_dimension(); ++trial)
          {
            // assemble velocity related part of momentum equation
            const typename dof_map<cell_type>::dof_t global_trial = boost::make_tuple(&velocity, cellIter->first, trial);
            Tensor<dimension, 2, double> trial_velocity_gradient = velocity.evaluate_gradient(cellIter->second, trial, quadIter->first);

            const double convective_term = (quadIter->second * test_velocity_gradient.colon_product(trial_velocity_gradient)).toScalar();
            //stiffness_matrix.addValues(1, 1, &global_test, &global_trial, &convective_term);
          }
          
          // Iterate over linear trial functions
          for(unsigned trial=0; trial<pressure.space_dimension(); ++trial)
          {
            // assemble pressure related part of momentum equation
            const typename dof_map<cell_type>::dof_t global_trial = boost::make_tuple(&pressure, cellIter->first, trial);
            Tensor<dimension, 0, double> trial_pressure = pressure.evaluate_tensor(cellIter->second, trial, quadIter->first);

            const double pressure_term = -(quadIter->second * trial_pressure * test_velocity_divergence).toScalar();
            //stiffness_matrix.addValues(1, 1, &global_test, &global_trial, &pressure_term);
          }
        }
          
        // Iterate over linear test functions 
        for(unsigned test=0; test<pressure.space_dimension(); ++test)
        {
          const typename dof_map<cell_type>::dof_t global_test = boost::make_tuple(&pressure, cellIter->first, test);
          Tensor<dimension, 0, double> test_pressure = pressure.evaluate_tensor(cellIter->second, test, quadIter->first);

          //Iterate over quadratic trial functions
          for(unsigned trial=0; trial<velocity.space_dimension(); ++trial)
          {
            // assemble continuity equation
            const typename dof_map<cell_type>::dof_t global_trial = boost::make_tuple(&velocity, cellIter->first, trial);
            Tensor<dimension, 0, double> trial_velocity_divergence = velocity.evaluate_divergence(cellIter->second, trial, quadIter->first);

            const double continuity_term = (quadIter->second * test_pressure * trial_velocity_divergence).toScalar(); 
            //stiffness_matrix.addValues(1, 1, &global_test, &global_trial, &continuity_term);
          }
        }
      }
    }*/   

    stiffness_matrix.assemble();
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
    typedef std::map<unsigned, std::set< boost::tuple<cell_id, unsigned> > > element_dof_map;
    const element_dof_map pressure_dofs(dofMap.getBoundaryDegreesOfFreedom(&pressure));
    const element_dof_map velocity_dofs(dofMap.getBoundaryDegreesOfFreedom(&velocity));

    // Assign a value for the pressure degrees around the edges
    for(element_dof_map::const_iterator pressureIter(pressure_dofs.begin()); pressureIter != pressure_dofs.end(); ++pressureIter)
    {
      const boost::tuple<cell_id, unsigned> dofInfo(*pressureIter->second.begin());
      const vertex_type position(pressure.getDofCoordinate(boost::get<0>(dofInfo), boost::get<1>(dofInfo)));
      const Location location = getLocation(position);

      if (location == BOTTOM_EDGE || location == TOP_EDGE)
      {
        const typename dof_map<cell_type>::dof_t pressureGlobalDof = boost::make_tuple(&pressure, boost::get<0>(dofInfo), boost::get<1>(dofInfo));
        stiffness_matrix.zeroRow(pressureGlobalDof, 1.0);

        const double rhs = 0.0;
        load_vector.setValues(1, &pressureGlobalDof, &rhs);

        //To help convergence
        unknown_vector.setValues(1, &pressureGlobalDof, &rhs);
      }
    }

    // Assign values for the x velocity degrees of freedom (x^2 + y^2 around the entire boundary)
    for(element_dof_map::const_iterator velocity_iter(velocity_dofs.begin()); velocity_iter!=velocity_dofs.end(); ++velocity_iter)
    {
      assert(!velocity_iter->second.empty());  // If this failed, it would mean a degree of freedom tied to no cell
      const boost::tuple<cell_id, unsigned> dofInfo(*velocity_iter->second.begin());

      const bool isXDof = velocity.isXDof(boost::get<0>(dofInfo), boost::get<1>(dofInfo));

      if (isXDof)
      {
        const vertex_type position(velocity.getDofCoordinate(boost::get<0>(dofInfo), boost::get<1>(dofInfo)));
        // Check this really is an edge cell
        assert(getLocation(position) != BODY);

        const typename dof_map<cell_type>::dof_t velocity_globalDof = boost::make_tuple(&velocity, boost::get<0>(dofInfo), boost::get<1>(dofInfo));
        stiffness_matrix.zeroRow(velocity_globalDof, 1.0);

        const double rhs = position[0] * position[0] + position[1] * position[1]; // x^2 + y^2
        load_vector.setValues(1, &velocity_globalDof, &rhs);

        // To help convergence
        unknown_vector.setValues(1, &velocity_globalDof, &rhs);
      }
    }

    // Assign values for the y velocity degrees of freedom (zero along top and bottom edges) 
    for(element_dof_map::const_iterator velocity_iter(velocity_dofs.begin()); velocity_iter!=velocity_dofs.end(); ++velocity_iter)
    {
      assert(!velocity_iter->second.empty());  // If this failed, it would mean a degree of freedom tied to no cell
      const boost::tuple<cell_id, unsigned> dofInfo(*velocity_iter->second.begin());

      const bool isYDof = velocity.isYDof(boost::get<0>(dofInfo), boost::get<1>(dofInfo));

      if (isYDof)
      {
        const vertex_type position(velocity.getDofCoordinate(boost::get<0>(dofInfo), boost::get<1>(dofInfo)));
        const Location location = getLocation(position);

        // Check this really is an edge cell
        assert(location != BODY);

        if (location == LEFT_EDGE || location == RIGHT_EDGE)
        {
          const typename dof_map<cell_type>::dof_t velocity_globalDof = boost::make_tuple(&velocity, boost::get<0>(dofInfo), boost::get<1>(dofInfo));
          stiffness_matrix.zeroRow(velocity_globalDof, 1.0);

          const double rhs = 0.0;
          load_vector.setValues(1, &velocity_globalDof, &rhs);
        
          // To help convergence
          unknown_vector.setValues(1, &velocity_globalDof, &rhs);
        }
      }
    }

    load_vector.assemble();
    unknown_vector.assemble();
  }

  void solve()
  {
    PETScKrylovSolver solver;
    solver.solve(stiffness_matrix.getMatrixHandle(), unknown_vector.getVectorHandle(), load_vector.getVectorHandle());
  }

  void print() const
  {
    std::cout << "Stiffness Matrix: " << std::endl;
    std::cout << std::endl << "Load Vector: " << std::endl;
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

        for(unsigned dof=0; dof<velocity.space_dimension(); ++dof)
        {
          Tensor<dimension, 1, double> velocity_basis = velocity.evaluate_tensor(cellIter->second, dof, vertex);
          const typename dof_map<cell_type>::dof_t velocityDof = boost::make_tuple(&velocity, cellIter->first, dof);

          double velocityCoeff;
          unknown_vector.getValues(1u, &velocityDof, &velocityCoeff);

          xVelocity += velocity_basis(0).toScalar() * velocityCoeff;
          yVelocity += velocity_basis(1).toScalar() * velocityCoeff;
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
