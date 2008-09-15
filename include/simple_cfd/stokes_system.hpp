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

template<typename TrialType, typename TestType>
class GradTrialInnerGradTest : public FEBinaryFunction<typename TrialType::cell_type>
{
public:
  typedef typename TrialType::cell_type cell_type;
  typedef typename FEBinaryFunction<cell_type>::vertex_type vertex_type;
  typedef typename FEBinaryFunction<cell_type>::finite_element_t finite_element_t;

private:
  const TrialType* const trial;
  const TestType* const test;

public:
  GradTrialInnerGradTest(const TrialType* const _trial, const TestType* const _test) : trial(_trial), test(_test)
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

  virtual double evaluate(const std::pair<cell_id, cell_type>& cell, const std::size_t testDof, const std::size_t trialDof, const vertex_type& location) const
  {
    typename TrialType::gradient_type trial_gradient = trial->evaluate_gradient(cell.second, trialDof, location);
    typename TestType::gradient_type test_gradient = test->evaluate_gradient(cell.second, testDof, location);
    const double result = (test_gradient.colon_product(trial_gradient)).toScalar();
    return result;
  }

  ScaledFEBinaryFunction<GradTrialInnerGradTest> operator*(const double s) const
  {
    return ScaledFEBinaryFunction<GradTrialInnerGradTest>(*this, s);
  }
};

template<typename TrialType, typename TestType>
class NegTrialInnerDivTest : public FEBinaryFunction<typename TrialType::cell_type>
{
public:
  typedef typename TrialType::cell_type cell_type;
  typedef typename FEBinaryFunction<cell_type>::vertex_type vertex_type;
  typedef typename FEBinaryFunction<cell_type>::finite_element_t finite_element_t;

private:
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

  virtual double evaluate(const std::pair<cell_id, cell_type>& cell, const std::size_t testDof, const std::size_t trialDof, const vertex_type& location) const
  {
    typename TrialType::value_type trial_value = trial->evaluate_tensor(cell.second, trialDof, location);
    typename TestType::divergence_type test_divergence = test->evaluate_divergence(cell.second, testDof, location);
    const double result = - (trial_value * test_divergence).toScalar();
    return result;
  }
};

template<typename TrialType, typename TestType>
class DivTrialInnerTest : public FEBinaryFunction<typename TrialType::cell_type>
{
public:
  typedef typename TrialType::cell_type cell_type;
  typedef typename FEBinaryFunction<cell_type>::vertex_type vertex_type;
  typedef typename FEBinaryFunction<cell_type>::finite_element_t finite_element_t;

private:
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

  virtual double evaluate(const std::pair<cell_id, cell_type>& cell, const std::size_t testDof, const std::size_t trialDof, const vertex_type& location) const
  {
    typename TrialType::divergence_type trial_divergence = trial->evaluate_divergence(cell.second, trialDof, location);
    typename TestType::value_type test_value = test->evaluate_tensor(cell.second, testDof, location);
    const double result = (test_value * trial_divergence).toScalar();
    return result;
  }
};

template<typename TrialType, typename TestType>
class TrialInnerTest : public FEBinaryFunction<typename TrialType::cell_type>
{
public:
  typedef typename TrialType::cell_type cell_type;
  typedef typename FEBinaryFunction<cell_type>::vertex_type vertex_type;
  typedef typename FEBinaryFunction<cell_type>::finite_element_t finite_element_t;

private:
  const TrialType* const trial;
  const TestType* const test;

public:
  TrialInnerTest(const TrialType* const _trial, const TestType* const _test) : trial(_trial), test(_test)
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

  virtual double evaluate(const std::pair<cell_id, cell_type>& cell, const std::size_t testDof, const std::size_t trialDof, const vertex_type& location) const
  {
    typename TrialType::value_type trial_value = trial->evaluate_tensor(cell.second, trialDof, location);
    typename TestType::value_type test_value = test->evaluate_tensor(cell.second, testDof, location);
    const double result = (test_value.colon_product(trial_value)).toScalar();
    return result;
  }
};

template<typename TrialType, typename TestType>
class TrialDotGradTrialInnerTest : public FEBinaryFunction<typename TrialType::cell_type>
{
public:
  typedef typename TrialType::cell_type cell_type;
  typedef typename FEBinaryFunction<cell_type>::vertex_type vertex_type;
  typedef typename FEBinaryFunction<cell_type>::finite_element_t finite_element_t;

private:
  const TrialType* const trial;
  const FEVector<cell_type> prevTrial;
  const TestType* const test;

public:
  TrialDotGradTrialInnerTest(const TrialType* const _trial, const FEVector<cell_type>& _prevTrial, const TestType* const _test) : 
                             trial(_trial), prevTrial(_prevTrial), test(_test)
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

  virtual double evaluate(const std::pair<cell_id, cell_type>& cell, const std::size_t testDof, const std::size_t trialDof, const vertex_type& location) const
  {
    boost::tuple<const finite_element_t*, cell_id, unsigned> trialDofTuple(trial, cell.first, trialDof);
    double prevTrialCoeff;
    prevTrial.getValues(1, &trialDofTuple, &prevTrialCoeff);

    typename TrialType::value_type trial_value = trial->evaluate_tensor(cell.second, trialDof, location);
    typename TrialType::gradient_type trial_gradient = trial->evaluate_gradient(cell.second, trialDof, location) * prevTrialCoeff;
    typename TestType::value_type test_value = test->evaluate_tensor(cell.second, testDof, location);
    const double result = (trial_value.inner_product(trial_gradient)).inner_product(test_value).toScalar();
    return result;
  }

  ScaledFEBinaryFunction<TrialDotGradTrialInnerTest> operator*(const double s) const
  {
    return ScaledFEBinaryFunction<TrialDotGradTrialInnerTest>(*this, s);
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

  dof_map<cell_type> systemDofMap;
  dof_map<cell_type> velocityDofMap;
  dof_map<cell_type> pressureDofMap;

  FEMatrix<cell_type> stiffness_matrix;
  FEVector<cell_type> unknown_vector;
  FEVector<cell_type> load_vector;

  GradTrialInnerGradTest<velocity_basis_t, velocity_basis_t> convective_term;
  NegTrialInnerDivTest<pressure_basis_t, velocity_basis_t> pressure_term;
  DivTrialInnerTest<velocity_basis_t, pressure_basis_t> continuity_term;
  TrialInnerTest<velocity_basis_t, velocity_basis_t> mass_term;

  FEVector<cell_type> prev_velocity_vector;
  FEVector<cell_type> prev_pressure_vector;

  const double k;
  const double theta;

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
                                             systemDofMap(buildDofMap(m, pressure, velocity)),
                                             velocityDofMap(systemDofMap.extractDofs(&velocity)),
                                             pressureDofMap(systemDofMap.extractDofs(&pressure)),
                                             stiffness_matrix(systemDofMap, systemDofMap), 
                                             unknown_vector(systemDofMap), load_vector(systemDofMap),
                                             convective_term(&velocity, &velocity),
                                             pressure_term(&pressure, &velocity),
                                             continuity_term(&velocity, &pressure),
                                             mass_term(&velocity, &velocity),
                                             prev_velocity_vector(velocityDofMap),
                                             prev_pressure_vector(pressureDofMap),
                                             k(1.0), theta(0.5)
  {  
    std::cout << "Size of dof map: " << systemDofMap.getMappingSize() << std::endl;
    std::cout << "Degrees of freedom: " << systemDofMap.getDegreesOfFreedomCount() << std::endl;
    std::cout << "Degrees of freedom on boundary: " << systemDofMap.getBoundaryDegreesOfFreedomCount() << std::endl;
  }

  void assemble()
  {
    // assemble velocity related part of momentum equation
    stiffness_matrix.addTerm(m, convective_term);

    // assemble pressure related part of momentum equation
    stiffness_matrix.addTerm(m, pressure_term);

    // assemble continuity equation
    stiffness_matrix.addTerm(m, continuity_term);

    stiffness_matrix.assemble();
  }

  void timeDependentAssembleAndSolve() 
  {
    std::cout << "Assembling linear terms..." << std::endl;

    // Add in all constant terms in the lhs matrix
    FEMatrix<cell_type> linear_stiffness_matrix(systemDofMap, systemDofMap);
    linear_stiffness_matrix.addTerm(m, mass_term);
    linear_stiffness_matrix.addTerm(m, convective_term * (theta*k));
    linear_stiffness_matrix.addTerm(m, pressure_term);
    linear_stiffness_matrix.addTerm(m, continuity_term);
    linear_stiffness_matrix.assemble();

    // Add in all constant terms in the rhs matrix
    TrialDotGradTrialInnerTest<velocity_basis_t, velocity_basis_t> nonLinearTermPrev(&velocity, prev_velocity_vector, &velocity);
    FEMatrix<cell_type> nonlinear_rhs_matrix(velocityDofMap, velocityDofMap);
    nonlinear_rhs_matrix.addTerm(m, mass_term);
    nonlinear_rhs_matrix.addTerm(m, convective_term * (-(1.0-theta)*k));
    nonlinear_rhs_matrix.addTerm(m, nonLinearTermPrev * (-(1.0-theta)*k));
    nonlinear_rhs_matrix.assemble();

    // Add non-linear term into rhs matrix then multiply to get rhs vector
    FEVector<cell_type> rhs_velocity(nonlinear_rhs_matrix*prev_velocity_vector);
    rhs_velocity.assemble();

    // This vector will hold the guesses for the unknowns each iteration
    FEVector<cell_type> unknown_guess(systemDofMap);
    double residual = 0.0;

    do
    {
      std::cout << "Assembling non-linear terms..." << std::endl;

      // Copy the existing matrices
      FEMatrix<cell_type> nonlinear_stiffness_matrix(linear_stiffness_matrix);

      // Create non-linear terms for calculating lhs part of system
      TrialDotGradTrialInnerTest<velocity_basis_t, velocity_basis_t> nonLinearTermCurrent(&velocity, unknown_guess, &velocity);

      // Add non-linear term into stiffness matrix
      nonlinear_stiffness_matrix.addTerm(m, nonLinearTermCurrent * (theta*k));
      nonlinear_stiffness_matrix.assemble();

      // Add rhs velocity-related vector into load vector
      load_vector.zero();
      load_vector.addSubvector(rhs_velocity);
      load_vector.assemble();

      stiffness_matrix = nonlinear_stiffness_matrix;

      std::cout << "Applying boundary conditions..." << std::endl;
      applyBoundaryConditions();

      residual = ((stiffness_matrix * unknown_guess) - load_vector).two_norm();
      std::cout << "Current non-linear residual: " << residual << std::endl;

      std::cout << "Starting solver..." << std::endl;
      solve();

      unknown_guess = unknown_vector;
    }
    while(residual > 3e-4);

    unknown_vector.extractSubvector(prev_velocity_vector);
    unknown_vector.extractSubvector(prev_pressure_vector);
    prev_velocity_vector.assemble();
    prev_pressure_vector.assemble();
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
    const element_dof_map pressure_dofs(systemDofMap.getBoundaryDegreesOfFreedom(&pressure));
    const element_dof_map velocity_dofs(systemDofMap.getBoundaryDegreesOfFreedom(&velocity));

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

      const bool isXDof = velocity.getTensorIndex(boost::get<0>(dofInfo), boost::get<1>(dofInfo)) == 0;

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

      const bool isYDof = velocity.getTensorIndex(boost::get<0>(dofInfo), boost::get<1>(dofInfo)) == 1;

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

    stiffness_matrix.assemble();
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
    render(20, 20, outFile);
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
