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
#include <cstdlib>
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

  virtual double evaluate(const mesh<TriangularCell>& m, const MeshEntity& entity, const std::size_t testDof, const std::size_t trialDof, const vertex_type& location) const
  {
    const std::size_t cid = m.getContainingCell(entity);
    typename TrialType::gradient_type trial_gradient = trial->evaluate_gradient(cid, trialDof, location);
    typename TestType::gradient_type test_gradient = test->evaluate_gradient(cid, testDof, location);
    const double result = (test_gradient.colon_product(trial_gradient)).toScalar();
    return result;
  }

  ScaledFEBinaryFunction<GradTrialInnerGradTest> operator*(const double s) const
  {
    return ScaledFEBinaryFunction<GradTrialInnerGradTest>(*this, s);
  }
};

template<typename TrialType, typename TestType>
class GradTrialInnerNormalMulTest : public FEBinaryFunction<typename TrialType::cell_type>
{
public:
  typedef typename TrialType::cell_type cell_type;
  typedef typename FEBinaryFunction<cell_type>::vertex_type vertex_type;
  typedef typename FEBinaryFunction<cell_type>::finite_element_t finite_element_t;

private:
  const TrialType* const trial;
  const TestType* const test;

public:
  GradTrialInnerNormalMulTest(const TrialType* const _trial, const TestType* const _test) : trial(_trial), test(_test)
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

  virtual double evaluate(const mesh<TriangularCell>& m, const MeshEntity& entity, const std::size_t testDof, const std::size_t trialDof, const vertex_type& location) const
  {
    assert(entity.getDimension() == m.getDimension() - 1);
    const std::size_t cid = m.getContainingCell(entity);
    typename TrialType::gradient_type trial_gradient = trial->evaluate_gradient(cid, trialDof, location);
    typename TestType::value_type test_value = test->evaluate_tensor(cid, testDof, location);
    Tensor<2, 1, double> facetNormal = m.getReferenceCell().getFacetNormal(m, cid, entity.getIndex(), location);
    const double result = test_value.colon_product(trial_gradient.inner_product(facetNormal)).toScalar();
    return result;
  }

  ScaledFEBinaryFunction<GradTrialInnerNormalMulTest> operator*(const double s) const
  {
    return ScaledFEBinaryFunction<GradTrialInnerNormalMulTest>(*this, s);
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

  virtual double evaluate(const mesh<TriangularCell>& m, const MeshEntity& entity, const std::size_t testDof, const std::size_t trialDof, const vertex_type& location) const
  {
    const std::size_t cid = m.getContainingCell(entity);
    typename TrialType::value_type trial_value = trial->evaluate_tensor(cid, trialDof, location);
    typename TestType::divergence_type test_divergence = test->evaluate_divergence(cid, testDof, location);
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

  virtual double evaluate(const mesh<TriangularCell>& m, const MeshEntity& entity, const std::size_t testDof, const std::size_t trialDof, const vertex_type& location) const
  {
    const std::size_t cid = m.getContainingCell(entity);
    typename TrialType::divergence_type trial_divergence = trial->evaluate_divergence(cid, trialDof, location);
    typename TestType::value_type test_value = test->evaluate_tensor(cid, testDof, location);
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

  virtual double evaluate(const mesh<TriangularCell>& m, const MeshEntity& entity, const std::size_t testDof, const std::size_t trialDof, const vertex_type& location) const
  {
    const std::size_t cid = m.getContainingCell(entity);
    typename TrialType::value_type trial_value = trial->evaluate_tensor(cid, trialDof, location);
    typename TestType::value_type test_value = test->evaluate_tensor(cid, testDof, location);
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

  virtual double evaluate(const mesh<TriangularCell>& m, const MeshEntity& entity, const std::size_t testDof, const std::size_t trialDof, const vertex_type& location) const
  {
    const std::size_t cid = m.getContainingCell(entity);
    boost::tuple<const finite_element_t*, cell_id, unsigned> trialDofTuple(trial, entity.getIndex(), trialDof);
    double prevTrialCoeff;
    prevTrial.getValues(1, &trialDofTuple, &prevTrialCoeff);

    typename TrialType::value_type trial_value = trial->evaluate_tensor(cid, trialDof, location);
    typename TrialType::gradient_type trial_gradient = trial->evaluate_gradient(cid, trialDof, location) * prevTrialCoeff;
    typename TestType::value_type test_value = test->evaluate_tensor(cid, testDof, location);
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

  GradTrialInnerGradTest<velocity_basis_t, velocity_basis_t> viscosity_term;
  GradTrialInnerNormalMulTest<velocity_basis_t, velocity_basis_t> viscosity_boundary_term;
  NegTrialInnerDivTest<pressure_basis_t, velocity_basis_t> pressure_term;
  DivTrialInnerTest<velocity_basis_t, pressure_basis_t> continuity_term;
  TrialInnerTest<velocity_basis_t, velocity_basis_t> mass_term;

  FEVector<cell_type> prev_velocity_vector;
  FEVector<cell_type> prev_pressure_vector;

  const double k;
  const double theta;
  const double kinematic_viscosity;

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
    return mapBuilder.getDofMap();
  }

public:
  stokes_system(const mesh<cell_type>& _m) : m(_m), pressure(m), velocity(m), 
                                             systemDofMap(buildDofMap(m, pressure, velocity)),
                                             velocityDofMap(systemDofMap.extractDofs(&velocity)),
                                             pressureDofMap(systemDofMap.extractDofs(&pressure)),
                                             stiffness_matrix(systemDofMap, systemDofMap), 
                                             unknown_vector(systemDofMap), load_vector(systemDofMap),
                                             viscosity_term(&velocity, &velocity),
                                             viscosity_boundary_term(&velocity, &velocity),
                                             pressure_term(&pressure, &velocity),
                                             continuity_term(&velocity, &pressure),
                                             mass_term(&velocity, &velocity),
                                             prev_velocity_vector(velocityDofMap),
                                             prev_pressure_vector(pressureDofMap),
                                             k(1.0), theta(0.5), kinematic_viscosity(1.0/250)
  {  
    std::cout << "Size of dof map: " << systemDofMap.getMappingSize() << std::endl;
    std::cout << "Degrees of freedom: " << systemDofMap.getDegreesOfFreedomCount() << std::endl;
    std::cout << "Degrees of freedom on boundary: " << systemDofMap.getBoundaryDegreesOfFreedomCount() << std::endl;
  }

  void assemble()
  {
    // assemble velocity related part of momentum equation
    stiffness_matrix.addTerm(m, viscosity_term);

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
    linear_stiffness_matrix.addTerm(m, viscosity_term * (theta * k * kinematic_viscosity));
    linear_stiffness_matrix.addBoundaryTerm(m, viscosity_boundary_term * (theta * k * kinematic_viscosity * -1.0));
    linear_stiffness_matrix.addTerm(m, pressure_term);
    linear_stiffness_matrix.addTerm(m, continuity_term);
    linear_stiffness_matrix.assemble();

    // Add in all constant terms in the rhs matrix
    TrialDotGradTrialInnerTest<velocity_basis_t, velocity_basis_t> nonLinearTermPrev(&velocity, prev_velocity_vector, &velocity);
    FEMatrix<cell_type> nonlinear_rhs_matrix(velocityDofMap, velocityDofMap);
    nonlinear_rhs_matrix.addTerm(m, mass_term);
    nonlinear_rhs_matrix.addTerm(m, viscosity_term * -((1.0-theta) * k * kinematic_viscosity));
    nonlinear_rhs_matrix.addBoundaryTerm(m, viscosity_term * ((1.0-theta) * k * kinematic_viscosity));
    nonlinear_rhs_matrix.addTerm(m, nonLinearTermPrev * (-(1.0-theta)*k));
    nonlinear_rhs_matrix.assemble();

    // Add non-linear term into rhs matrix then multiply to get rhs vector
    FEVector<cell_type> rhs_velocity(nonlinear_rhs_matrix*prev_velocity_vector);
    rhs_velocity.assemble();

    // This vector will hold the guesses for the unknowns each iteration
    FEVector<cell_type> unknown_guess(unknown_vector);
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
      applyCylinderBoundaryConditions();

      std::cout << "Starting solver..." << std::endl;
      solve();

      unknown_guess = unknown_vector;
      residual = ((stiffness_matrix * unknown_guess) - load_vector).two_norm();
      std::cout << "Current non-linear residual: " << residual << std::endl;
    }
    while(residual > 1e-3);

    unknown_vector.extractSubvector(prev_velocity_vector);
    unknown_vector.extractSubvector(prev_pressure_vector);
    prev_velocity_vector.assemble();
    prev_pressure_vector.assemble();
  }

  Location getLocation(const vertex_type& v)
  {
    // FIXME: make me into a function that doesn't depend on the specifics of the mesh generation

    Location location = BODY;

    if (v[0] == 3.0)
      location = RIGHT_EDGE;

    if (v[1] == 0.0)
      location = BOTTOM_EDGE;

    if (v[1] == 1.0)
      location = TOP_EDGE;

    if (v[0] == 0.0)
      location = LEFT_EDGE;

    return location;
  }

  void applyBoundaryConditions()
  {
    typedef std::map<unsigned, std::set< boost::tuple<cell_id, unsigned> > > element_dof_map;
    const element_dof_map velocity_dofs(systemDofMap.getBoundaryDegreesOfFreedom(&velocity));

    // Assign values for the x velocity degrees of freedom on the left hand side
    for(element_dof_map::const_iterator velocity_iter(velocity_dofs.begin()); velocity_iter!=velocity_dofs.end(); ++velocity_iter)
    {
      assert(!velocity_iter->second.empty());  // If this failed, it would mean a degree of freedom tied to no cell
      const boost::tuple<cell_id, unsigned> dofInfo(*velocity_iter->second.begin());

      const bool isXDof = velocity.getTensorIndex(boost::get<0>(dofInfo), boost::get<1>(dofInfo)) == 0;
      const vertex_type position(velocity.getDofCoordinate(boost::get<0>(dofInfo), boost::get<1>(dofInfo)));
      
      if (getLocation(position) == LEFT_EDGE)
      {
        const typename dof_map<cell_type>::dof_t velocity_globalDof = boost::make_tuple(&velocity, boost::get<0>(dofInfo), boost::get<1>(dofInfo));
        stiffness_matrix.zeroRow(velocity_globalDof, 1.0);

        // Set x velocity to same value on inflow and outflow boundary, and y velocity to 0
        const double rhs = isXDof ? 5.0 : 0.0;
        load_vector.setValues(1, &velocity_globalDof, &rhs);

        // To help convergence
        unknown_vector.setValues(1, &velocity_globalDof, &rhs);
      }
    }

    // Assign values for the velocity degrees of freedom along top and bottom edges
    for(element_dof_map::const_iterator velocity_iter(velocity_dofs.begin()); velocity_iter!=velocity_dofs.end(); ++velocity_iter)
    {
      assert(!velocity_iter->second.empty());  // If this failed, it would mean a degree of freedom tied to no cell
      const boost::tuple<cell_id, unsigned> dofInfo(*velocity_iter->second.begin());

      const bool isXDof = velocity.getTensorIndex(boost::get<0>(dofInfo), boost::get<1>(dofInfo)) == 0;
      const vertex_type position(velocity.getDofCoordinate(boost::get<0>(dofInfo), boost::get<1>(dofInfo)));
      const Location location = getLocation(position);

      // Check this really is an edge cell
      assert(location != BODY);

      if ((location == TOP_EDGE || location == BOTTOM_EDGE) && !isXDof)
      {
        const typename dof_map<cell_type>::dof_t velocity_globalDof = boost::make_tuple(&velocity, boost::get<0>(dofInfo), boost::get<1>(dofInfo));
        stiffness_matrix.zeroRow(velocity_globalDof, 1.0);

        const double rhs = 0.0;
        load_vector.setValues(1, &velocity_globalDof, &rhs);
        
        // To help convergence
        unknown_vector.setValues(1, &velocity_globalDof, &rhs);
      }
    }

    stiffness_matrix.assemble();
    load_vector.assemble();
    unknown_vector.assemble();
  }

  void applyCylinderBoundaryConditions()
  {
    const unsigned velocitySpaceDimension = velocity.space_dimension();
    const vertex_type centre(1.0, 0.5);
    const double radius = 0.15;

    for(typename mesh<cell_type>::global_iterator cellIter(m.global_begin(dimension)); cellIter!=m.global_end(dimension); ++cellIter)
    {
      for(unsigned dof=0; dof<velocitySpaceDimension; ++dof)
      {
        const vertex_type dofLocation = velocity.getDofCoordinate(cellIter->getIndex(), dof);
        const vertex_type offset = dofLocation - centre;

        if((offset[0] * offset[0] + offset[1] * offset[1]) < radius * radius)
        {
          const typename dof_map<cell_type>::dof_t velocity_globalDof = boost::make_tuple(&velocity, cellIter->getIndex(), dof);
          stiffness_matrix.zeroRow(velocity_globalDof, 1.0);
          const double rhs = 0.0;
          load_vector.setValues(1, &velocity_globalDof, &rhs);
        }
      }
    }
  }

  void solve()
  {
    PETScKrylovSolver solver;
    solver.setMaxIterations(25000);
    solver.solve(stiffness_matrix.getMatrixHandle(), unknown_vector.getVectorHandle(), load_vector.getVectorHandle());

    if (!solver.converged())
    {
      std::cout << "Convergence failure: " << solver.getConvergedReason() << std::endl;
      exit(EXIT_FAILURE);
    }
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
    //renderOld(3.0, 1.0, 90, 30, outFile);
    render(outFile);
    outFile.close();
  }

  std::pair<double, double> getVelocityVector(const vertex_type& vertex) const
  {
    for(typename mesh<cell_type>::global_iterator cellIter(m.global_begin(dimension)); cellIter!=m.global_end(dimension); ++cellIter)
    {
      if (m.getReferenceCell().contains(m, cellIter->getIndex(), vertex))
      {
        return getVelocityVector(MeshEntity(m.getDimension(), cellIter->getIndex()), vertex);
      }
    }
    return std::make_pair(0.0, 0.0);
  }

  std::pair<double, double> getVelocityVector(const MeshEntity& entity, const vertex_type& vertex) const
  {
    // From any given mesh entity, we determine a cell we can compute a value from
    std::size_t cid;
    if (entity.getDimension() == m.getDimension())
    {
      cid = entity.getIndex();
    }
    else
    {
      const std::vector<std::size_t> cellIndices(m.getIndices(entity, m.getDimension()));
      assert(cellIndices.size() > 0);
      cid = cellIndices.front();
    }

    //assert(m.getReferenceCell().contains(m, cid, vertex));
    double xVelocity(0.0), yVelocity(0.0);

    for(unsigned dof=0; dof<velocity.space_dimension(); ++dof)
    {
      Tensor<dimension, 1, double> velocity_basis = velocity.evaluate_tensor(cid, dof, vertex);
      const typename dof_map<cell_type>::dof_t velocityDof = boost::make_tuple(&velocity, cid, dof);

      double velocityCoeff;
      unknown_vector.getValues(1u, &velocityDof, &velocityCoeff);

      xVelocity += velocity_basis(0).toScalar() * velocityCoeff;
      yVelocity += velocity_basis(1).toScalar() * velocityCoeff;
     }
     return std::make_pair(xVelocity, yVelocity);
  }

  void renderOld(const double width, const double height, const unsigned xPoints, const unsigned yPoints, std::ostream& out)
  {
    const double xSpacing = width/(xPoints-1);
    const double ySpacing = height/(yPoints-1);

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
        const double x = xPoint * xSpacing;
        const double y = yPoint * ySpacing;

        const std::pair<double, double> velocity(getVelocityVector(vertex_type(x, y)));
        out << velocity.first << " " << velocity.second << " " << 0.0 << std::endl;
      }
    }
  }

  void render(std::ostream& out)
  {
    std::vector< vertex<2> > vertices;
    std::vector< std::pair<double, double> > velocities;

    for(typename mesh<cell_type>::global_iterator vIter(m.global_begin(0)); vIter!=m.global_end(0); ++vIter)
    {
      const vertex<2> v(m.getVertex(vIter->getIndex()));
      vertices.push_back(v);
      velocities.push_back(getVelocityVector(v));
    }

    for(typename mesh<cell_type>::global_iterator eIter(m.global_begin(1)); eIter!=m.global_end(1); ++eIter)
    {
      const std::vector<std::size_t> vertexIndices(m.getIndices(*eIter, 0));
      const vertex<2> v((m.getVertex(vertexIndices[0]) + m.getVertex(vertexIndices[1]))/2);
      vertices.push_back(v);
      velocities.push_back(getVelocityVector(*eIter, v));
    }

    out << "# vtk DataFile Version 2.0" << std::endl;
    out << "Simple Navier-Stokes Solver" << std::endl;
    out << "ASCII" << std::endl;
    out << "DATASET POLYDATA" << std::endl;
    out << "POINTS " << vertices.size() << " DOUBLE " << std::endl;

    for(std::size_t point = 0; point < vertices.size(); ++point)
    {
      const vertex<2> v(vertices[point]);
      out << v[0] << " " << v[1] << " 0" << std::endl;
    }

    out << "POLYGONS " << m.numEntities(dimension) << " " << m.numRelations(dimension, 0)+m.numEntities(dimension)  << std::endl; 
    for(typename mesh<cell_type>::global_iterator cIter(m.global_begin(dimension)); cIter!=m.global_end(dimension); ++cIter)
    {
      const std::vector<std::size_t> vIndices(m.getIndices(*cIter, 0));
      out << vIndices.size();

      for(std::size_t v=0; v<vIndices.size(); ++v)
      {
        out << " " << vIndices[v];
      }
      out << std::endl;
    }

    out << "POINT_DATA " << vertices.size() << std::endl;
    out << "VECTORS velocity_field DOUBLE" << std::endl;

    for(std::size_t point = 0; point < velocities.size(); ++point)
    {
      out << velocities[point].first << " " << velocities[point].second << " 0" << std::endl;
    }
  }
};

}

#endif
