#ifndef SIMPLE_CFD_STOKES_SYSTEM
#define SIMPLE_CFD_STOKES_SYSTEM

#include <iostream>
#include <map>
#include <algorithm>
#include <ostream>
#include <fstream>
#include <utility>
#include <cstdlib>
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
#include "subdomain.hpp"
#include "function.hpp"
#include "boundary_condition.hpp"
#include "boundary_condition_handler.hpp"
#include "cell_vertices.hpp"
#include "dof.hpp"

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

  virtual double evaluate(const CellVertices<cell_type::dimension>& vertices, const MeshEntity& gEntity, 
    const MeshEntity& lEntity, const std::size_t testDof, const std::size_t trialDof, const vertex_type& location) const
  {
    typename TrialType::gradient_type trial_gradient = trial->evaluate_gradient(vertices, trialDof, location);
    typename TestType::gradient_type test_gradient = test->evaluate_gradient(vertices, testDof, location);
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

  virtual double evaluate(const CellVertices<cell_type::dimension>& vertices, const MeshEntity& gEntity, 
    const MeshEntity& lEntity, const std::size_t testDof, const std::size_t trialDof, const vertex_type& location) const
  {
    assert(gEntity.getDimension() == cell_type::dimension - 1);
    typename TrialType::gradient_type trial_gradient = trial->evaluate_gradient(vertices, trialDof, location);
    typename TestType::value_type test_value = test->evaluate_tensor(vertices, testDof, location);
    Tensor<2, 1> facetNormal = cell_type().getFacetNormal(vertices, lEntity.getIndex(), location);
    const double result = test_value.colon_product(trial_gradient.inner_product(facetNormal)).toScalar();
    return result;
  }

  ScaledFEBinaryFunction<GradTrialInnerNormalMulTest> operator*(const double s) const
  {
    return ScaledFEBinaryFunction<GradTrialInnerNormalMulTest>(*this, s);
  }
};


template<typename TrialType, typename TestType>
class TrialInnerDivTest : public FEBinaryFunction<typename TrialType::cell_type>
{
public:
  typedef typename TrialType::cell_type cell_type;
  typedef typename FEBinaryFunction<cell_type>::vertex_type vertex_type;
  typedef typename FEBinaryFunction<cell_type>::finite_element_t finite_element_t;

private:
  const TrialType* const trial;
  const TestType* const test;

public:
  TrialInnerDivTest(const TrialType* const _trial, const TestType* const _test) : trial(_trial), test(_test)
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

  virtual double evaluate(const CellVertices<cell_type::dimension>& vertices, const MeshEntity& gEntity, 
    const MeshEntity& lEntity, const std::size_t testDof, const std::size_t trialDof, const vertex_type& location) const
  {
    typename TrialType::value_type trial_value = trial->evaluate_tensor(vertices, trialDof, location);
    typename TestType::divergence_type test_divergence = test->evaluate_divergence(vertices, testDof, location);
    const double result = (trial_value * test_divergence).toScalar();
    return result;
  }

  ScaledFEBinaryFunction<TrialInnerDivTest> operator*(const double s) const
  {
    return ScaledFEBinaryFunction<TrialInnerDivTest>(*this, s);
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

  virtual double evaluate(const CellVertices<cell_type::dimension>& vertices, const MeshEntity& gEntity, 
    const MeshEntity& lEntity, const std::size_t testDof, const std::size_t trialDof, const vertex_type& location) const
  {
    typename TrialType::divergence_type trial_divergence = trial->evaluate_divergence(vertices, trialDof, location);
    typename TestType::value_type test_value = test->evaluate_tensor(vertices, testDof, location);
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

  virtual double evaluate(const CellVertices<cell_type::dimension>& vertices, const MeshEntity& gEntity, 
    const MeshEntity& lEntity, const std::size_t testDof, const std::size_t trialDof, const vertex_type& location) const
  {
    typename TrialType::value_type trial_value = trial->evaluate_tensor(vertices, trialDof, location);
    typename TestType::value_type test_value = test->evaluate_tensor(vertices, testDof, location);
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

  virtual double evaluate(const CellVertices<cell_type::dimension>& vertices, const MeshEntity& gEntity, 
    const MeshEntity& lEntity, const std::size_t testDofIndex, const std::size_t trialDofIndex, const vertex_type& location) const
  {
    typedef Dof<cell_type::dimension> dof_t;
    const dof_t trialDof(trial, gEntity.getIndex(), trialDofIndex);
    double prevTrialCoeff;
    prevTrial.getValues(1, &trialDof, &prevTrialCoeff);

    typename TrialType::value_type trial_value = trial->evaluate_tensor(vertices, trialDofIndex, location);
    typename TrialType::gradient_type trial_gradient = trial->evaluate_gradient(vertices, trialDofIndex, location) * prevTrialCoeff;
    typename TestType::value_type test_value = test->evaluate_tensor(vertices, testDofIndex, location);
    const double result = (trial_value.inner_product(trial_gradient)).inner_product(test_value).toScalar();
    return result;
  }

  ScaledFEBinaryFunction<TrialDotGradTrialInnerTest> operator*(const double s) const
  {
    return ScaledFEBinaryFunction<TrialDotGradTrialInnerTest>(*this, s);
  }
};

class Edges : public SubDomain<2>
{
private:
  const double EPSILON;

public:
  Edges() : EPSILON(1e-5)
  {
  }

  bool inside(const vertex<dimension>& v) const
  {
    return v[0] < EPSILON || v[1] < EPSILON || v[1] > (1.0 - EPSILON);
  }
};

class Cylinder : public SubDomain<2>
{
public:
  bool inside(const vertex<dimension>& v) const
  {
    const vertex<dimension> centre(1.0, 0.5);
    const double radius = 0.15;
    const vertex<dimension> offset(centre - v);

    return (offset[0]*offset[0] + offset[1]*offset[1]) < (radius * radius);
  }
};

class Zero : public Function<2, 1>
{
public:
  virtual Tensor<2, 1> evaluate(const vertex<dimension>& v) const
  {
    return Tensor<2, 1>();
  }
};

class EdgeConditions : public Function<2, 1>
{
private:
  const double EPSILON;

public:
  EdgeConditions() : EPSILON(1e-5)
  {
  }

  virtual Tensor<2, 1> evaluate(const vertex<dimension>& v) const
  {
    Tensor<2, 1> t;

    if (v[0] < EPSILON && v[1] > EPSILON && v[1] < 1.0 - EPSILON)
    {
      t(0) = 5.0;
    }

    return t;
  }
};

template<typename C>
class StokesSystem
{
public:
  typedef C cell_type;

private:
  static const std::size_t dimension = cell_type::dimension;
  typedef typename DofMap<cell_type>::dof_t dof_t;
  typedef vertex<dimension> vertex_type;
  typedef FiniteElement<dimension> finite_element_t;
  typedef LagrangeTriangleLinear<0> pressure_basis_t;
  typedef LagrangeTriangleQuadratic<1> velocity_basis_t;

  const Mesh<dimension>& m;
  pressure_basis_t pressure;
  velocity_basis_t velocity;

  DofMap<cell_type> systemDofMap;
  DofMap<cell_type> velocityDofMap;
  DofMap<cell_type> pressureDofMap;

  DofMap<cell_type> velocityDofMapHomogeneous;
  DofMap<cell_type> velocityDofMapDirichlet;

  GradTrialInnerGradTest<velocity_basis_t, velocity_basis_t> viscosity_term;
  GradTrialInnerNormalMulTest<velocity_basis_t, velocity_basis_t> viscosity_boundary_term;
  TrialInnerDivTest<pressure_basis_t, velocity_basis_t> pressure_term;
  DivTrialInnerTest<velocity_basis_t, pressure_basis_t> continuity_term;
  TrialInnerTest<velocity_basis_t, velocity_basis_t> mass_term;

  FEVector<cell_type> prev_velocity_vector;
  FEVector<cell_type> prev_pressure_vector;
  FEVector<cell_type> velocity_vector;
  FEVector<cell_type> pressure_vector;

  const double k;
  const double theta;
  const double kinematic_viscosity;

  // Subdomains
  Edges edges;
  Cylinder cylinder;

  // Functions
  EdgeConditions edgeConditions;
  Zero zero;

  // Boundary Conditions
  BoundaryCondition<2,1> edgeVelocities;
  BoundaryCondition<2,1> cylinderVelocities;

  enum Location
  {
    TOP_EDGE,
    BOTTOM_EDGE,
    LEFT_EDGE,
    RIGHT_EDGE,
    BODY
  };

  static DofMap<cell_type> buildDofMap(const Mesh<dimension>& m,
                                        const pressure_basis_t& pressure, 
                                        const velocity_basis_t& velocity)
  {
    DofMapBuilder<cell_type> mapBuilder(m);
    mapBuilder.addFiniteElement(pressure);
    mapBuilder.addFiniteElement(velocity);
    return mapBuilder.getDofMap();
  }

public:
  StokesSystem(const Mesh<dimension>& _m) : m(_m), systemDofMap(buildDofMap(m, pressure, velocity)),
                                             velocityDofMap(systemDofMap.extractDofs(&velocity)),
                                             pressureDofMap(systemDofMap.extractDofs(&pressure)),
                                             viscosity_term(&velocity, &velocity),
                                             viscosity_boundary_term(&velocity, &velocity),
                                             pressure_term(&pressure, &velocity),
                                             continuity_term(&velocity, &pressure),
                                             mass_term(&velocity, &velocity),
                                             prev_velocity_vector(velocityDofMap),
                                             prev_pressure_vector(pressureDofMap),
                                             velocity_vector(velocityDofMap),
                                             pressure_vector(pressureDofMap),
                                             k(0.01), theta(0.5), kinematic_viscosity(1.0/250),
                                             edgeVelocities(edges, edgeConditions), 
                                             cylinderVelocities(cylinder, zero)
  {  
    std::cout << "Size of dof map: " << systemDofMap.getMappingSize() << std::endl;
    std::cout << "Degrees of freedom: " << systemDofMap.getDegreesOfFreedomCount() << std::endl;

    std::vector< std::pair<const finite_element_t*, const SubDomain<dimension>*> > boundaryConditions;
    boundaryConditions.push_back(std::make_pair(&velocity, &edges));
    boundaryConditions.push_back(std::make_pair(&velocity, &cylinder));

    const std::pair< DofMap<cell_type>, DofMap<cell_type> > splitDofs(velocityDofMap.splitHomogeneousDirichlet(boundaryConditions));
    velocityDofMapHomogeneous = splitDofs.first;
    velocityDofMapDirichlet = splitDofs.second;
  }

  FEMatrix<cell_type> getLumpedInverse(const FEMatrix<cell_type>& matrix)
  {
    FEVector<cell_type> diagonal(matrix.getLumpedDiagonal());
    diagonal.reciprocal();
    FEMatrix<cell_type> invertedMatrix(matrix.getRowMappings(), matrix.getColMappings());
    invertedMatrix.addToDiagonal(diagonal);
    return invertedMatrix;
  }

  void coupledSolve()
  {
    prev_velocity_vector = velocity_vector;
    prev_pressure_vector = pressure_vector;

    std::cout << "Assembling linear terms..." << std::endl;

    // Add in all constant terms in the lhs matrix
    FEMatrix<cell_type> linear_stiffness_matrix(systemDofMap, systemDofMap);
    linear_stiffness_matrix.addTerm(m, mass_term);
    linear_stiffness_matrix.addTerm(m, viscosity_term * (theta * k * kinematic_viscosity));
    linear_stiffness_matrix.addBoundaryTerm(m, viscosity_boundary_term * (theta * k * kinematic_viscosity * -1.0));
    linear_stiffness_matrix.addTerm(m, pressure_term * -1.0);
    linear_stiffness_matrix.addTerm(m, continuity_term);
    linear_stiffness_matrix.assemble();

    // Add in all constant terms in the rhs matrix
    TrialDotGradTrialInnerTest<velocity_basis_t, velocity_basis_t> nonLinearTermPrev(&velocity, prev_velocity_vector, &velocity);
    FEMatrix<cell_type> nonlinear_rhs_matrix(velocityDofMap, velocityDofMap);
    nonlinear_rhs_matrix.addTerm(m, mass_term);
    nonlinear_rhs_matrix.addTerm(m, viscosity_term * -((1.0-theta) * k * kinematic_viscosity));
    nonlinear_rhs_matrix.addBoundaryTerm(m, viscosity_boundary_term * ((1.0-theta) * k * kinematic_viscosity));
    nonlinear_rhs_matrix.addTerm(m, nonLinearTermPrev * (-(1.0-theta)*k));
    nonlinear_rhs_matrix.assemble();

    // Add non-linear term into rhs matrix then multiply to get rhs vector
    FEVector<cell_type> rhs_velocity(nonlinear_rhs_matrix*prev_velocity_vector);
    rhs_velocity.assemble();

    // This vector will hold the guesses for the unknowns each iteration
    FEVector<cell_type> load_vector(systemDofMap);
    FEVector<cell_type> unknown_vector(systemDofMap);
    FEVector<cell_type> unknown_guess(unknown_vector);
    double residual = 0.0;

    while(true)
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

      std::cout << "Applying boundary conditions..." << std::endl;
      applyEdgeVelocityBoundaryConditions(nonlinear_stiffness_matrix, unknown_vector, load_vector);
      applyCylinderVelocityBoundaryConditions(nonlinear_stiffness_matrix, unknown_vector, load_vector);

      // Calculate residual here
      residual = ((nonlinear_stiffness_matrix * unknown_guess) - load_vector).two_norm();
      std::cout << "Current non-linear residual: " << residual << std::endl;

      if (residual <= 1e-3)
        break;

      std::cout << "Starting solver..." << std::endl;
      solve(nonlinear_stiffness_matrix, unknown_vector, load_vector, false);

      unknown_guess = unknown_vector;
    }

    unknown_vector.extractSubvector(pressure_vector);
    unknown_vector.extractSubvector(velocity_vector);
    pressure_vector.assemble();
    velocity_vector.assemble();
  }

  void projectionSolve() 
  {
    prev_velocity_vector = velocity_vector;
    prev_pressure_vector = pressure_vector;

    std::cout << "Assembling linear terms..." << std::endl;

    FEVector<cell_type> homogeneous_prev_velocity_vector(velocityDofMapHomogeneous);
    prev_velocity_vector.extractSubvector(homogeneous_prev_velocity_vector);
    homogeneous_prev_velocity_vector.assemble();

    // Add in all constant terms in the lhs matrix
    FEMatrix<cell_type> linear_lhs_matrix(velocityDofMapHomogeneous, velocityDofMapHomogeneous);
    FEMatrix<cell_type> linear_dirichlet_rhs_matrix(velocityDofMapHomogeneous, velocityDofMapDirichlet);

    linear_lhs_matrix.addTerm(m, mass_term);
    linear_lhs_matrix.addTerm(m, viscosity_term * (theta * k * kinematic_viscosity));
    linear_lhs_matrix.assemble();

    linear_dirichlet_rhs_matrix.addTerm(m, mass_term);
    linear_dirichlet_rhs_matrix.addTerm(m, viscosity_term * (theta * k * kinematic_viscosity));
    linear_dirichlet_rhs_matrix.assemble();

    // Add in all constant terms in the rhs matrix
    TrialDotGradTrialInnerTest<velocity_basis_t, velocity_basis_t> nonLinearTermPrev(&velocity, prev_velocity_vector, &velocity);
    FEMatrix<cell_type> nonlinear_rhs_matrix(velocityDofMapHomogeneous, velocityDofMapHomogeneous);
    nonlinear_rhs_matrix.addTerm(m, mass_term);
    nonlinear_rhs_matrix.addTerm(m, viscosity_term * (-(1.0-theta) * k * kinematic_viscosity));
    nonlinear_rhs_matrix.addTerm(m, nonLinearTermPrev * (-(1.0-theta) * k));
    nonlinear_rhs_matrix.assemble();

    FEMatrix<cell_type> pressure_matrix(velocityDofMapHomogeneous, pressureDofMap);
    pressure_matrix.addTerm(m, pressure_term);
    pressure_matrix.assemble();

    FEMatrix<cell_type> velocity_mass_matrix(velocityDofMapHomogeneous, velocityDofMapHomogeneous);
    velocity_mass_matrix.addTerm(m, mass_term);
    velocity_mass_matrix.assemble();

    FEMatrix<cell_type> inverted_mass_matrix(getLumpedInverse(velocity_mass_matrix));

    // Add non-linear term into rhs matrix then multiply to get rhs vector
    FEVector<cell_type> rhs_velocity(nonlinear_rhs_matrix*homogeneous_prev_velocity_vector);

    FEVector<cell_type> dirichletValues(velocityDofMapDirichlet);
    edgeVelocities.populateDirichletValues(dirichletValues, velocity);
    cylinderVelocities.populateDirichletValues(dirichletValues, velocity);
    dirichletValues.assemble();

    // This vector will hold the guesses for the unknowns each iteration
    FEVector<cell_type> velocity_guess(prev_velocity_vector);
    FEVector<cell_type> unknown_velocity(velocityDofMapHomogeneous);
    FEVector<cell_type> pressure_guess(prev_pressure_vector);
    double residual = 0.0;

    for(int picard_iteration=0; picard_iteration<2; ++picard_iteration)
    {
      std::cout << "Assembling non-linear terms..." << std::endl;

      pressure_guess = prev_pressure_vector;

      // Copy the existing matrices
      FEMatrix<cell_type> nonlinear_lhs_matrix(linear_lhs_matrix);
      FEMatrix<cell_type> nonlinear_dirichlet_rhs_matrix(linear_dirichlet_rhs_matrix);

      // Create non-linear terms for calculating lhs part of system
      TrialDotGradTrialInnerTest<velocity_basis_t, velocity_basis_t> nonLinearTermCurrent(&velocity, velocity_guess, &velocity);

      // Add non-linear term into stiffness matrix
      nonlinear_lhs_matrix.addTerm(m, nonLinearTermCurrent * (theta*k));
      nonlinear_lhs_matrix.assemble();

      nonlinear_dirichlet_rhs_matrix.addTerm(m, nonLinearTermCurrent * (theta*k));
      nonlinear_dirichlet_rhs_matrix.assemble();

      //std::cout << "Applying boundary conditions..." << std::endl;
      //applyEdgeVelocityBoundaryConditions(nonlinear_lhs_matrix, unknown_velocity, rhs_velocity);
      //applyCylinderVelocityBoundaryConditions(nonlinear_lhs_matrix, unknown_velocity, rhs_velocity);

      std::cout << "Starting solver..." << std::endl;

      FEVector<cell_type> modified_rhs_velocity(rhs_velocity - nonlinear_dirichlet_rhs_matrix*dirichletValues + pressure_matrix * pressure_guess * k);
      solve(nonlinear_lhs_matrix, unknown_velocity, modified_rhs_velocity);

      std::cout << "L2-norm of homogeneous velocity vector: " << unknown_velocity.two_norm() << std::endl;
      std::cout << "L2-norm of C^Tu: " << pressure_matrix.trans_mult(unknown_velocity).two_norm() << std::endl;

      residual = ((nonlinear_lhs_matrix * unknown_velocity) - modified_rhs_velocity).two_norm();
      std::cout << "Current non-linear residual in momentum equation: " << residual << std::endl;

      // Now solve mass-lumped continuity equation
      FEMatrix<cell_type> continuity_lhs(pressure_matrix.trans_mult(inverted_mass_matrix)*pressure_matrix);
      FEVector<cell_type> continuity_rhs(pressure_matrix.trans_mult(inverted_mass_matrix*(velocity_mass_matrix*unknown_velocity*-1.0)));
      FEVector<cell_type> phi(pressureDofMap);
      std::cout << "Solving for phi..." << std::endl;
      solve(continuity_lhs, phi, continuity_rhs);
      std::cout << "L2-norm of phi: " << phi.two_norm() << std::endl;

      //pressure_guess = pressure_guess - (phi * (1/k));

      FEVector<cell_type> velocity_correction_rhs(velocity_mass_matrix*unknown_velocity + pressure_matrix*phi);
      std::cout << "Solving velocity correction..." << std::endl;
      solve(velocity_mass_matrix, unknown_velocity, velocity_correction_rhs); 
      std::cout << "L2-norm of corrected homogeneous velocity vector: " << unknown_velocity.two_norm() << std::endl;
      std::cout << "L2-norm of C^Tu with modified velocity vector: " << pressure_matrix.trans_mult(unknown_velocity).two_norm() << std::endl;

      velocity_guess.zero();
      velocity_guess.addSubvector(unknown_velocity);
      velocity_guess.addSubvector(dirichletValues);
      velocity_guess.assemble();
    }

    velocity_vector.zero();
    velocity_vector.addSubvector(unknown_velocity);
    velocity_vector.addSubvector(dirichletValues);
    velocity_vector.assemble();

    pressure_vector = pressure_guess;
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

  void applyEdgeVelocityBoundaryConditions(FEMatrix<cell_type>& stiffness_matrix, FEVector<cell_type>& unknown_vector, FEVector<cell_type>& load_vector)
  {
    const unsigned velocitySpaceDimension = velocity.spaceDimension();

    for(typename Mesh<dimension>::global_iterator cellIter(m.global_begin(dimension)); cellIter!=m.global_end(dimension); ++cellIter)
    {
      for(unsigned dof=0; dof<velocitySpaceDimension; ++dof)
      {
        const dof_t velocity_globalDof = dof_t(&velocity, cellIter->getIndex(), dof);
        const vertex_type dofLocation = velocity.getDofCoordinateGlobal(m, cellIter->getIndex(), dof);
        const bool isXDof = velocity.getTensorIndex(m, cellIter->getIndex(), dof) == 0;
        const Location location = getLocation(dofLocation);

        if (location == LEFT_EDGE)
        {
          stiffness_matrix.zeroRow(velocity_globalDof, 1.0);

          // Set x velocity to same value on inflow and outflow boundary, and y velocity to 0
          const double rhs = isXDof ? 5.0 : 0.0;
          load_vector.setValues(1, &velocity_globalDof, &rhs);

          // To help convergence
          unknown_vector.setValues(1, &velocity_globalDof, &rhs);
        }
        else if ((location == TOP_EDGE || location == BOTTOM_EDGE) && !isXDof)
        {
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

  void applyCylinderVelocityBoundaryConditions(FEMatrix<cell_type>& stiffness_matrix, FEVector<cell_type>& unknown_vector, FEVector<cell_type>& load_vector)
  {
    const unsigned velocitySpaceDimension = velocity.spaceDimension();
    const vertex_type centre(1.0, 0.5);
    const double radius = 0.15;

    for(typename Mesh<dimension>::global_iterator cellIter(m.global_begin(dimension)); cellIter!=m.global_end(dimension); ++cellIter)
    {
      for(unsigned dof=0; dof<velocitySpaceDimension; ++dof)
      {
        const vertex_type dofLocation = velocity.getDofCoordinateGlobal(m, cellIter->getIndex(), dof);
        const vertex_type offset = dofLocation - centre;

        if((offset[0] * offset[0] + offset[1] * offset[1]) < radius * radius)
        {
          const dof_t velocity_globalDof = dof_t(&velocity, cellIter->getIndex(), dof);
          stiffness_matrix.zeroRow(velocity_globalDof, 1.0);
          const double rhs = 0.0;
          load_vector.setValues(1, &velocity_globalDof, &rhs);
        }
      }
    }
  }

  void solve(FEMatrix<cell_type>& stiffness_matrix, FEVector<cell_type>& unknown_vector, FEVector<cell_type>& load_vector, const bool usePreconditioner = true)
  {
    PETScKrylovSolver solver;
    solver.setMaxIterations(25000);
    solver.setAbsoluteTolerance(1e-4);
    solver.setRelativeTolerance(0.0);
    solver.enablePreconditioner(usePreconditioner);
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
    render(outFile);
    outFile.close();
  }

  std::pair<double, double> getVelocityVector(const std::size_t cid, const vertex_type& vertex) const
  {
    const CellVertices<dimension> vertices(m.getCoordinates(cid));
    double xVelocity(0.0), yVelocity(0.0);

    for(unsigned dof=0; dof<velocity.spaceDimension(); ++dof)
    {
      Tensor<dimension, 1> velocity_basis = velocity.evaluate_tensor(vertices, dof, vertex);
      const dof_t velocityDof = dof_t(&velocity, cid, dof);

      double velocityCoeff;
      velocity_vector.getValues(1u, &velocityDof, &velocityCoeff);

      xVelocity += velocity_basis(0).toScalar() * velocityCoeff;
      yVelocity += velocity_basis(1).toScalar() * velocityCoeff;
     }
     return std::make_pair(xVelocity, yVelocity);
  }

  void render(std::ostream& out)
  {
    std::vector< vertex<dimension> > vertices;
    std::vector< std::pair<double, double> > velocities;

    for(typename Mesh<dimension>::global_iterator vIter(m.global_begin(0)); vIter!=m.global_end(0); ++vIter)
    {
      const vertex<dimension> v(m.getVertex(vIter->getIndex()));
      vertices.push_back(v);

      const std::size_t cid = m.getContainingCell(*vIter);
      const MeshEntity localVertexEntity = m.getLocalEntity(cid, *vIter);
      const vertex<dimension> localVertex = m.getLocalCoordinate(cid, localVertexEntity.getIndex());
      velocities.push_back(getVelocityVector(cid, localVertex));
    }

    for(typename Mesh<dimension>::global_iterator eIter(m.global_begin(1)); eIter!=m.global_end(1); ++eIter)
    {
      const std::size_t cid = m.getContainingCell(*eIter);
      const std::vector<std::size_t> vertexIndices(m.getIndices(*eIter, 0));
      assert(vertexIndices.size() == 2);
      
      const MeshEntity v1Entity = m.getLocalEntity(cid, MeshEntity(0, vertexIndices[0]));
      const MeshEntity v2Entity = m.getLocalEntity(cid, MeshEntity(0, vertexIndices[1]));

      const vertex<dimension> localVertex((m.getLocalCoordinate(cid, v1Entity.getIndex()) + m.getLocalCoordinate(cid, v2Entity.getIndex()))/2.0);
      vertices.push_back(m.referenceToPhysical(cid, localVertex));
      velocities.push_back(getVelocityVector(cid, localVertex));
    }

    out << "# vtk DataFile Version 2.0" << std::endl;
    out << "Simple Navier-Stokes Solver" << std::endl;
    out << "ASCII" << std::endl;
    out << "DATASET POLYDATA" << std::endl;
    out << "POINTS " << vertices.size() << " DOUBLE " << std::endl;

    for(std::size_t point = 0; point < vertices.size(); ++point)
    {
      const vertex<dimension> v(vertices[point]);
      out << v[0] << " " << v[1] << " 0" << std::endl;
    }

    out << "POLYGONS " << m.numEntities(dimension) << " " << m.numRelations(dimension, 0)+m.numEntities(dimension)  << std::endl; 
    for(typename Mesh<dimension>::global_iterator cIter(m.global_begin(dimension)); cIter!=m.global_end(dimension); ++cIter)
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
