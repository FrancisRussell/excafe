#include <dolfin.h>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "SteadyStokes.h"

using namespace dolfin;

class VelocityXBC : public Expression
{
public:
  VelocityXBC(Mesh& mesh) : Expression()
  {
  }

  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = x[0]*x[0] + x[1]*x[1];
  }
};

class ZeroBC : public Expression
{
public:
  ZeroBC(Mesh& mesh) : Expression()
  {
  }

  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = 0.0;
  }
};

class Boundary : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return on_boundary;
  }
};

class LeftRightBoundary : public SubDomain
{
private:
  const double EPSILON;
public:
  LeftRightBoundary() : EPSILON(1e-9)
  {
  }

  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return on_boundary && (std::abs(x[0]) < EPSILON || std::abs(x[0]-1.0) < EPSILON);
  }
};

class TopBottomBoundary : public SubDomain
{
private:
  const double EPSILON;
public:
  TopBottomBoundary() : EPSILON(1e-9)
  {
  }

  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return on_boundary && (std::abs(x[1]) < EPSILON || std::abs(x[1]-1.0) < EPSILON);
  }
};

int main(int argc, char* argv[])
{
  // Create Mesh
  UnitSquareMesh mesh(50,50);

  // Create functions for boundary conditions
  VelocityXBC velocityXBC(mesh);
  ZeroBC zeroBC(mesh);

  // Functions defining boundaries
  Boundary boundary;
  LeftRightBoundary leftRightBoundary;
  TopBottomBoundary topBottomBoundary;

  // Define sub systems for boundary conditions
  SteadyStokes::FunctionSpace mixedSpace(mesh);
  SubSpace velocity(mixedSpace, 0);
  SubSpace velocity_x(velocity, 0);
  SubSpace velocity_y(velocity, 1);
  SubSpace pressure(mixedSpace, 1);

  // Velocity boundary condition
  DirichletBC velocity_x_bc(velocity_x, velocityXBC, boundary, "topological");
  DirichletBC velocity_y_bc(velocity_y, zeroBC, leftRightBoundary, "topological");
  DirichletBC pressure_bc(pressure, zeroBC, topBottomBoundary, "topological");

  // Set up PDE
  Constant f(0.0, 0.0);
  SteadyStokes::BilinearForm a(mixedSpace, mixedSpace);
  SteadyStokes::LinearForm L(mixedSpace);
  L.f = f;

  std::vector<const BoundaryCondition*> bcs;
  bcs.push_back(&velocity_x_bc);
  bcs.push_back(&velocity_y_bc);
  bcs.push_back(&pressure_bc);

  Function w(mixedSpace);
  LinearVariationalProblem pde(a, L, w, bcs);
  
  // Solve PDE
  LinearVariationalSolver solver(pde);
  solver.solve();
  Function u = w[0];
  Function p = w[1];

  // Plot solution
  const boost::shared_ptr<VTKPlotter> velocityPlotter = plot(u);
  interactive(true);
  const boost::shared_ptr<VTKPlotter> pressurePlotter = plot(p);
  interactive(true);

  // Save solution in VTK format
  File ufile_pvd("velocity.pvd");
  ufile_pvd << u;
  File pfile_pvd("pressure.pvd");
  pfile_pvd << p;

  return EXIT_SUCCESS;
}
