#include<dolfin.h>
#include<SteadyStokes.h>
#include<cstdlib>
#include<iostream>

using namespace dolfin;

class VelocityXBC : public Function
{
public:
  VelocityXBC(Mesh& mesh) :  Function(mesh)
  {
  }

  void eval(real* values, const real* x) const
  {
    (*values) = x[0]*x[0] + x[1]*x[1];
  }
};

class ZeroBC : public Function
{
public:
  ZeroBC(Mesh& mesh) :  Function(mesh)
  {
  }

  void eval(real* values, const real* x) const
  {
    (*values) = 0.0;
  }
};

class Boundary : public SubDomain
{
  bool inside(const real* x, bool on_boundary) const
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

  bool inside(const real* x, bool on_boundary) const
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

  bool inside(const real* x, bool on_boundary) const
  {
    return on_boundary && (std::abs(x[1]) < EPSILON || std::abs(x[1]-1.0) < EPSILON);
  }
};

int main(int argc, char* argv[])
{
  // Create Mesh
  UnitSquare mesh(10,10);

  // Create functions for boundary conditions
  VelocityXBC velocityXBC(mesh);
  ZeroBC zeroBC(mesh);

  // Functions defining boundaries
  Boundary boundary;
  LeftRightBoundary leftRightBoundary;
  TopBottomBoundary topBottomBoundary;
  
  // Define sub systems for boundary conditions
  SubSystem velocity_x(0, 0);
  SubSystem velocity_y(0, 1);
  SubSystem pressure(1);

  // Velocity boundary condition
  DirichletBC velocity_x_bc(velocityXBC, mesh, boundary, velocity_x);
  DirichletBC velocity_y_bc(zeroBC, mesh, leftRightBoundary, velocity_y);
  DirichletBC pressure_bc(zeroBC, mesh, topBottomBoundary, pressure);

  // Set up PDE
  Function f(mesh, 0.0);
  SteadyStokesBilinearForm a;
  SteadyStokesLinearForm L(f);
  Array<BoundaryCondition*> bcs(&velocity_x_bc, &velocity_y_bc, &pressure_bc);
  LinearPDE pde(a, L, mesh, bcs);
  
  // Solve PDE
  Function u;
  Function p;
  pde.set("PDE linear solver", "direct");
  pde.solve(u, p);

  // Plot solution
  plot(u);
  plot(p);

  // Save solution in VTK format
  File ufile_pvd("velocity.pvd");
  ufile_pvd << u;
  File pfile_pvd("pressure.pvd");
  pfile_pvd << p;

  return EXIT_SUCCESS;
}
