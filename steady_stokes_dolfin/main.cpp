#include<dolfin.h>
#include<SteadyStokes.h>
#include<cstdlib>
#include<iostream>

using namespace dolfin;

class VelocityXBC : public Function
{
public:
  VelocityXBC(Mesh& mesh)
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
  ZeroBC(Mesh& mesh)
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
  UnitSquare mesh(50,50);

  // Create functions for boundary conditions
  VelocityXBC velocityXBC(mesh);
  ZeroBC zeroBC(mesh);

  // Functions defining boundaries
  Boundary boundary;
  LeftRightBoundary leftRightBoundary;
  TopBottomBoundary topBottomBoundary;
  
  // Define sub systems for boundary conditions
  SteadyStokesFunctionSpace mixedSpace(mesh);
  SubSpace velocity(mixedSpace, 0);
  SubSpace velocity_x(velocity, 0);
  SubSpace velocity_y(velocity, 1);
  SubSpace pressure(mixedSpace, 1);

  // Velocity boundary condition
  DirichletBC velocity_x_bc(velocity_x, velocityXBC, boundary, topological);
  DirichletBC velocity_y_bc(velocity_y, zeroBC, leftRightBoundary, topological);
  DirichletBC pressure_bc(pressure, zeroBC, topBottomBoundary, topological);

  // Set up PDE
  Constant f(0.0);
  SteadyStokesBilinearForm a(mixedSpace, mixedSpace);
  SteadyStokesLinearForm L(mixedSpace);
  L.f = f;
  Array<BoundaryCondition*> bcs(&velocity_x_bc, &velocity_y_bc, &pressure_bc);
  VariationalProblem pde(a, L, bcs);
  
  // Solve PDE
  Function w;
  //Function u;
  //Function p;
  pde.set("PDE linear solver", "direct");
  pde.solve(w);
  Function u = w[0];
  Function p = w[1];

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
