using namespace cfd;

int main(int argc, char* argv[])
{
  PETScManager::instance().init(argc, argv);

  typedef TriangularCell cell_type;
  TriangularMeshBuilder meshBuilder(3.0, 1.0, 1.0/900.0);
  meshBuilder.addPolygon(Polygon(vertex<2>(0.5, 0.5), 16, 0.148, 0), 5);
  Mesh<cell_type::dimension> m(meshBuilder.buildMesh());

  FESystem<cell_type::dimension> system(m);
  Element velocity = system.addElement(new LagrangeTriangleQuadratic());
  Element pressure = system.addElement(new LagrangeTriangleLinear());

  FunctionSpace velocitySpace = FunctionSpace(velocity, m);
  FunctionSpace pressureSpace = FunctionSpace(pressure, m);
  FunctionSpace coupledSpace = velocitySpace + pressureSpace;

  m.addField("velocity", velocitySpace);
  m.addField("pressure", pressureSpace);

  Problem timestep = buildNavierStokesSolve();

  for(int i=0; i<6000; ++i)
  {
    timestep.execute();
    std::stringstream filename;
    filename << "./navier_stokes_" << boost::format("%|04|") % i << ".vtk";
    system.outputToFile(filename.str());
  }
}

Problem buildNavierStokesSolve()
{
  Problem p = system.newProblem();
  Field velocityField = p.getField("velocity");
  Field pressureField = p.getField("pressure");

  Operator systemMatrix(coupledSpace, coupledSpace);
  systemMatrix +=
    B(velocity, velocity)*dx +
    B(scalar(theta * k * kinematic_viscosity) * grad(velocity), grad(velocity))*dx +
    B(scalar(theta * k * kinematic_viscosity * -1.0) * inner(grad(velocity), n), velocity)*ds +
    B(scalar(-1.0 * k) * pressure, div(velocity))*dx +
    B(div(velocity), pressure)*dx;

  Operator nonLinearRhs(velocitySpace, velocitySpace);
    nonLinearRhs +=
    B(velocity, velocity)*dx +
    B(scalar(-(1.0-theta) * k * kinematic_viscosity) * grad(velocity), grad(velocity))*dx +
    B(scalar((1.0-theta) * k * kinematic_viscosity) * inner(grad(velocity), n), velocity)*ds +
    B(scalar(-(1.0-theta)*k) * inner(velocityField, grad(velocity)), velocity)*dx;

  Field velocityRhs = nonLinearRhs * velocityField;
  Field load(coupledSpace);

  TemporalIndex n;
  Scalar residual(n);
  Field unknownVelocity(velocitySpace, n);
  Field unknownGuess(coupledSpace, n);
  Field coupledSolution(coupledSpace, n);

  unknownGuess[0] = velocityField;
  unknownGuess[n] = coupledSolution[n-1];
  unknownVelocity[n] = unknown[n];

  Operator linearisedSystem[n] = systemMatrix;
  linearisedSystem[n] += B(scalar(theta*k) * inner(unknownVelocity[n], grad(velocity)), velocity)*dx;

  /* boundary condition magic */

  residual[n] = ((linearisedSystem[n] * unknownGuess[n]) - load).two_norm();
  coupledSolution[n] = linear_solve(linearisedSystem[n], unknownGuess[n], load);

  n.setTermination(residual[n] < 1e-3);

  p.setNewValue("velocity", coupledSolution[final].subField(velocitySpace));
  p.setNewValue("pressure", coupledSolution[final].subField(pressureSpace));

  p.finish();
  return p;
}
