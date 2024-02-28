/*
 * main_testLMat.cpp
 *
 *  Created on: Aug 30, 2022
 *      Author: forma
 */
#include "Eigen/Core"
#include "FemElement.hpp"
#include "FemSolver.hpp"
#include <iostream>
#include <fstream>
#include <cmath>

#include "FemData.hpp"
int
main()
{
  using namespace apsc::Fem1d;
  // P1
  std::cout << " P1\n";
  {
	  std::vector<double> nodesP1 = {0, 1};

	  FemElement1D p1(nodesP1);
	  std::cout << " Mass\n";
	  AssembleMass massAssembler{apsc::NumericalIntegration::GaussLegendre2p{},1.0};
	  std::cout << massAssembler(p1);
	  std::cout << std::endl;
	  std::cout << " Stiff\n";
	  AssembleStiff stiffAssembler{apsc::NumericalIntegration::GaussLegendre2p{},1.0};
	  std::cout << stiffAssembler(p1);
	  std::cout << std::endl;
	  std::cout << " React\n";
	  AssembleReact reactAssembler{apsc::NumericalIntegration::GaussLegendre2p{},
		  [](const double &x) { return 1.0; }};
	  std::cout << reactAssembler(p1);
	  std::cout << std::endl;
	  std::cout << " Source\n";
	  AssembleSource sourceAssembler{apsc::NumericalIntegration::GaussLegendre2p{},[](const double &x) { return x * x; }};
	  std::cout << sourceAssembler(p1);
	  std::cout << std::endl;

  }
  std::cout << " P3\n";

  std::vector<double> nodesP3 = {0., 1./3, 2./3, 1.};
  {
	  FemElement1D p1(nodesP3);

	  std::cout << " Mass\n";
	  AssembleMass massAssembler{apsc::NumericalIntegration::GaussLegendre3p{},1.0};
	  std::cout << massAssembler(p1);
	  std::cout << std::endl;
	  std::cout << " Stiff\n";
	  AssembleStiff stiffAssembler{apsc::NumericalIntegration::GaussLegendre3p{},1.0};
	  std::cout << stiffAssembler(p1);
	  std::cout << std::endl;
	  std::cout << " React\n";
	  AssembleReact reactAssembler{apsc::NumericalIntegration::GaussLegendre3p{},[](const double &x) { return 1.0; }};
	  std::cout << reactAssembler(p1);
	  std::cout << std::endl;
	  std::cout << " Source\n";
	  AssembleSource sourceAssembler{apsc::NumericalIntegration::GaussLegendre3p{},[](const double &x) { return 1.0; }};
	  std::cout << sourceAssembler(p1);
	  std::cout << std::endl;

  }

  std::cout << " Fem mesh\n";
  std::vector<double>  meshVertex = {0, 0.1, 0.2, 0.5, 0.6, 0.8, 0.9, 1.0};
  apsc::Fem1d::FemMesh femMesh = apsc::Fem1d::make_FemMesh(5, meshVertex,apsc::Fem1d::Chebishev);

  Bc bc;
  bc.types={Dirichlet,Neumann};
  bc.h={0.0,0.0};
  bc.c={0.0,1.0};
  double mu0=0.5;
  auto mu = [mu0](double const & x){return mu0+std::abs(x*x*x);};
  FemSolverSteady solve(femMesh,bc,apsc::NumericalIntegration::GaussLegendre4p{},1.0,mu,1.0);
  solve.compute();
  solve.print(std::cout);
  std::ofstream file("result.dat");
  solve.print(file);
  file.close();
  return 0;
}
