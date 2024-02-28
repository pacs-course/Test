/*
 * FemSolver.cpp
 *
 *  Created on: Sep 2, 2022
 *      Author: forma
 */
#include "FemSolver.hpp"
#include "Eigen/SparseLU"
#include <iostream>
#include <iomanip>

void
apsc::Fem1d::FemSolverSteady::compute() const
{
	auto n=femMesh.system_size();
	std::vector<Eigen::Triplet<double>> triplets;
	/* in theory elements may have different degree...
	 * But for simplicity a FemMesh is a mesh of elements with the same degree.
	 */
	SpMat A(n,n);
	Eigen::VectorXd rhs;
	rhs.resize(n);
	rhs.fill(0.0);
	auto num_elements=femMesh.num_elements();
	// Assemble matrix and rhs
	for (auto k=0u;k<num_elements;++k)
	{
		FemElement1D e=femMesh.elements[k];
		auto global_num=femMesh.global_numbering_map[k];
		LocalMatrix lM=
				this->stiffAssembler(e)+this->reactAssembler(e);

		LocalVector rL=this->sourceAssembler(e);
		for (auto il =0u;il<global_num.size();il++)
		{
			auto i = global_num[il];
			rhs[i] += rL[il];
			for (auto jl =0u;jl<global_num.size();jl++)
			{
				auto j = global_num[jl];
				triplets.emplace_back(i,j,lM(il,jl));
			}
		}
	}
	A.setFromTriplets(triplets.begin(),triplets.end());
	triplets.clear(); // free some memory
	triplets.shrink_to_fit();
	// Correct BC
	correctBC(A,rhs);
	// Solve
	Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int> >   solver;
	solver.analyzePattern(A);
	solver.factorize(A);
	solution_=solver.solve(rhs);
	// get the nodes for printing
	nodes_=apsc::Fem1d::get_nodes(femMesh);
}

std::ostream &
apsc::Fem1d::FemSolverSteady::print(std::ostream & stream) const
{
	stream<<"#Coordinate\t Value\n";
	for (auto i=0u;i<nodes_.size();++i)
	{
		stream << std::setprecision(10)<<std::setw(14)<<nodes_[i]<<"\t"<<solution_[i]<<"\n";
	}
	return stream;
}

void
apsc::Fem1d::FemSolverSteady::correctBC(SpMat &A, Eigen::VectorXd &rhs) const
{

	for (int ib=0;ib<2;++ib)
	{
		std::size_t i= (ib==0)?0:rhs.size()-1u;
		std::size_t e= (ib==0)?0:femMesh.num_elements()-1u;
		if (this->bc.types[ib]==apsc::Fem1d::Dirichlet)
		{
			// I DO THE SIMPLEST THING: set the row to 1 0 0...
			// I assume FemMesh is correctly ordered (this is a perequisite)
			auto nodesToZero=femMesh.global_numbering_map[e];
			//I assume the first node is 0
			for (std::size_t j=0u;j<nodesToZero.size();++j)
			{
				if (nodesToZero[j]!=i)
				{
					A.coeffRef(i,nodesToZero[j])=0.; // zero off diagonal elements
				}
				else
				{
					rhs[i]=A.coeff(i,i)*bc.h[ib];
				}
			}
			A.makeCompressed(); // not needed but better be sure
		}
		else
		{
			rhs[i]+=bc.h[ib];// add Nuunam
			if(bc.types[ib]==apsc::Fem1d::Robin)
				A.coeffRef(i,i)+=bc.c[ib];

		}
	}
	A.makeCompressed();// not needed, just to be sure.
}
