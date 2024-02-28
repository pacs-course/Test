/*
 * FemSolver.hpp
 *
 *  Created on: Aug 30, 2022
 *      Author: forma
 */

#ifndef FEM1D_FEMSOLVER_HPP_
#define FEM1D_FEMSOLVER_HPP_
#include "FemElement.hpp"
#include "FemMesh.hpp"
#include "FemTraits.hpp"
#include "Eigen/SparseCore"
#include <variant>
#include <iosfwd>
#include <tuple>
#include "FemData.hpp"
namespace apsc::Fem1d
{
/*!
 * Equation
 * \f[
 * - (\mu u^{\prime}^\prime + a u =s
 * \f]
 * where \f$ \mu \f$, \f$ a \f$ and \f$ s \f$ can be either a value or a function
 *
 */
	class FemSolverSteady: public FemTraits
	{
	// A misplaced comment
	public:
		using SpMat=Eigen::SparseMatrix<double,Eigen::ColMajor>;
		FemSolverSteady() = default;
		/*!
		 * Set the solver for
		 *
		 * \f[
		 * - (\mu u^{\prime})^\prime + a u =s, \quad a <x <b
		 * \f]
		 * where \f$ \mu \f$, \f$ a \f$ and \f$ s \f$ can be either a value or a function
		 *
		 * The boudnary conditions have the general form
		 *
		 * Dirichlet
		 * \f[
		 * u=h
		 * \f]
		 * Neuman
		 * \f[
		 * \mu u^\prime(a)=h_a, \quad  -\mu u^\prime(b) = h in x=b
		 * \f]
		 * Robin   mu u' + c u=h in x=a -mu u' + cu =h in x=b
		 *
		 * @param mesh The mesh holding the finite elements
		 * @param bc The boundary conditions
		 * @param mu The coefficient for the diffusion operator. A double or a (const double &) -> double
		 * @param s The coefficient for the source term. A double or a (const double &) -> double
		 * @param a The coefficient for the reaction operator. A double or a (const double &) -> double
		 */
		FemSolverSteady(FemMesh const & mesh, Bc const & bc, QuadRule const & rule,
				CoefficientHolder mu, CoefficientHolder s, CoefficientHolder a={0.0}):
			femMesh{mesh}, bc{bc}, stiffAssembler{rule,mu}, reactAssembler{rule,a},sourceAssembler{rule,s} {}

		//! @todo Allows for different rules
		//	FemSolverSteady(FemMesh const & mesh, Bc const & bc, QuadRule const & rule,
		//			CoefficientHolder mu, CoefficientHolder s, CoefficientHolder a={0.0}):
		//				femMesh{mesh}, bc{bc}, stiffAssembler{rule,mu}, reactAssembler{rule,a},sourceAssembler{rule,s} {}

		void
		setA(const CoefficientHolder &a)
		{
			reactAssembler.set_coeff(a);
		}

		void
		setMu(const CoefficientHolder &mu)
		{
			stiffAssembler.set_coeff(mu);
		}

		void
		setS(const CoefficientHolder &s)
		{
			sourceAssembler.set_coeff(s);
		}

		void
		setFemMesh(const FemMesh &femMesh)
		{
			this->femMesh = femMesh;
		}

		auto
		nodesAndSolution() const
		{
			return std::make_tuple(this->nodes_,this->solution_);
		}

		void compute() const;

		void
		setBc(const Bc &bc)
		{
			this->bc = bc;
		}

		std::ostream &
		print(std::ostream & stream) const;
	private:
		FemMesh  femMesh;
		Bc bc;
		AssembleStiff stiffAssembler;
		AssembleReact reactAssembler;
		AssembleSource sourceAssembler;
		mutable Eigen::VectorXd solution_;
		mutable std::vector<double> nodes_;
		void correctBC(SpMat & A, Eigen::VectorXd & rhs) const;
	};

} // namespace apsc::Fem1d

#endif /* FEM1D_FEMSOLVER_HPP_ */
