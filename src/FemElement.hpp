/*
 * FemElement.hpp
 *
 *  Created on: Aug 30, 2022
 *      Author: forma
 */

#ifndef FEM1D_FEMELEMENT_HPP_
#define FEM1D_FEMELEMENT_HPP_
#include <vector>
#include <array>
#include <memory>
#include <functional>
#include "Lagrange1D.hpp"
#include "Gauss_rule.hpp"
#include "FemTraits.hpp"
#include "Eigen/Core"
namespace apsc::Fem1d
{
	/*!
     * @brief 1D finite element class. Straight from the book
	 *
	 * @todo be able to specify different quadrules for the different operators.
	 * In fact we can separate the quadrature rule from the finite element class by specifying it
	 * as parameter in the methds that build the local matrices/source term
	 */
	class FemElement1D: public Lagrange1D, public apsc::Fem1d::FemTraits
	{
	public:
		/*!
		 * Constructor
		 * @param nodes a vector containing the local nodes coordinates.
		 * @param quad A quadratture rule following the API of NumericalIntegration::StandardQuadratureRule
		 */
		FemElement1D(std::vector<double> const & nodes)
		{
			nodes_=nodes;
		}
		/*!
		 * Not much uses, you can state the local values.
		 * @param v a vector of values
		 */
		void set_values(std::vector<double> const & v){
			values_=v;
		}
		/*!
		 * Sets the local nodes
		 * @param nodes
		 */
		void set_nodes(std::vector<double> const & nodes)
		{
			nodes_=nodes;
		}

		/*!
		 * Returns the ith node
		 * @param i index
		 * @return the node coordinate
		 */
		double operator[](std::size_t i) const {return nodes_[i];}
		/*!
		 * Returns the ith node
		 * @param i index
		 * @return the node coordinate
		 */
		double & operator[](std::size_t i){return nodes_[i];}

		/*!
		 * Compute local mass matrix
		 */
		LocalMatrix mass(QuadRule const & quadRule)const ;
		/*!
		 * Compute Stiffness matrix
		 * \f[
		 * \int \mu \psi_i^\prime\psi_j^\prime
		 * \f]
		 * @param mu The \f$\mu \f$ constant coefficient
		 * @return The local stiffness matrix
		 */
		LocalMatrix stiff(double const & mu,QuadRule const & quadRule) const;
		/*!
		 * Compute Stiffness matrix
		 * @param mu The \f$\mu \f$ function coefficient
		 * @return The local stiffness matrix
		 */
		LocalMatrix stiff(Function const & mu,QuadRule const & quadRule)const ;
		/*!
		 * Compute reaction matrix
		 ** \f[
		 * \int c \psi_i\psi_j
		 * \f]
		 * @param c  The \f$ c \f$ constant coefficient
		 * @return The local reaction matrix
		 */
		LocalMatrix react(double const & c,QuadRule const & quadRule)const ;
		/*!
		 * Compute reaction matrix
		 *\f[
		 * \int c \psi_i\psi_j
		 * \f]
		 * @param c  The \f$ c \f$ function coefficient
		 * @return The local reaction matrix
		 */
		LocalMatrix react(Function const & c,QuadRule const & quadRule)const ;
		/*!
		 * Compute source term
		 *
		 *\f[
		 * \int s \psi_i
		 * \f]
		 * @param s  The \f$ s \f$ constant source
		 * @return The local source term
		 */
		LocalVector source(double const & s,QuadRule const & quadRule)const ;
		/*!
		 * Compute source term
		 *
		 *\f[
		 * \int s \psi_i
		 * \f]
		 * @param s  The \f$ s \f$ function source
		 * @return The local source term
		 */
		LocalVector source(Function const & s,QuadRule const & quadRule)const ;
		/*!
		 * The element size
		 * @return the element size
		 */
		double h()const{ return nodes_.back()-nodes_.front();}
		/*!
		 * The element polynomial order
		 * @return The order
		 */
		int p()const{ return nodes_.size()-1;}
		/*!
		 * A helper function to compute integral of a function over the finite element
		 * @param f function
		 * @return The integral
		 */
		double compute_integral(Function const & f, QuadRule const & quadRule)const
		{
		  return quadRule.apply(f, nodes_.front(), nodes_.back());
		}
	private:
	};
}



#endif /* FEM1D_FEMELEMENT_HPP_ */
