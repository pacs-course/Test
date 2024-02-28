/*
 * Fem1dTraits.hpp
 *
 *  Created on: Aug 30, 2022
 *      Author: forma
 */

#ifndef FEM1D_FEMTRAITS_HPP_
#define FEM1D_FEMTRAITS_HPP_
#include <variant>
#include "Gauss_rule.hpp"
#include "Eigen/Core"

namespace apsc::Fem1d
{
	struct FemTraits
	{
		using QuadRule=apsc::NumericalIntegration::QuadratureRuleBase;
		using LocalMatrix = Eigen::MatrixXd;
		using Function = std::function<double (const double &)>;
		using LocalVector = Eigen::Matrix<double,Eigen::Dynamic,1>;
		using CoefficientHolder = std::variant<double,Function>;
	};
}



#endif /* FEM1D_FEMTRAITS_HPP_ */
