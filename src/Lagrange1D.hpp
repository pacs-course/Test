/*
 * Lagrange1D.hpp
 *
 *  Created on: Aug 30, 2022
 *      Author: forma
 */

#ifndef FEM1D_LAGRANGE1D_HPP_
#define FEM1D_LAGRANGE1D_HPP_
#include <array>
#include <vector>
#include <algorithm>
namespace apsc::Fem1d
{
	struct Lagrange1D
	{
		/*!
		*/
		double operator()(std::size_t i, double const & x) const
		{
			double res=1.0;
			auto n=nodes_.size();
			for (std::size_t j=1u;j<n;++j)
			{
				auto k = (i+j) % n;
				res *=(x-nodes_[k])/(nodes_[i]-nodes_[k]);
			}
			return res;
		}
		double der(std::size_t i, double const & x) const
		{
			double res=0.0;
			auto n=nodes_.size();
			for (std::size_t l=1u;l<n;++l)
			{
				auto j = (i+l) % n;
				auto fact =1.0/(nodes_[i]-nodes_[j]);
				for (std::size_t k=0u;k<n;++k)
				{
					if(k==i or k==j) continue;
					fact *=(x-nodes_[k])/(nodes_[i]-nodes_[k]);
				}
				res+=fact;
			}
			return res;
		}
		double der(double const & x) const
		{
			double res=0.0;
			auto n=nodes_.size();
			for (std::size_t i=0u;i<n;++i)
			{
				res+=this->der(i,x)*values_[i];
			}
			return res;
		}
		double operator()(double const & x) const
		{
			double res=0.0;
			auto n=nodes_.size();
			for (std::size_t i=0u;i<n;++i)
			{
				res+=this->operator()(i,x)*values_[i];
			}
			return res;
		}

		std::vector<double> nodes_;
		std::vector<double> values_;
	};


}



#endif /* FEM1D_LAGRANGE1D_HPP_ */
