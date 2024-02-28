/*
 * FemMesh.hpp
 *
 *  Created on: Sep 1, 2022
 *      Author: forma
 */

#ifndef FEM1D_FEMMESH_HPP_
#define FEM1D_FEMMESH_HPP_
#include <vector>
#include "FemElement.hpp"
//#include <unordered_map>
#include <cmath>
#include <numbers>
namespace apsc::Fem1d
{
	/*!

	*/
	struct FemMesh
	{
		//! all mesh elements
		std::vector<FemElement1D> elements;
		//! elements -> global node number
		std::vector<std::vector<std::size_t>> global_numbering_map;
		/*!
		 Returns the number of elements
		 @return The number f elements in the mesh			
		*/
		std::size_t num_elements() const{return elements.size();}
		unsigned int degree;
		std::size_t system_size() const
		{
			// Mesh is assumed to be ordered, so the size is
			// 1 + number of last node
			return global_numbering_map.back().back()+1u;
		}
	};

	enum MeshType{Lagrange=0,Chebishev=1};

	template <class MeshVertex>
	FemMesh
	make_FemMesh(unsigned int degree,
			MeshVertex const & meshVertex,
			MeshType type=Lagrange)
	{
		FemMesh res;
		res.degree=degree;
		std::size_t localDofNumber=degree+1;
		std::vector<double> localNodes;
		std::vector<std::size_t> localGlobalNum;
		localNodes.resize(localDofNumber);
		localGlobalNum.resize(localDofNumber);
		std::size_t number_of_elements=meshVertex.size()-1u;
		res.elements.reserve(number_of_elements);
		res.global_numbering_map.reserve(number_of_elements);
		std::size_t node_counter=0u;
		for (std::size_t i=0;i<meshVertex.size()-1u;++i)
		{
			auto left  = meshVertex[i];
			auto right = meshVertex[i+1];
			auto h = (right -left)/degree;
			if (type==Lagrange)
			{
				for (auto k=0u; k<degree;++k)
				{
					localNodes[k]=left +k*h;
					localGlobalNum[k]=node_counter++;
				}
			}
			else
			{
				constexpr double pi=std::numbers::pi_v<double>;
				localNodes.front()=left;
				localGlobalNum[0]=node_counter++;
				for (auto k=1u; k<degree;++k)
				{
					localNodes[k]=(left+right)/2.
							-h*std::cos(pi*k/static_cast<double>(degree))/2.;
					localGlobalNum[k]=node_counter++;
				}
			}
			localNodes.back()=right;
			localGlobalNum.back()=node_counter;
			res.elements.emplace_back(localNodes);
			res.global_numbering_map.emplace_back(localGlobalNum);
		}
		return res;
	}
	/*!
	 * Since in principle we may have elements with different degree in the same mesh
	 * for generality I create a simple general tool to extract the global nodes coordinates
	 * from a FemMesh
	 * @tparam FEMMESH A FemMesh type, it provides
	 * @param femMesh
	 * @return the vector with the (global) nodes coordinates
	 */
	inline auto get_nodes(FemMesh const & femMesh)
	{
		std::vector<double> nodes;
		nodes.resize(femMesh.system_size());
		for (auto i=0u;i<femMesh.num_elements();++i)
		{
			auto const & nl = femMesh.global_numbering_map[i];
			auto const & e =femMesh.elements[i];
			for (auto j=0u;j<e.nodes_.size();++j)
			{
				auto gn = nl[j];
				auto coord = e[j];
				nodes[gn]=coord;
			}
		}
		return nodes;
	}
}




#endif /* FEM1D_FEMMESH_HPP_ */
