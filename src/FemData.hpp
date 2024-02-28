/*
 * Fem1D_Data.hpp
 *
 *  Created on: Aug 30, 2022
 *      Author: forma
 */

#ifndef FEM1D_FEMDATA_HPP_
#define FEM1D_FEMDATA_HPP_
#include "FemElement.hpp"
#include "FemTraits.hpp"
#include <memory>

namespace apsc::Fem1d{


/*!
 * @enum BCType Boundary condition type
 */
	enum BCType{Dirichlet=0, Neumann=1, Robin};

/*!
 * It holds the boundary
 *
 * Dirichlet u=h
 * Neuman  mu u'=h in x=a  -mu u' = h in x=b
 * Robin   mu u' + c u=h in x=a -mu u' + cu =h in x=b

 */
	struct Bc
	{
		//! The type
		std::array<BCType,2> types;
		//! The coefficient h
		std::array<double,2> h={0.,0.};
		//! The coefficient c
		std::array<double,2> c={0.,0.};
	};

struct Coefficients
{
		/*!
		 * Variants
		 */
		FemTraits::CoefficientHolder mu;
		FemTraits::CoefficientHolder a;
		FemTraits::CoefficientHolder s;

};
	/*!
	 * The base class for the basic operations to construct local matrices
	 *
	 * With this trick we can handle the fact the the coefficients can be a double
	 * or a function double-> double without the need of testing the variant insside tha main loop over the elements.
	 *
	 * We use CRTP technique
	 * @tparam Derived The actual operation
	 */
	template<class Derived>
	class BasicMatrixOperator: public FemTraits
	{
	public:
		/*!
		 * We need to pass the coefficient in the forma of a std::variant that encapsulates
		 * a double or a std::function<double (const double &)>. Depending on the actual content
		 * the operator will call the appropriate overload.
		 *
		 * @note this wrapper around the actual operator to construct the local matrix has the role
		 * to avoid testing the actual content of the variant inside the elemental loop, thus increasing performance.
		 *
		 * @param coeff The coefficient. It is a variant that may take a double or a std::function
		 */
		BasicMatrixOperator(QuadRule const & rule, FemTraits::CoefficientHolder coeff={0.0}):coeff_(coeff),
				quadRule_(rule.clone()){};
/*		BasicMatrixOperator(BasicMatrixOperator const & other):
			coeff_{other.coeff_},quadRule_{other.quadRule_->clone()},f{other.f}{};

		BasicMatrixOperator & operator =(BasicMatrixOperator const & other)
		{
			coeff_=other.coeff_;
			quadRule_=other.quadRule_->clone();
			f=other.f;
		}
*/
		LocalMatrix operator()(FemElement1D const & fem) const
		{
			return as_leaf().operator()(fem);
		}
		void set_coeff(FemTraits::CoefficientHolder coeff){coeff_=coeff;}
		void set_rule(QuadRule const & rule){quadRule_=rule.clone();}
	protected:
		//! The classic CRTP trick
		Derived & as_leaf(){return static_cast<Derived &>(*this);}
		//! The classic CRTP trick
		Derived const & as_leaf()const {return static_cast<const Derived &>(*this);}
		//! The coefficient (used only by some of the derived classes)
		FemTraits::CoefficientHolder coeff_;
		apsc::PointerWrapper<QuadRule> quadRule_;
		/*! The actual function to call.
		 * Heeded to resolve the overload. Used only by some of the derived classes.
		 */
		std::function<LocalMatrix(FemElement1D const &)> f;
	};

	/*!
	 * Encapsulates the null operator.
	 *
	 * It may be used when a coefficient is zero, avoiding the explicit computation
	 * of the array
	 */
	struct NullMatrixOperator: public BasicMatrixOperator<NullMatrixOperator>
	{
		using BasicMatrixOperator<NullMatrixOperator>::BasicMatrixOperator;
		LocalMatrix operator()(FemElement1D const & fem) const
		{
			int siz=fem.nodes_.size();
			return FemTraits::LocalMatrix::Zero(siz,siz);
		}
	};

	/*!
	 * Assembles the local mass matrix
	 */
	struct AssembleMass: public BasicMatrixOperator<AssembleMass>
	{
		using BasicMatrixOperator<AssembleMass>::BasicMatrixOperator;
		/*!
		 * The operator that computes the local mass matrix
		 * @param fem A finite element
		 * @return The mass matrix
		 */
		LocalMatrix operator()(FemElement1D const & fem) const
		{
			return fem.mass(*quadRule_);
		}

	};

	/*!
	 * Assembles the stiffness matrix.
	 *
	 * Here we see the use of the internal function to resolve the overload.
	 *
	 */
	struct AssembleStiff: public BasicMatrixOperator<AssembleStiff>
	{
		//using Function=FemTraits::Function;
		//using LocalMatrix=FemTraits::LocalMatrix;
		/*!
		 * Constructor
		 *
		 * according to the value held by the coefficient we choose the correct form
		 *
		 * @param coeff The coefficient
		 */
		AssembleStiff(QuadRule const & rule, FemTraits::CoefficientHolder coeff={0.0}):
			BasicMatrixOperator<AssembleStiff>(rule,coeff)
		{
			if (std::holds_alternative<double>(coeff_))
			{
				if (std::get<double>(coeff_)!=0)
				f=[this](FemElement1D const & fem)
									{
					return fem.stiff(std::get<double>(this->coeff_),*quadRule_);
									};
				else
					f=[this](FemElement1D const & fem)
														{
										return FemTraits::LocalMatrix::Zero(fem.nodes_.size(),fem.nodes_.size());
														};
			}
			else
			{
				f=[this](FemElement1D const & fem)
										{
					return fem.stiff(std::get<Function>(this->coeff_),*quadRule_);
										};
			}
		}
		/*!
		 * The operator that computes the local stiffness matrix
		 * @param fem A finite element
		 * @return The stiffness matrix
		 */
		LocalMatrix operator()(FemElement1D const & fem) const
		{
			return f(fem);
		}
	};

	/*!
	 * Assembles the reaction term
	 *
	 * Here we see the use of the internal function to resolve the overload.
	 *
	 * */

	struct AssembleReact: public BasicMatrixOperator<AssembleStiff>
	{
		//using Function =FemTraits::Function;
		//using LocalMatrix=FemTraits::LocalMatrix;
		AssembleReact(QuadRule const & rule,FemTraits::CoefficientHolder coeff={1.0}):
			BasicMatrixOperator(rule, coeff)
		{
			if (std::holds_alternative<double>(coeff_))
			{
				if (std::get<double>(coeff_)!=0.0)
				f=[this](FemElement1D const & fem)
									{
					return fem.react(std::get<double>(this->coeff_),*quadRule_);
									};
				else
					f=[this](FemElement1D const & fem)
														{
										return FemTraits::LocalMatrix::Zero(fem.nodes_.size(),fem.nodes_.size());
														};
			}
			else
			{
				f=[this](FemElement1D const & fem)
										{
					return fem.react(std::get<Function>(this->coeff_),*quadRule_);
										};
			}
		}
		/*!
		 * The operator that computes the local reaction term matrix
		 * @param fem A finite element
		 * @return The reaction term matrix
		 */
		LocalMatrix operator()(FemElement1D const & fem) const
		{
			return f(fem);
		}
	};

	/*!
	 * Assembles the source term
	 *
	 * Here we see the use of the internal function to resolve the overload.
	 *
	 * */
	struct AssembleSource: public FemTraits
	{
		AssembleSource(QuadRule const & rule, FemTraits::CoefficientHolder coeff={1.0}):
			coeff_{coeff},quadRule_{rule.clone()}
		{
			if (std::holds_alternative<double>(coeff_))
			{
				if(std::get<double>(coeff_)!=0.0)
				f=[this](FemElement1D const & fem)
									{
					return fem.source(std::get<double>(this->coeff_),*quadRule_);
									};
				else
					f=[this](FemElement1D const & fem)
														{
										return FemTraits::LocalVector::Zero(fem.nodes_.size());
														};

			}
			else
			{
				f=[this](FemElement1D const & fem)
										{
					return fem.source(std::get<Function>(this->coeff_),*quadRule_);
										};

			}
		}
		/*!
		 * The operator that computes the local reaction term matrix
		 * @param fem A finite element
		 * @return The reaction term matrix
		 */
		LocalVector operator()(FemElement1D const & fem) const
		{
			return f(fem);
		}
		void set_coeff(FemTraits::CoefficientHolder coeff){coeff_=coeff;}
		void set_rule(QuadRule const & rule){quadRule_=rule.clone();}
	private:
		FemTraits::CoefficientHolder coeff_;
		apsc::PointerWrapper<QuadRule> quadRule_;
		/*! The actual function to call.
		 * Needed to resolve the overload. Used only by some of the derived classes.
		 */
		std::function<LocalVector (FemElement1D const &)> f;

	};


	/*!
	 * Encapsulates the null operator for the source term.
	 *
	 * It may be used when a coefficient is zero, avoiding the explicit computation
	 * of the vector
	 */
	struct NullVectorOperator: public FemTraits
	{
	public:
		/*
		 * Constuctor takes a coefficient only for compatibility.
		 */
		void set_coeff(FemTraits::CoefficientHolder){};
		void set_rule(QuadRule const &){};

		NullVectorOperator(QuadRule const &, FemTraits::CoefficientHolder){};

			LocalVector operator()(FemElement1D const & fem) const
			{
				return FemTraits::LocalVector::Zero(fem.nodes_.size());
			}
	};

}// end namespace
#endif

