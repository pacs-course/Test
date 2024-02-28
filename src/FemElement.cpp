#include "FemElement.hpp"

apsc::Fem1d::FemElement1D::LocalMatrix
apsc::Fem1d::FemElement1D::mass(QuadRule const & quadRule) const
{
  using namespace apsc::Fem1d;
  // A different comment
  FemElement1D::LocalMatrix m{nodes_.size(), nodes_.size()};
  for(std::size_t i = 0u; i < nodes_.size(); ++i)
    {
      for(std::size_t j = 0u; j <= i; ++j)
        {
          FemElement1D::Function f = [i, j, this](const double &x) {
            auto const &self = *this;
            return self(i, x) * self(j, x);
          };
          m(i, j) = this->compute_integral(f,quadRule);
          if(i != j)
            m(j, i) = m(i, j);
        }
    }
  return m;
}

apsc::Fem1d::FemElement1D::LocalMatrix
apsc::Fem1d::FemElement1D::stiff(const double &mu,QuadRule const & quadRule) const
{
  using namespace apsc::Fem1d;
  FemElement1D::LocalMatrix m{nodes_.size(), nodes_.size()};
  for(std::size_t i = 0u; i < nodes_.size(); ++i)
    {
      for(std::size_t j = 0u; j <= i; ++j)
        {
          FemElement1D::Function f = [i, j, this](const double &x) {
            auto const &self = *this;
            return self.der(i, x) * self.der(j, x);
          };
          m(i, j) = this->compute_integral(f,quadRule);
          if(i != j)
            m(j, i) = m(i, j);
        }
    }
  return mu * m;
}

apsc::Fem1d::FemElement1D::LocalMatrix
apsc::Fem1d::FemElement1D::stiff(const Function &mu,QuadRule const & quadRule) const
{
  using namespace apsc::Fem1d;
  FemElement1D::LocalMatrix m{nodes_.size(), nodes_.size()};
  for(std::size_t i = 0u; i < nodes_.size(); ++i)
    {
      for(std::size_t j = 0u; j <= i; ++j)
        {
          FemElement1D::Function f = [i, j, this, &mu](const double &x) {
            auto const &self = *this;
            return mu(x) * self.der(i, x) * self.der(j, x);
          };
          m(i, j) = this->compute_integral(f,quadRule);
          if(i != j)
            m(j, i) = m(i, j);
        }
    }
  return m;
}

apsc::Fem1d::FemElement1D::LocalMatrix
apsc::Fem1d::FemElement1D::react(const double &c,QuadRule const & quadRule) const
{
  FemElement1D::LocalMatrix m = mass(quadRule);
  return c * m;
}

apsc::Fem1d::FemElement1D::LocalMatrix
apsc::Fem1d::FemElement1D::react(const Function &c,QuadRule const & quadRule) const
{
  using namespace apsc::Fem1d;
  FemElement1D::LocalMatrix m{nodes_.size(), nodes_.size()};
  for(std::size_t i = 0u; i < nodes_.size(); ++i)
    {
      for(std::size_t j = 0u; j <= i; ++j)
        {
          FemElement1D::Function f = [i, j, this, &c](const double &x) {
            auto const &self = *this;
            return c(x) * self(i, x) * self(j, x);
          };
          m(i, j) = this->compute_integral(f,quadRule);
          if(i != j)
            m(j, i) = m(i, j);
        }
    }
  return m;
}

apsc::Fem1d::FemElement1D::LocalVector
apsc::Fem1d::FemElement1D::source(const double &s,QuadRule const & quadRule) const
{
  using namespace apsc::Fem1d;
  FemElement1D::LocalVector v{nodes_.size()};
  for(std::size_t i = 0u; i < nodes_.size(); ++i)
    {
      FemElement1D::Function f = [i, this](const double &x) {
        auto const &self = *this;
        return self(i, x);
      };
      v[i] = this->compute_integral(f,quadRule);
    }
  return s * v;
}

apsc::Fem1d::FemElement1D::LocalVector
apsc::Fem1d::FemElement1D::source(const Function &s, QuadRule const & quadRule) const
{
  using namespace apsc::Fem1d;
  FemElement1D::LocalVector v{nodes_.size()};
  for(std::size_t i = 0u; i < nodes_.size(); ++i)
    {
      FemElement1D::Function f = [i, this, &s](const double &x) {
        auto const &self = *this;
        return s(x) * self(i, x);
      };
      v[i] = this->compute_integral(f,quadRule);
    }
  return v;
}
/*
 * FemElement.cpp
 *
 *  Created on: Aug 30, 2022
 *      Author: forma
 */
