/**
 * @file Material.cpp
 * Implementation of the Material class.
 *
 * @author Haohang Huang
 * @date May 21, 2019
 */

#include "Material.h"

Material::Material(const bool & Anisotropy, const bool & Nonlinearity)
  : anisotropy(Anisotropy), nonlinearity(Nonlinearity), E_(MatrixXd::Zero(6,6)), thermalStrain_(VectorXd::Zero(6))
{
}

Material::~Material()
{
}

const double & Material::modulus() const
{
    return M_;
}

const double & Material::modulusX() const
{
    return Mx_;
}

const double & Material::modulusY() const
{
    return My_;
}

const double & Material::modulusZ() const
{
    return Mz_;
}

const double & Material::modulusG() const
{
    return G_;
}

const MatrixXd & Material::EMatrix() const
{
    return E_;
}

MatrixXd Material::EMatrix(const VectorXd & modulus) const
{
    (void)modulus; // silence warning
    return MatrixXd::Zero(6,6); // to silent warning
}

VectorXd Material::stressDependentModulus(const VectorXd & stress) const
{
    (void)stress; // silence warning
    return VectorXd::Zero(4); // to silent warning
}

const Vector3d & Material::bodyForce() const
{
    return bodyForce_;
}

void Material::setBodyForce(const Vector3d & force)
{
    bodyForce_ = force;
}

const VectorXd & Material::thermalStrain() const
{
    return thermalStrain_;
}

void Material::setThermalStrain(const VectorXd & thermalStrain)
{
    thermalStrain_ = thermalStrain;
}
