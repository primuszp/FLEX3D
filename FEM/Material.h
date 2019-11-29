/**
 * @file Material.h
 * Material class for the element properties.
 *
 * @author Haohang Huang
 * @date May 21, 2019
 */

#ifndef Material_h
#define Material_h

#include "Eigen/Eigen"
#include <vector>

using namespace Eigen;

/** Material class for storing the engineering properties of elements. */
class Material
{
  public:

    /**
     * Custom constructor to create an element material.
     *
     * @param anisotropy A sign for material isotropy. 0 if isotropic, 1 if cross-anisotropic.
     * @param nonlinearity A sign for material linearity. 0 if linear elastic, 1 if nonlinear elastic.
     */
    Material(const bool & anisotropy, const bool & nonlinearity);

    /**
     * Destructor.
     *
     * @note As an abstract class, the destructor must be virtual.
     */
    virtual ~Material();

    /**
     * Get the modulus of the element.
     *
     * @return The modulus value.
     */
    const double & modulus() const;
    const double & modulusX() const;
    const double & modulusY() const;
    const double & modulusZ() const;
    const double & modulusG() const;

    /**
     * Get the stress-strain constitutive matrix of the element. Linear elastic
     * material will use this to get the constant E matrix.
     *
     * @return The 6-by-6 E matrix.
     */
    const MatrixXd & EMatrix() const;

    /**
     * Compute the stress-dependent E matrix from resilient modulus. Nonlinear
     * elastic material will use this to get the stress-dependent E matrix.
     *
     * @param modulus The stress-dependent resilient modulus. If isotropic, it's
     * just a single value (but wrapped in a VectorXd); if anisotropic, it's a 4-by-1 vector for
     * X, Y, Z modulus & shear modulus.
     * @return The constitutive matrix.
     */
    virtual MatrixXd EMatrix(const VectorXd & modulus) const;

    /**
     * Compute the stress-dependent resilient modulus of the element. Used in nonlinear scheme.
     *
     * @param stress The principal stresses in sigma3, sigma2, sigma1 order.
     * @return The stress-dependent resilient modulus computed from models. If
     * isotropic, it's just a single value (but wrapped in VectorXd); if anisotropic, it's a 4-by-1 vector
     * for X, Y, Z modulus & shear modulus.
     */
    virtual VectorXd stressDependentModulus(const VectorXd & stress) const;

    /**
     * Get the body force to be used in the load condition.
     *
     * @return The body force as a 3-by-1 vector for 3D problem.
     */
    const Vector3d & bodyForce() const;

    /**
     * Assign the body force to allow incremental loading in nonlinear scheme.
     *
     * @param force The incremental body force to be assigned.
     */
    void setBodyForce(const Vector3d & force);

    /**
     * Get the thermal strain to be used in the stress computation.
     *
     * @return The thermal strain as a 6-by-1 vector for 3D problem.
     */
    const VectorXd & thermalStrain() const;

    /**
     * Assign the thermal strain to allow incremental loading in nonlinear scheme.
     *
     * @param thermalStrain The incremental thermal strain to be assigned.
     */
    void setThermalStrain(const VectorXd & thermalStrain);

    /** A sign for material isotropy. 0 if isotropic, 1 if cross-anisotropic. */
    bool anisotropy;

    /** A sign for material linearity. 0 if linear elastic, 1 if nonlinear elastic. */
    bool nonlinearity;

  protected:

    /** The Young's/Resilient modulus. Will be assigned in derived class for both linear and nonlinear (the initial guess Modulus) materials */
    double M_, Mx_, My_, Mz_, G_; // Isotropic modulus, X modulus, Y modulus, Z modulus, shear modulus

    /** Poisson ratio */
    double v_, vxy_, vyz_, vzx_; // Isotropic Poission, XY Poisson, YZ Poisson, ZX Poisson

    /** The 6-by-6 stress-strain constitutive matrix sigma = E * e */
    MatrixXd E_;

    /** The body force (unit weight) */
    Vector3d bodyForce_;

    /** The thermal strain */
    VectorXd thermalStrain_;

};

#endif /* Material_h */
