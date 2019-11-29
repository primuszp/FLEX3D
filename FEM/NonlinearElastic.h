/**
 * @file NonlinearElastic.h
 * Derived class from Material for nonlinear elastic materials.
 *
 * @author Haohang Huang
 * @date May 21, 2019
 */

#ifndef NonlinearElastic_h
#define NonlinearElastic_h

#include "Material.h"

class NonlinearElastic : public Material
{
  public:
    /** See the documentation of base class Material. */
    NonlinearElastic(const bool & anisotropy, const bool & nonlinearity, const std::vector<double> & properties, const int & model, const std::vector<double> & parameters);
    ~NonlinearElastic();

    VectorXd stressDependentModulus(const VectorXd & stress) const;
    MatrixXd EMatrix(const VectorXd & modulus) const;

  protected:
    int modelNo; /* Designator for resilient model used */
    std::vector<double> coeff; /** Regression coefficients used in the resilient model */

};

#endif /* NonlinearElastic_h */
