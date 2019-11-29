/**
 * @file Nonlinear.h
 * Derived class from Analysis for nonlinear elastic problems.
 *
 * @author Haohang Huang
 * @date May 21, 2019
 */

#ifndef Nonlinear_h
#define Nonlinear_h

#include "Analysis.h"

/* Derived class for solving nonlinear elastic problems.
 */
class Nonlinear : public Analysis
{
  public:
    /* See the documentation of base class Analysis.
     */
    Nonlinear(Mesh & meshInfo);  // ctor cannot be inherited, should explicitly call base class's ctor in derived class's ctor
    ~Nonlinear();
    void solve();

    /**
     * Compute principal stresses at Gaussian points, and update the modulus and E matrix for next iterations.
     *
     * @return A boolean value incidating the convergence status at this iteration.
     */
    bool nonlinearIteration(double damping);

    /**
     * Helper function for computing the principal stresses at current stress state.
     *
     * @param stress [sigma_x, sigma_y, sigma_z, tau_xy, tau_yz, tau_zx]
     * @return The principal stresses in sigma3, sigma2, sigma1 order.
     */
    VectorXd principalStress(const VectorXd & stress) const;

  private:
    int gravityIncrementNum; /* No. of body load (gravity & temperature & residual) increments */
    int loadIncrementNum; /* No. of traffic load (point & edge & face) increments */
    double gravityDamping; /* Damping ratio lambda for body force incremental loading */
    double loadDamping; /* Damping ratio lambda for traffic incremental loading */

};

#endif /* Nonlinear_h */
