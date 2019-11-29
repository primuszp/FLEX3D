/**
 * @file NonlinearElastic.cpp
 * Implementation of NonlinearElastic class.
 *
 * @author Haohang Huang
 * @date May 21, 2019
 */

#include "NonlinearElastic.h"
#include <cmath>
#include <iostream>
NonlinearElastic::NonlinearElastic(const bool & anisotropy, const bool & nonlinearity, const std::vector<double> & properties, const int & model, const std::vector<double> & parameters)
  : Material(anisotropy, nonlinearity), modelNo(model), coeff(parameters)
{
    // Just copy from LinearElastic constructor
    // For linear elastic problem, just store a constant E matrix is sufficient, so initialize it in constructor
    int i = 0;
    if (!anisotropy) {
        // Isotropic: Modulus, Poisson's ratio, body force (x,y,z), thermal coefficient, temperature change
        double M = properties[i++]; M_ = M; // initialize members of base class Material
        double v = properties[i++]; v_ = v;

        // Compute the stress-strain constitutive matrix
        E_ << 1 - v, v, v, 0, 0, 0,
              v,   1-v, v, 0, 0, 0,
              v,  v,  1-v, 0, 0, 0,
              0, 0, 0, (1-2*v)/2, 0, 0,
              0, 0, 0, 0, (1-2*v)/2, 0,
              0, 0, 0, 0, 0, (1-2*v)/2;
        E_ = E_ * M / (1+v) /(1-2*v);
    }
    else {
        // Orthotropic: https://knowledge.autodesk.com/support/moldflow-insight/learn-explore/caas/CloudHelp/cloudhelp/2018/ENU/MoldflowInsight-Analyses/files/GUID-29FDEC3E-0D52-4712-A98C-540228A4C33B-htm.html
        // Anisotropic: Modulus X, Modulus Y, Modulus Z, Poisson's ratio XY, Poisson's ratio YZ, Poisson's ratio ZX, Shear Modulus G, body force (x,y,z), thermal coefficient, temperature change
        double Mx = properties[i++]; Mx_ = Mx;
        double My = properties[i++]; My_ = My;
        double Mz = properties[i++]; Mz_ = Mz;
        double vxy = properties[i++]; vxy_ = vxy;
        double vyz = properties[i++]; vyz_ = vyz;
        double vzx = properties[i++]; vzx_ = vzx;
        // note the relation vxy/My = vyx/Mx
        double vyx = Mx/My * vxy;
        double vzy = My/Mz * vyz;
        double vxz = Mz/Mx * vzx;
        double G = properties[i++]; G_ = G; // here I simplify a little bit by set G_xy = G_yz = G_zx = G

        // Helper coefficient
        double A = (1 - vxy*vyx - vyz*vzy - vxz*vzx - 2*vxy*vyz*vzx) / (Mx*My*Mz);

        // Compute the stress-strain constitutive matrix
        E_ << (1 - vyz*vzy) / (My*Mz*A), (vyz + vzx*vyz) / (My*Mz*A), (vzx + vyz*vzy) / (My*Mz*A), 0, 0, 0,
              (vxy + vzx*vxz) / (Mz*Mx*A), (1 - vzx*vxz) / (Mz*Mx*A), (vzy + vzx*vxy) / (Mz*Mx*A), 0, 0, 0,
              (vxz + vxy*vyz) / (Mx*My*A), (vzy + vxz*vyz) / (Mx*My*A), (1 - vxy*vyx) / (Mx*My*A), 0, 0, 0,
              0, 0, 0, G, 0, 0,
              0, 0, 0, 0, G, 0,
              0, 0, 0, 0, 0, G;
    }

    // Assign body force (unit weight)
    bodyForce_ << properties[i], properties[i+1], properties[i+2];
    i += 3;

    // Assign thermal parameters
    double alpha  = properties[i++];
    double deltaT = properties[i++]; /** The temperature change (assume same across the entire element) */
    double strain = alpha * deltaT;
    thermalStrain_ << strain, strain, strain, 0, 0, 0;
}

NonlinearElastic::~NonlinearElastic()
{
}

VectorXd NonlinearElastic::stressDependentModulus(const VectorXd & stress) const
{
    // stress(0)-sigma3; stress(1)-sigma2; stress(2)-sigma1
    // Bulk stress: theta = sigma1 + sigma2 + sigma3
    // Deviatoric stress: sigma_d = sigma1 - sigma3
    // Octahedral shear stress (usually meaningful in 3D): tau_oct = sqrt((1-2)^2 + (2-3)^2 + (1-3)^2) / 3, for 2D := sqrt(2)/3 * sigma_d
    double bulk = std::abs( stress(2) + stress(1) + stress(0) );
    double deviator = std::abs( stress(2) - stress(0) );
    double octahedral = std::sqrt(std::pow(stress(2) - stress(1), 2) + std::pow(stress(1) - stress(0), 2) + std::pow(stress(2) - stress(0), 2)) / 3;
    double atm = 14.696; // atmospheric pressure, 101.325 kPa or 14.7 psi

    // Calculated new stress-dependent resilient modulus from different models
    double Mx = 0, My = 0, Mz = 0, G = 0; // X modulus, Y modulus, Z modulus, shear modulus
    switch (modelNo) {
        case 1 :
            // Hicks and Monismith (1971) K-theta model
            // M = k1 * theta^k2
            if (!anisotropy) { // if isotropy, just take the vertical (Z) modulus
                Mz = coeff[0] * std::pow(bulk, coeff[1]);
            } else {
                Mx = coeff[0] * std::pow(bulk, coeff[1]); // although K-theta only has 2 parameters, we fill 0 to complement the triplet
                My = coeff[3] * std::pow(bulk, coeff[4]);
                Mz = coeff[6] * std::pow(bulk, coeff[7]);
                G = coeff[9] * std::pow(bulk, coeff[10]);
            }
            break;
        case 2 :
            // Uzan (1985) model:
            // M = k1 * theta^k2 * sigma_d^k3
            if (!anisotropy) {
                Mz = coeff[0] * std::pow(bulk, coeff[1]) * std::pow(deviator, coeff[2]);
            } else {
                Mx = coeff[0] * std::pow(bulk, coeff[1]) * std::pow(deviator, coeff[2]);
                My = coeff[3] * std::pow(bulk, coeff[4]) * std::pow(deviator, coeff[5]);
                Mz = coeff[6] * std::pow(bulk, coeff[7]) * std::pow(deviator, coeff[8]);
                G = coeff[9] * std::pow(bulk, coeff[10]) * std::pow(deviator, coeff[11]);
            }
            break;
        case 3 :
            // Pezo (1993) UT-Austin model
            // M = k1 * sigma_d^k2 * sigma_3^k3
            if (!anisotropy) {
                Mz = coeff[0] * std::pow(deviator, coeff[1]) * std::pow(std::abs(stress(0)), coeff[2]);
            } else {
                Mx = coeff[0] * std::pow(deviator, coeff[1]) * std::pow(std::abs(stress(0)), coeff[2]);
                My = coeff[3] * std::pow(deviator, coeff[4]) * std::pow(std::abs(stress(0)), coeff[5]);
                Mz = coeff[6] * std::pow(deviator, coeff[7]) * std::pow(std::abs(stress(0)), coeff[8]);
                G = coeff[9] * std::pow(deviator, coeff[10]) * std::pow(std::abs(stress(0)), coeff[11]);
            }
            break;
        case 4 :
            // NCHRP 1-37A (2004) MEPDG model
            // M = k1 * p_a * (theta/p_a)^k2 * (tau_oct/p_a + 1)^k3
            if (!anisotropy) {
                Mz = coeff[0] * atm * std::pow(bulk/atm, coeff[1]) * std::pow(octahedral/atm + 1, coeff[2]);
            } else {
                Mx = coeff[0] * atm * std::pow(bulk/atm, coeff[1]) * std::pow(octahedral/atm + 1, coeff[2]);
                My = coeff[3] * atm * std::pow(bulk/atm, coeff[4]) * std::pow(octahedral/atm + 1, coeff[5]);
                Mz = coeff[6] * atm * std::pow(bulk/atm, coeff[7]) * std::pow(octahedral/atm + 1, coeff[8]);
                G = coeff[9] * atm * std::pow(bulk/atm, coeff[10]) * std::pow(octahedral/atm + 1, coeff[11]);
            }
            break;
    }
    VectorXd result(4);
    result << Mx, My, Mz, G;
    return result;
}

MatrixXd NonlinearElastic::EMatrix(const VectorXd & modulus) const
{
    MatrixXd E(6,6);
    if (!anisotropy) {
        double v = v_;
        E <<  1 - v, v, v, 0, 0, 0,
              v,   1-v, v, 0, 0, 0,
              v,  v,  1-v, 0, 0, 0,
              0, 0, 0, (1-2*v)/2, 0, 0,
              0, 0, 0, 0, (1-2*v)/2, 0,
              0, 0, 0, 0, 0, (1-2*v)/2;
        E = E * modulus(0) / (1+v) /(1-2*v); // for isotropic case, modulus(0) is the constant modulus
    } else {
        double Mx = modulus(0);
        double My = modulus(1);
        double Mz = modulus(2);
        double vxy = vxy_;
        double vyz = vyz_;
        double vzx = vzx_;
        double vyx = Mx/My * vxy;
        double vzy = My/Mz * vyz;
        double vxz = Mz/Mx * vzx;
        double G = modulus(3);

        // Helper coefficient
        double A = (1 - vxy*vyx - vyz*vzy - vxz*vzx - 2*vxy*vyz*vzx) / (Mx*My*Mz);

        // Compute the stress-strain constitutive matrix
        E <<  (1 - vyz*vzy) / (My*Mz*A), (vyz + vzx*vyz) / (My*Mz*A), (vzx + vyz*vzy) / (My*Mz*A), 0, 0, 0,
              (vxy + vzx*vxz) / (Mz*Mx*A), (1 - vzx*vxz) / (Mz*Mx*A), (vzy + vzx*vxy) / (Mz*Mx*A), 0, 0, 0,
              (vxz + vxy*vyz) / (Mx*My*A), (vzy + vxz*vyz) / (Mx*My*A), (1 - vxy*vyx) / (Mx*My*A), 0, 0, 0,
              0, 0, 0, G, 0, 0,
              0, 0, 0, 0, G, 0,
              0, 0, 0, 0, 0, G;
    }
    return E;
}
