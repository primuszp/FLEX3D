/**
 * @file LinearElastic.cpp
 * Implementation of LinearElastic class.
 *
 * @author Haohang Huang
 * @date May 21, 2019
 */

#include "LinearElastic.h"

LinearElastic::LinearElastic(const bool & anisotropy, const bool & nonlinearity, const std::vector<double> & properties)
  : Material(anisotropy, nonlinearity)
{
    // Elasticity matrix: Taylor's book Chapter 2
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

LinearElastic::~LinearElastic()
{
}
