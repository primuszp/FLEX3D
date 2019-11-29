/**
 * @file Linear.cpp
 * Implementation of Linear class.
 *
 * @author Haohang Huang
 * @date May 21, 2019
 */

#include "Linear.h"
#include <chrono>
#include <iostream>

Linear::Linear(Mesh & meshInfo) : Analysis(meshInfo)
{
}

Linear::~Linear()
{
}

void Linear::solve()
{
    applyForce();
    assembleStiffness();

    // Option1: SimplicialLDLT <SparseMatrix<double> > solver;
    // Option2: ConjugateGradient <SparseMatrix<double> > solver;
    // SimplicialLDLT is direct solver: Recommended for very sparse and not too large problems (e.g., 2D Poisson eq.)
    // ConjugateGradient is iterative solver: Recommended for large symmetric problems (e.g., 3D Poisson eq.)
    // Ref: https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html

    // SimplicialLDLT <SparseMatrix<double> > solver;
    // ConjugateGradient <SparseMatrix<double> > solver;
    SimplicialLDLT <SparseMatrix<double> > solver;
    solver.compute(globalStiffness);
    nodalDisp = solver.solve(nodalForce);

    // Compute strain and stress and accumulate at each node
    computeStrainAndStress();

    // Average strain and stress at each node
    averageStrainAndStress();
}
