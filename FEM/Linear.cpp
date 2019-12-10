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
	auto start = std::chrono::high_resolution_clock::now();
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start; // in second

    applyForce();

	finish = std::chrono::high_resolution_clock::now();
	elapsed = finish - start; // in second
	// std::cout << "apply: " << elapsed.count() << " s" << std::endl;
	start = std::chrono::high_resolution_clock::now();

    assembleStiffness();

	finish = std::chrono::high_resolution_clock::now();
	elapsed = finish - start; // in second
	// std::cout << "assemble: " << elapsed.count() << " s" << std::endl;
	start = std::chrono::high_resolution_clock::now();

    // Option1: SimplicialLDLT <SparseMatrix<double> > solver;
	// SimplicialLDLT is a direct solver: Recommended for very sparse and not too large problems (e.g., 2D Poisson eq.). We used this in FEM2D.
    // Option2: ConjugateGradient <SparseMatrix<double> > solver;
    // ConjugateGradient is an iterative solver: Recommended for large symmetric problems (e.g., 3D Poisson eq.). Here we have 3D problem and the stiffness matrix is symmetric. CG is very good.
	// Option3: BiCGSTAB <SparseMatrix<double> > solver;
	// Bi-CG, a little slower than classic CG.
	// After test, for large 3D problem, LDLT is extremely slow. Two CG methods ~100s for a 150000 node system. Other methods either slow or wrong answer.
    // Ref:
	// https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html
	// https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html

	// SimplicialLDLT <SparseMatrix<double> > solver;
    ConjugateGradient <SparseMatrix<double> > solver;
	// BiCGSTAB <SparseMatrix<double> > solver;
    solver.compute(globalStiffness);
	finish = std::chrono::high_resolution_clock::now();
	elapsed = finish - start; // in second
	// std::cout << "compute: " << elapsed.count() << " s" << std::endl;
	start = std::chrono::high_resolution_clock::now();

    nodalDisp = solver.solve(nodalForce);

	finish = std::chrono::high_resolution_clock::now();
	elapsed = finish - start; // in second
	// std::cout << "solve: " << elapsed.count() << " s" << std::endl;
	start = std::chrono::high_resolution_clock::now();

    // Compute strain and stress and accumulate at each node
    computeStrainAndStress();

	finish = std::chrono::high_resolution_clock::now();
	elapsed = finish - start; // in second
	// std::cout << "compute stress: " << elapsed.count() << " s" << std::endl;
	start = std::chrono::high_resolution_clock::now();

    // Average strain and stress at each node
    averageStrainAndStress();

	finish = std::chrono::high_resolution_clock::now();
	elapsed = finish - start; // in second
	// std::cout << "avg stress: " << elapsed.count() << " s" << std::endl;
	start = std::chrono::high_resolution_clock::now();
}
