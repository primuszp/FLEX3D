/**
 * @file Nonlinear.cpp
 * Implementation of Nonlinear class.
 *
 * @author Haohang Huang
 * @date May 21, 2019
 */

#include "Nonlinear.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <functional>

Nonlinear::Nonlinear(Mesh & meshInfo) : Analysis(meshInfo)
{
    gravityIncrementNum = (mesh.iterations)[0];
    loadIncrementNum = (mesh.iterations)[1];
    gravityDamping = (mesh.iterations)[2];
    loadDamping = (mesh.iterations)[3];
}

Nonlinear::~Nonlinear()
{
}

void Nonlinear::solve()
{
// Several options that can affect the results:
// 1. Whether to do one last run after the convergence
// 2. Change damping ratio
// 3. Whether to do incremental loading, and how many increments
// 4. A strict (check every individual Gaussian) or loose (check only at center) convergence criteria, also to use different modulus for every Gaussian or for every element
// 5. Convergence criteria to compare with old modulus or new modulus
// no increment, all Gaussian but only check convergence at center, 0.3 damping, no final run --> 2.57e-4 result

    bool incremental = true;
    if (gravityIncrementNum == 0 && loadIncrementNum == 0) incremental = false; // input two zeros means I don't want to have incremental loading

    if (incremental) {
        // -----------------------------------------------------------------------------
        // --------------- Start of Incremental Loading Scheme -------------------------
        // -----------------------------------------------------------------------------

        // Gravity, temperature, and residual stress increments
        // Idea: for each material, re-assign the body force, residual stress and
        // thermal strain incrementally. At the beginning we should calculate the increment value
        // A good observation: with gravity load only, the stress is independent with the modulus,
        // so any arbitrary initial guess of the modulus won't affect the stress-dependent modulus.
        const std::vector<Material*> & materials = mesh.materialList;
        std::vector<Vector3d> gravityIncrement;
        std::vector<VectorXd> thermalIncrement;
        gravityIncrement.reserve(materials.size());
        thermalIncrement.reserve(materials.size());

        // Compute the gravity/thermal/residual load increment for each material
        for (auto & m : materials) {
            gravityIncrement.push_back(m->bodyForce() / gravityIncrementNum);
            thermalIncrement.push_back(m->thermalStrain() / gravityIncrementNum);
        }

        // Workflow: apply ith increment --> achieve convergence --> apply (i+1)th increment, repeat
        for (int ic = 1; ic <= gravityIncrementNum; ic++) {
            // Apply the load incrementally for each material
            for (unsigned m = 0; m < materials.size(); m++) {
                materials[m]->setBodyForce(gravityIncrement[m] * ic);
                materials[m]->setThermalStrain(thermalIncrement[m] * ic);
            }

            // Achieve the modulus convergence at each increment
            bool nonlinearConvergence = false;
            int count = 0;
            while (!nonlinearConvergence) { // convergence criteria
                // Assemble the K and F based on the mesh information (without applying any
                // load, this is for body force and temperature load incremental only!)
                // @BUG(solved) Normally applyForce() will initialize the global force vector
                // and the assembleStiffness() function below will always do += for nodalForce. But
                // in the body force increment stage, applyForce() should not be called, therefore we
                // need to manually reset the nodalForce otherwise it will keeps accumulating.
                nodalForce = VectorXd::Zero(3 * mesh.nodeCount());
                assembleStiffness();

                // Solve K U = F
                SimplicialLDLT <SparseMatrix<double> > solver;
                solver.compute(globalStiffness);
                nodalDisp = solver.solve(nodalForce);

                // Traverse each element, compute stress at Gaussian points, and update the modulus for the next (i + 1) iteration (if current iteration is i)
                nonlinearConvergence = nonlinearIteration(gravityDamping);

                count++;
            }
            std::cout << "Body Force Increment No." << ic << ", No. of iterations = " << count << std::endl;
            std::cout << "-----------------------------------------" << std::endl;

            // For the exit iteration, the new converged modulus is updated, but the nodalDisp
            // is for the last iteration, so we should do one more solve to match the modulus & displacment
            nodalForce = VectorXd::Zero(3 * mesh.nodeCount());
            assembleStiffness();
            SimplicialLDLT <SparseMatrix<double> > solver;
            solver.compute(globalStiffness);
            nodalDisp = solver.solve(nodalForce);
        }

        // Traffic load increments (point/edge/face load)
        std::vector<double> pointLoadIncrement = mesh.loadValue;
        std::vector<std::vector<double> > edgeLoadIncrement = mesh.edgeLoadValue;
        std::vector<std::vector<double> > faceLoadIncrement = mesh.faceLoadValue;

        // Compute the point/edge/face load increment
        // std::transform is the way to do element-wise std::vector operation
        std::transform(pointLoadIncrement.begin(), pointLoadIncrement.end(), pointLoadIncrement.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, 1.0 / loadIncrementNum)); // in-place change
        for (auto & e : edgeLoadIncrement)
            std::transform(e.begin(), e.end(), e.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, 1.0 / loadIncrementNum));
        for (auto & e : faceLoadIncrement)
            std::transform(e.begin(), e.end(), e.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, 1.0 / loadIncrementNum));

        // Workflow: apply ith increment --> achieve convergence --> apply (i+1)th increment, repeat
        for (int ic = 1; ic <= loadIncrementNum; ic++) {
            // Apply the load incrementally
            std::transform(pointLoadIncrement.begin(), pointLoadIncrement.end(), mesh.loadValue.begin(), std::bind(std::multiplies<int>(), std::placeholders::_1, ic));
            for (unsigned e = 0; e < edgeLoadIncrement.size(); e++)
                std::transform(edgeLoadIncrement[e].begin(), edgeLoadIncrement[e].end(), (mesh.edgeLoadValue)[e].begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, ic));
            for (unsigned e = 0; e < faceLoadIncrement.size(); e++)
                std::transform(faceLoadIncrement[e].begin(), faceLoadIncrement[e].end(), (mesh.faceLoadValue)[e].begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, ic));

            // Achieve the modulus and tension convergence at each increment
            bool nonlinearConvergence = false;
            int count = 0;
            while (!nonlinearConvergence) { // convergence criteria
                // Assemble the K and F based on the mesh information (with traffic load applied)
                applyForce();
                assembleStiffness();

                // Solve K U = F
                SimplicialLDLT <SparseMatrix<double> > solver;
                solver.compute(globalStiffness);
                nodalDisp = solver.solve(nodalForce);

                // Traverse each element, compute stress at Gaussian points, and update the modulus for the next (i + 1) iteration (if current iteration is i)
                nonlinearConvergence = nonlinearIteration(loadDamping);

                count++;
            }
            std::cout << "Traffic Load Increment No." << ic << ", No. of iterations = " << count << std::endl;
            std::cout << "-----------------------------------------" << std::endl;

            // For the exit iteration, the new converged modulus is updated, but the nodalDisp
            // is for the last iteration, so we should do one more solve to match the modulus & displacment
            applyForce();
            assembleStiffness();
            SimplicialLDLT <SparseMatrix<double> > solver;
            solver.compute(globalStiffness);
            nodalDisp = solver.solve(nodalForce);

        }
        // -----------------------------------------------------------------------------
        // ----------------- End of Incremental Loading Scheme -------------------------
        // -----------------------------------------------------------------------------
    }
    else { // non-incremental loading
        // -------------------------------------------------------------------------
        // --------------- Start of Direct Loading Scheme --------------------------
        // -------------------------------------------------------------------------
        bool nonlinearConvergence = false;
        int i = 0; // for debug print only
        while (!nonlinearConvergence) { // convergence criteria
            // Assemble the K and F based on the mesh information. At 1st iteration, the initial guess modulus M0 will be used; later on at ith iteration, the stress-dependent modulus updated from (i - 1) iteratiion will be used
            applyForce();
            assembleStiffness();

            // Solve K U = F
            SimplicialLDLT <SparseMatrix<double> > solver;
            solver.compute(globalStiffness);
            nodalDisp = solver.solve(nodalForce);

            // Traverse each element, compute stress at Gaussian points, and update the modulus for the next (i + 1) iteration (if current iteration is i)
            nonlinearConvergence = nonlinearIteration(0.3);

            i++;
        }
        std::cout << "No. of Nonlinear Iterations: " << i << std::endl;

        applyForce();
        assembleStiffness();
        SimplicialLDLT <SparseMatrix<double> > solver;
        solver.compute(globalStiffness);
        nodalDisp = solver.solve(nodalForce);
        // For the exit iteration, the new converged modulus is updated, but the nodalDisp
        // is for the last iteration, so we should do one more solve to match the modulus & displacment
        // -------------------------------------------------------------------------
        // -------------------- End of Direct Loading Scheme -----------------------
        // -------------------------------------------------------------------------
    }

    // After the nonlinear scheme converge, compute nodal strain & stress from the final displacment results
    computeStrainAndStress();
    averageStrainAndStress();

}

bool Nonlinear::nonlinearIteration(double damping)
{
    bool convergence = true;
    double sumError = 0;
    double sumModulus = 0;
    double criteria1 = 0.05;
    double criteria2 = 0.002;

    Element* curr;
    int numNodes; // number of nodes belong to the element
    int numGaussianPt; // number of Gaussian points of the element
    for (int i = 0; i < mesh.elementCount(); i++) {
        curr = mesh.elementArray()[i];
        Material* material = curr->material();
        if (material->nonlinearity) { // compute stress for nonlinear elastic element only, skip all linear elastic ones
            const VectorXi & nodeList = curr->getNodeList();
            numNodes = curr->getSize();
            numGaussianPt = (int)curr->shape()->gaussianPt().size();

            // Assemble the nodal displacement vector for this element
            VectorXd nodeDisp(3 * numNodes);
            for (int j = 0; j < numNodes; j++) {
                nodeDisp(3 * j) = nodalDisp(3 * nodeList(j));
                nodeDisp(3 * j + 1) = nodalDisp(3 * nodeList(j) + 1);
                nodeDisp(3 * j + 2) = nodalDisp(3 * nodeList(j) + 2);
            }

            // Step 1: Compute stress at gaussian points based on cached MR & E from last iteration
            // Step 2: Update new modulus based on the stress from step 1 and mix with old modulus via damping ratio
            // Step 3: Cache the modulus to be used in the next iteration
            // Note: EMatrix(VectorXd), for isotropic case, modulus should be wrapped into a VectorXd variable; for anisotropic case, it is already a VectorXd
            // Note: Tutu's approach only use the center Gaussian point for the whole element, as follows
            // MatrixXd B = curr->BMatrix(curr->shape()->gaussianPt(4));
            // VectorXd strain = B * nodeDisp; // e = B * u
            // double modulus_old = (curr->modulusAtGaussPt)(4); // M_(i-1)
            // VectorXd modulus_old_vec << modulus_old;
            // VectorXd stress = material->EMatrix(modulus_old_vec) * (strain - curr->thermalStrain()); // sigma = E_(i-1) * (e - e0), note that the M and E are both from previous iteration
            // VectorXd modulus_vec = material->stressDependentModulus(principalStress(stress));
            // double modulus_new = modulus_vec(1); // M_i, 1 for vertical modulus
            // double modulus = (1 - damping) * modulus_old + damping * modulus_new; // true M_i after applying damping ratio

            // For isotropy case (or a simplified anisotropy case), iterate only on the single modulus (vertical modulus for anisotropy)
            if (!material->anisotropy) {
                for (int g = 0; g < numGaussianPt; g++) {
                    // Step 1: Compute stress
                    MatrixXd B = curr->BMatrix(curr->shape()->gaussianPt(g));
                    VectorXd strain = B * nodeDisp; // e = B * u
                    double modulus_old = (curr->modulusAtGaussPt)(g); // M_(i-1)
                    VectorXd modulus_old_vec(1);
                    modulus_old_vec << modulus_old;
                    VectorXd stress = material->EMatrix(modulus_old_vec) * (strain - curr->thermalStrain()); // sigma = E_(i-1) * (e - e0), note that the M and E are both from previous iteration
                    // Step 2: Compute new modulus and mix by damping ratio
                    VectorXd modulus_vec = material->stressDependentModulus(principalStress(stress));
                    double modulus_new = modulus_vec(2); // [M_X, M_Y, M_Z, G], 2 for vertical (Z) modulus
                    double modulus = (1 - damping) * modulus_old + damping * modulus_new; // true M_i after applying damping ratio
                    // Step 3: Cache new modulus at Gaussian points
                    (curr->modulusAtGaussPt)(g) = modulus; // for isotropic, new modulus is a double

                    // Convergence criteria
                    // Criteria 1: modulus error within 5%
                    // Criteraia 2: accumulative modulus error within 0.2%
                    double error = std::abs(modulus - modulus_new);
                    bool strict = false;
                    if (strict) {
                        // Strict version (check criteria 1 & 2 at all Gaussian points)
                        // Criteria 1
                        if (error / modulus_old > criteria1) // tutu uses modulus_old, but I want to use modulus
                            convergence = false;
                            // should I just return false here? No! Because you still need to update all elements' modulus synchronously
                        // Criteria 2
                        sumError += error * error;
                        sumModulus += modulus_old * modulus_old; // tutu uses modulus_old, but I want to use modulus
                    } else {
                        // Loose version (check criteria 1 & 2 at center Gaussian only)
                        int centerGaussian = (numGaussianPt - 1) / 2;
                        // Criteria 1
                        if (g == centerGaussian && error / modulus_old > criteria1) // tutu uses modulus_old, but I want to use modulus
                            convergence = false;
                            // should I just return false here? No! Because you still need to update all elements' modulus synchronously
                        // Criteraia 2
                        if (g == centerGaussian) {
                            sumError += error * error;
                            sumModulus += modulus_old * modulus_old; // tutu uses modulus_old, but I want to use modulus
                        }
                    }

                }
            }
            // For anisotropy case, iterate on all 4 moduli (X, Y, Z, G)
            else {
                for (int g = 0; g < numGaussianPt; g++) {
                    // Step 1: Compute stress
                    MatrixXd B = curr->BMatrix(curr->shape()->gaussianPt(g));
                    VectorXd strain = B * nodeDisp; // e = B * u
                    VectorXd modulus_old = (curr->modulusAtGaussPt).row(g); // M_(i-1)
                    VectorXd stress = material->EMatrix(modulus_old) * (strain - curr->thermalStrain()); // sigma = E_(i-1) * (e - e0), note that the M and E are both from previous iteration
                    // Step 2: Compute new modulus and mix by damping ratio
                    VectorXd modulus_new = material->stressDependentModulus(principalStress(stress)); // M_i
                    VectorXd modulus = (1 - damping) * modulus_old + damping * modulus_new; // true M_i after applying damping ratio
                    // Step 3: Cache new modulus at Gaussian points
                    (curr->modulusAtGaussPt).row(g) = modulus; // for isotropic, new modulus is a VectorXd

                    // Convergence criteria
                    // Criteria 1: modulus error within 5%
                    // Criteraia 2: accumulative modulus error within 0.2%
                    bool strict = false;
                    VectorXd error = (modulus - modulus_new).array() / modulus_old.array(); // tutu uses modulus_old, but I want to use modulus
                    error = error.array().abs(); // or error.cwiseAbs(), https://stackoverflow.com/questions/25340940/how-do-i-compute-the-absolute-value-of-a-vector-in-eigen
                    if (strict) {
                        // Strict version (check criteria 1 & 2 at all Gaussian points)
                        // Criteria 1
                        if (error(0) > criteria1 && error(1) > criteria1 && error(2) > criteria1 && error(3) > criteria1)
                            convergence = false;
                        // Criteria 2
                        error = error.array().square();
                        modulus_old = modulus_old.array().square();
                        sumError += error.sum();
                        sumModulus += modulus_old.sum(); // tutu uses modulus_old, but I want to use modulus
                    } else {
                        // Loose version (check criteria 1 & 2 at center Gaussian only)
                        int centerGaussian = (numGaussianPt - 1) / 2;
                        // Criteria 1
                        if (g == centerGaussian && error(0) > criteria1 && error(1) > criteria1 && error(2) > criteria1 && error(3) > criteria1)
                            convergence = false;
                        // Criteraia 2
                        error = error.array().square();
                        modulus_old = modulus_old.array().square();
                        if (g == centerGaussian) {
                            sumError += error.sum();
                            sumModulus += modulus_old.sum();
                        }
                    }

                }
            }

        }

    }

    std::cout << "Sum of Error Percentage (%): " << sumError / sumModulus * 100 << std::endl;

    return (sumError / sumModulus < criteria2 && convergence) ? true : false;

}

VectorXd Nonlinear::principalStress(const VectorXd & stress) const
{
    // http://academic.uprm.edu/pcaceres/Courses/MMII/IMoM-5A.pdf
    // VectorXd stress = [sigma_x, sigma_y, sigma_z, tau_xy, tau_yz, tau_zx]
    // stress tensor =
    // [sigma_x, tau_yx, tau_zx]
    // [tau_xy, sigma_y, tau_zy]
    // [tau_xz, tau_yz, sigma_z]
    // Eigenvalues of stress tensor is the principal stress sigma3, sigma2, sigma1 (Eigen returns in increasing order)
    // In our coordinates, -:compression +:tension
    MatrixXd tensor(3,3);
    tensor << stress(0), stress(3), stress(5),
              stress(3), stress(1), stress(4),
              stress(5), stress(4), stress(2);
    SelfAdjointEigenSolver<MatrixXd> es(tensor, EigenvaluesOnly);
    return es.eigenvalues();

}
