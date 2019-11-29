/**
 * @file Analysis.cpp
 * Implementation of Analysis class.
 *
 * @author Haohang Huang
 * @date May 20, 2019
 */

// #define _USE_MATH_DEFINES
// #include <cmath>
#include "Analysis.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <unordered_map>

Analysis::Analysis(Mesh & meshInfo)
  :
  mesh(meshInfo), globalStiffness(3 * mesh.nodeCount(), 3 * mesh.nodeCount()),
  nodalDisp(VectorXd::Zero(3 * mesh.nodeCount())), nodalForce(VectorXd::Zero(3 * mesh.nodeCount())),
  nodalStrain(MatrixXd::Zero(mesh.nodeCount(), 6)), nodalStress(MatrixXd::Zero(mesh.nodeCount(), 6))
{
}

Analysis::~Analysis()
{
}

void Analysis::applyForce()
{
    // @BUG (solved) previous miss this initialization step, therefore in the nonlinear analysis the force will accumulate in every iteration and blow up!!!
    // Other variable such as nodelDisp, globalStiffness, nodalStrain, nodalStress will be rewritten every time, so doesn't matter
    nodalForce = VectorXd::Zero(3 * mesh.nodeCount());

    // Apply point load
    for (unsigned i = 0; i < mesh.loadNodeList.size(); i++)
        nodalForce(mesh.loadNodeList[i]) += mesh.loadValue[i];

    // Traverse all elements with edge and face load
    Element* curr;

    // Apply edge load
    for (unsigned i = 0; i < mesh.edgeLoadElementList.size(); i++) {
        curr = mesh.getElement(mesh.edgeLoadElementList[i]);
        const VectorXi & nodeList = curr->getNodeList();
        int elementType = curr->getSize();
        // The shape function for each edge in an isoparametric element is the same
        const std::vector<double> & gaussianPoint = curr->shape()->edgeGaussianPt(); // length 3 vector
        const std::vector<double> & gaussianWeight = curr->shape()->edgeGaussianWt(); // length 3 vector
        int numGaussianPts = (int)gaussianPoint.size();

        // The edge load information
        const std::vector<int> & edges = mesh.loadEdgeList[i];
        const std::vector<double> & loads = mesh.edgeLoadValue[i];

        VectorXd forceVec(3 * elementType);
        forceVec.setZero();
        // Traverse all loaded edges in this element
        for (unsigned j = 0; j < edges.size(); j++) {

            Vector3d load(loads[3 * j], loads[3 * j + 1], loads[3 * j + 2]); // 3x1 load vector

            // The global coordinates of the edge nodes
            const std::vector<int> & edgeNodes = curr->shape()->edge(edges[j]);
            int numEdgeNodes = (int)edgeNodes.size();
            MatrixXd nodeCoord(3, numEdgeNodes);
            for (int n = 0; n < numEdgeNodes; n++)
                nodeCoord.col(n) = mesh.getNode(nodeList(edgeNodes[n]))->getGlobalCoord();

            // Integration at gaussian points
            // sum: N^T * F * |J| * W(i)
            VectorXd result(3 * numEdgeNodes); result.setZero(); // 9x1 vector
            for (int g = 0; g < numGaussianPts; g++) {
                MatrixXd N = curr->shape()->edgeFunctionMat(g); // 3x9 shape matrix
                VectorXd globalDeriv = nodeCoord * curr->shape()->edgeFunctionDeriv(g); // 3x3 matrix * 3x1 vector = 3x1 vector [dx/d(xi); dy/d(xi); dz/d(xi)]
                double dx = globalDeriv(0);
                double dy = globalDeriv(1);
                double dz = globalDeriv(2);
                double jacobianDet = std::sqrt(dx * dx + dy * dy + dz * dz);
                result += N.transpose() * load * jacobianDet * gaussianWeight[g]; // 9x1 vector, (Fx,Fy,Fz) at 3 nodes, so 3 * 3 = 9
            }

            // Assemble edge force vector to the element force vector by placing the edge load at correct node index
            for (int n = 0; n < numEdgeNodes; n++) {
                forceVec(3 * edgeNodes[n]) += result(3 * n);
                forceVec(3 * edgeNodes[n] + 1) += result(3 * n + 1);
                forceVec(3 * edgeNodes[n] + 2) += result(3 * n + 2);
            }

        } // After traverse all loaded edges of this element

        // Assemble element force vector to the global force vector
        for (int k = 0; k < elementType; k++) {
            nodalForce(3 * nodeList(k)) += forceVec(3 * k);
            nodalForce(3 * nodeList(k) + 1) += forceVec(3 * k + 1);
            nodalForce(3 * nodeList(k) + 2) += forceVec(3 * k + 2);
        }

    }

    // Apply face load
    for (unsigned i = 0; i < mesh.faceLoadElementList.size(); i++) {
        curr = mesh.getElement(mesh.faceLoadElementList[i]);
        const VectorXi & nodeList = curr->getNodeList();
        int elementType = curr->getSize();
        // The shape function for each face in an isoparametric element is the same
        const std::vector<Vector2d> & gaussianPoint = curr->shape()->faceGaussianPt(); // length 9 vector of Vector2d
        const std::vector<double> & gaussianWeight = curr->shape()->faceGaussianWt(); // length 9 vector
        int numGaussianPts = (int)gaussianPoint.size();

        // The face load information
        const std::vector<int> & faces = mesh.loadFaceList[i];
        const std::vector<double> & loads = mesh.faceLoadValue[i];

        VectorXd forceVec(3 * elementType);
        forceVec.setZero();
        // Traverse all loaded faces in this element
        for (unsigned j = 0; j < faces.size(); j++) {

            Vector3d load(loads[3 * j], loads[3 * j + 1], loads[3 * j + 2]); // 3x1 load vector

            // The global coordinates of the face nodes
            const std::vector<int> & faceNodes = curr->shape()->face(faces[j]);
            int numFaceNodes = (int)faceNodes.size();
            MatrixXd nodeCoord(3, numFaceNodes);
            for (int n = 0; n < numFaceNodes; n++)
                nodeCoord.col(n) = mesh.getNode(nodeList(faceNodes[n]))->getGlobalCoord();

            // Integration at gaussian points
            // sum: N^T * F * |J| * W(i)
            VectorXd result(3 * numFaceNodes); result.setZero(); // 24x1 vector
            for (int g = 0; g < numGaussianPts; g++) {
                MatrixXd N = curr->shape()->faceFunctionMat(g); // 3x24 shape matrix
                MatrixXd globalDeriv = nodeCoord * curr->shape()->faceFunctionDeriv(g); // 3x8 matrix * 8x2 vector = 3x2 matrix [dx/d(xi) dx/d(eta); dy/d(xi) dy/d(eta); dz/d(xi) dz/d(eta)]
                // Here I got into trouble. The jacobian is not a square matrix so there is no determinant for it.
                // How to do surface integral? Check these: https://en.wikipedia.org/wiki/Surface_integral & https://scicomp.stackexchange.com/questions/26886/evaluating-the-surface-integral-in-an-fem-finite-elements-method-procedure & https://cloud.tencent.com/developer/article/1082443
                // |J| = Magnitude (2-norm) of the cross-product vector
                Vector3d v1 = globalDeriv.col(0);
                Vector3d v2 = globalDeriv.col(1);
                Vector3d v3 = v1.cross(v2);
                double jacobianDet = v3.norm();
                result += N.transpose() * load * jacobianDet * gaussianWeight[g]; // 24x1 vector, (Fx,Fy,Fz) at 8 nodes, so 3 * 8 = 24
            }

            // Assemble face force vector to the element force vector by placing the face load at correct node index
            for (int n = 0; n < numFaceNodes; n++) {
                forceVec(3 * faceNodes[n]) += result(3 * n);
                forceVec(3 * faceNodes[n] + 1) += result(3 * n + 1);
                forceVec(3 * faceNodes[n] + 2) += result(3 * n + 2);
            }

        } // After traverse all loaded faces of this element

        // Assemble element force vector to the global force vector
        for (int k = 0; k < elementType; k++) {
            nodalForce(3 * nodeList(k)) += forceVec(3 * k);
            nodalForce(3 * nodeList(k) + 1) += forceVec(3 * k + 1);
            nodalForce(3 * nodeList(k) + 2) += forceVec(3 * k + 2);
        }

    }
}

void Analysis::assembleStiffness()
{
    // Sparse matrix operation notes:
    //     m.setZero() to remove all non-zero coefficients
    //     m.rows() to get number of rows
    //     m.cols() to get number of columns
    //     m.coeffRef(i,j) = k to set value to the element already exists;
    //     m.insert(i,j) = k to set value to the element does not already exist;

    // Reserve the space by estimating the maximum number of non-zero elements
    // in the sparse matrix. This should actually be No. of elements * 3n * 3n.
    // 3n x 3n is the dimension of local stiffness matrix. Note n should be the
    // size of the element type with the MOST nodes (depends on hybrid types). Here n = 20.
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    unsigned long int maxNum = mesh.elementCount() * 60 * 60;
    tripletList.reserve(maxNum);

    const std::vector<int> & DOFList = mesh.boundaryNodeList;
    const std::vector<double> & boundaryValue = mesh.boundaryValue;
    std::unordered_map<unsigned long int, double> boundaryHash;
    boundaryHash.reserve(DOFList.size());
    for (unsigned i = 0; i < DOFList.size(); i++)
        boundaryHash[DOFList[i]] = boundaryValue[i];

    // Assemble global matrix from local matrix of each element, meanwhile modify
    // the stiffness matrix based on boundary condition and adjust the force vector
    // accordingly
    // Old version will first assemble all local stiffness matrix and then deal with
    // boundary conditions, which is slow b/c it needs to redundantly access the
    // sparse global stiffness matrix.
    // New version does things incrementally. When we assemble the local stiffness
    // matrix, we check if it contains a boundary DOF and do subtraction accordingly.
    // Note that Eigen sparse matrix you only need to set the non-zeros. What you
    // don't set will by default be 0.
    Element* curr;
    for (int i = 0; i < mesh.elementCount(); i++) {
        curr = mesh.elementArray()[i];
        int size = curr->getSize();// element type (No. of nodes)
        const VectorXi & nodeList = curr->getNodeList();// the index of nodes belong to this element, e.g., for elementQ4, it will give you a vector contain (10,11,15,14), use this to locate row & column in globalStiffness matrix

        // Bootstrap the computation of local stiffness matrix and force vector. After calling this function, the member variables are all computed
        // @BUG (solved) previous this bootstrap step is in the ctor of derived class ElementQ8, so in the nonlinear analysis, the localStiffness and body & temp force are only computed once at the beginning!
        curr->computeStiffnessAndForce();

        // Traverse each node and assemble the local stiffness values to global stiffness matrix
        const MatrixXd & localStiffness = curr->localStiffness(); // 3n-by-3n matrix
        for (int j = 0; j < size; j++) { // j to access row
            for (int k = 0; k < size; k++) { // k to access col
                // Since the j,k are the node index, actually for every (j,k) we
                // should access a 3x3 block in the local stiffness matrix (3nx3n):
                // (3j,3k), (3j,3k+1), (3j,3k+2) | (3j+1,3k), (3j+1,3k+1), (3j+1,3k+2) | (3j+2,3k), (3j+2,3k+1), (3j+2,3k+2)
                // and place the block into the global stiffness matrix.
                // Note that each block (j,k) in local stiffness matrix will index to
                // (3*nodeList(j)+0/+1/+2, 3*nodeList(k)+0/+1/+2) in global stiffness matrix.

                // For applying the boundary condition, we subtract the column
                // multiplied by the boundary value, and then cross out the column
                // and row at the fixed DOF. So rows are useless and we simply don't
                // assign them into sparse matrix. Columns are subtracted incrementally
                // from the force vector.
                //                           K * U = F
                // for boundary fixed at U(i), we want to zero out ith row & ith column in K & put 1 at K(i,i)
                // 1. to zero out row, just don't assign values in K
                // 2. to zero out column, subtract "col .* boundary value" from F

                // A hash table is used to check the boundary DOF location in constant time.

                // Think from a column-wise perspective
                // First column in local stiffness block (col = 3k, row = 3j & 3j + 1 & 3j + 2)
                if (boundaryHash.find(3 * nodeList(k)) != boundaryHash.end()) { // if the corresponding DOF is on the crossed-out column, subtract it from force vector F
                    // only Local(row = 3j, col = 3k) can possibly be at crossing K(3*nodeList(j), 3*nodeList(k)). And the only way for 3*nodeList(j) == 3*nodeList(k) is by condition j == k
                    // if at crossing, assign 1 at the end, so do nothing here; if not at crossing, subtract
                    if (j != k) // same as if (3 * nodeList(j) != 3 * nodeList(k))
                        nodalForce(3 * nodeList(j)) -= localStiffness(3 * j , 3 * k) * boundaryHash[3 * nodeList(k)];
                    nodalForce(3 * nodeList(j) + 1) -= localStiffness(3 * j + 1, 3 * k) * boundaryHash[3 * nodeList(k)];
                    nodalForce(3 * nodeList(j) + 2) -= localStiffness(3 * j + 2, 3 * k) * boundaryHash[3 * nodeList(k)];
                } else { // not on the crossed-out column, check row index
                    if (boundaryHash.find(3 * nodeList(j)) == boundaryHash.end()) // if not on the crossed-out row,  assemble into K; otherwise do nothing (by default 0)
                        tripletList.push_back(T(3 * nodeList(j), 3 * nodeList(k), localStiffness(3 * j , 3 * k)));
                    if (boundaryHash.find(3 * nodeList(j) + 1) == boundaryHash.end())
                        tripletList.push_back(T(3 * nodeList(j) + 1, 3 * nodeList(k), localStiffness(3 * j + 1 , 3 * k)));
                    if (boundaryHash.find(3 * nodeList(j) + 2) == boundaryHash.end())
                        tripletList.push_back(T(3 * nodeList(j) + 2, 3 * nodeList(k), localStiffness(3 * j + 2 , 3 * k)));
                }

                // Second column in local stiffness block (col = 3k + 1, row = 3j & 3j + 1 & 3j + 2)
                if (boundaryHash.find(3 * nodeList(k) + 1) != boundaryHash.end()) {
                    nodalForce(3 * nodeList(j)) -= localStiffness(3 * j, 3 * k + 1) * boundaryHash[3 * nodeList(k) + 1];
                    // only Local(row = 3j + 1, col = 3k + 1) can possibly be at crossing K(3*nodeList(j)+1, 3*nodeList(k)+1).
                    if (j != k) // same as if (3 * nodeList(j) + 1 != 3 * nodeList(k) + 1)
                        nodalForce(3 * nodeList(j) + 1) -= localStiffness(3 * j + 1, 3 * k + 1) * boundaryHash[3 * nodeList(k) + 1];
                    nodalForce(3 * nodeList(j) + 2) -= localStiffness(3 * j + 2, 3 * k + 1) * boundaryHash[3 * nodeList(k) + 1];
                } else {
                    if (boundaryHash.find(3 * nodeList(j)) == boundaryHash.end())
                        tripletList.push_back(T(3 * nodeList(j), 3 * nodeList(k) + 1, localStiffness(3 * j , 3 * k + 1)));
                    if (boundaryHash.find(3 * nodeList(j) + 1) == boundaryHash.end())
                        tripletList.push_back(T(3 * nodeList(j) + 1, 3 * nodeList(k) + 1, localStiffness(3 * j + 1 , 3 * k + 1)));
                    if (boundaryHash.find(3 * nodeList(j) + 2) == boundaryHash.end())
                        tripletList.push_back(T(3 * nodeList(j) + 2, 3 * nodeList(k) + 1, localStiffness(3 * j + 2 , 3 * k + 1)));
                }

                // Third column in local stiffness block (col = 3k + 2, row = 3j & 3j + 1 & 3j + 2)
                if (boundaryHash.find(3 * nodeList(k) + 2) != boundaryHash.end()) {
                    nodalForce(3 * nodeList(j)) -= localStiffness(3 * j, 3 * k + 2) * boundaryHash[3 * nodeList(k) + 2];
                    nodalForce(3 * nodeList(j) + 1) -= localStiffness(3 * j + 1, 3 * k + 2) * boundaryHash[3 * nodeList(k) + 2];
                    // only Local(row = 3j + 2, col = 3k + 2) can possibly be at crossing K(3*nodeList(j)+2, 3*nodeList(k)+2).
                    if (j != k) // same as if (3 * nodeList(j) + 2 != 3 * nodeList(k) + 2)
                        nodalForce(3 * nodeList(j) + 2) -= localStiffness(3 * j + 2, 3 * k + 2) * boundaryHash[3 * nodeList(k) + 2];
                } else {
                    if (boundaryHash.find(3 * nodeList(j)) == boundaryHash.end())
                        tripletList.push_back(T(3 * nodeList(j), 3 * nodeList(k) + 2, localStiffness(3 * j , 3 * k + 2)));
                    if (boundaryHash.find(3 * nodeList(j) + 1) == boundaryHash.end())
                        tripletList.push_back(T(3 * nodeList(j) + 1, 3 * nodeList(k) + 2, localStiffness(3 * j + 1 , 3 * k + 2)));
                    if (boundaryHash.find(3 * nodeList(j) + 2) == boundaryHash.end())
                        tripletList.push_back(T(3 * nodeList(j) + 2, 3 * nodeList(k) + 2, localStiffness(3 * j + 2 , 3 * k + 2)));
                }

            }
        }

        // Traverse each node and assemble the local nodal force values to global force vector
        const VectorXd & forceVec = curr->nodalForce(); // 3n-by-1 vector
        for (int j = 0; j < size; j++) {
            // Assemble the element-wise body force and temperature load vector
            // at non-boundary DOFs into the global force vector nodalForce;
            // For boundary DOFs, the nodalForce entry will be set to boundary value, i.e. 1 x U = U, at the end along with the K(i,i) = 1 crossing
            // point/edge/face load are already applied in nodalForce by Analysis::applyForce(),
            // so here we apply only the body force and temperature load.
            if (boundaryHash.find(3 * nodeList(j)) == boundaryHash.end()) // check non-boundary DOF
                nodalForce(3 * nodeList(j)) += forceVec(3 * j);
            if (boundaryHash.find(3 * nodeList(j) + 1) == boundaryHash.end())
                nodalForce(3 * nodeList(j) + 1) += forceVec(3 * j + 1);
            if (boundaryHash.find(3 * nodeList(j) + 2) == boundaryHash.end())
                nodalForce(3 * nodeList(j) + 2) += forceVec(3 * j + 2);
        }

    }

    // Edit boundary DOF entries
    // Assign the crossing of boundary DOFs K(i,i) to 1 & F(i) = U(i).
    // s.t. the corresponding location in U and F are both assigned the boundary value, we have 1 x U = U
    // which will always have the exact answer at the restricted boundary DOF
    for (unsigned i = 0; i < DOFList.size(); i++) {
        tripletList.push_back(T(DOFList[i], DOFList[i], 1));
        nodalForce(DOFList[i]) = boundaryValue[i]; // or boundaryHash[DOFList[i]]
    }

    // Convert into sparse matrix
    // Property of setFromTriplets: it will automatically accumulate at repeated locations, so we can feel free to keep push_back above
    // https://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html#acc35051d698e3973f1de5b9b78dbe345
    globalStiffness.setFromTriplets(tripletList.begin(), tripletList.end());
    globalStiffness.makeCompressed();

    // Free the memory
    std::vector<T>().swap(tripletList);

}

void Analysis::computeStrainAndStress()
{
    // Idea:
    // 1. For each element, compute the strain at Gaussian points using nodal
    // displacments by e = Bu. (Note: the strain computed at Gaussian points are accurate!)
    // 2. Let the nodal strain be unknowns, from equation we can solve reversely for
    // e(Gaussian) = N e(nodal) --> e(nodal) = N^-1 e(Gaussian)
    // where N^-1 is the pesudo inverse of the N matrix
    // 3. Compute stress from strain by sigma = E e
    // 4. Cumulate strain and stress value at each node, and average at the end
    // In detail of step 2:
    // 9 Gaussian points & 8 nodes, we have e(Gaussian) as a 9-by-4 matrix.
    // N is the stacked shape function [N1 N2 ... N8] at 9 Gaussian points, i.e.
    // a 9-by-8 matrix. And e(nodal) is the stacked strain vector of 8 nodes as
    // a 8-by-4 matrix. Since N is not square, pesudo inverse is used to solve
    // the least square system.

    // 27 Gaussian points & 20 nodes, we have e(Gaussian) as a 27-by-6 matrix.
    // N is the stacked shape function [N1 N2 ... N20] at 27 Gaussian points, i.e.
    // a 27-by-20 matrix. And e(nodal) is the stacked strain vector of 20 nodes as
    // a 20-by-6 matrix. Since N is not square, pesudo inverse is used to solve
    // the least square system.

    // Traverse each element
    Element* curr;
    int numNodes; // number of nodes belong to the element
    int numGaussianPt; // number of Gaussian points of the element
    for (int i = 0; i < mesh.elementCount(); i++) {
        curr = mesh.elementArray()[i];
        const VectorXi & nodeList = curr->getNodeList();
        numNodes = curr->getSize();
        numGaussianPt = (int)curr->shape()->gaussianPt().size();

        // Assemble the nodal displacement vector for an element (directly from the solved displacement vector)
        VectorXd nodeDisp(3 * numNodes);
        for (int j = 0; j < numNodes; j++) {
            nodeDisp(3 * j) = nodalDisp(3 * nodeList(j));
            nodeDisp(3 * j + 1) = nodalDisp(3 * nodeList(j) + 1);
            nodeDisp(3 * j + 2) = nodalDisp(3 * nodeList(j) + 2);
        }

        MatrixXd strainAtGaussPt(numGaussianPt, 6); // e(Gaussian), 6 for 3D problem
        MatrixXd stressAtGaussPt(numGaussianPt, 6); // sigma(Gaussian)
        MatrixXd shapeAtGaussPt(numGaussianPt, numNodes); // 27x20 matrix N(Gaussian) for element B20, each row is the shape functions [N1...N20] for each Gaussian point

        // Compute strain and stress at gaussian points from e = Bu, sigma = Ee & Assemble into matrix form
        for (int g = 0; g < numGaussianPt; g++) {
            MatrixXd B = curr->BMatrix(curr->shape()->gaussianPt(g));
            VectorXd e = B * nodeDisp; // e(Gaussian) = B * u
            strainAtGaussPt.row(g) = e.transpose();
            VectorXd modulus = (curr->modulusAtGaussPt).row(g); // for nonlinear, this is the stabilized modulus at the Gaussian point; for linear elastic, it's just the constant modulus M
            stressAtGaussPt.row(g) = (curr->EMatrix(modulus) * (e - curr->thermalStrain())).transpose(); // subtract thermal strain, stress = E * (strain - thermal strain)
            shapeAtGaussPt.row(g) = curr->shape()->functionVec(g).transpose();
        }

        // Solve/extrapolate for nodal strain value via a least square linear system using pesudo inverse
        // Current solution: use SVD decomposition
        MatrixXd strainAtNodes = shapeAtGaussPt.bdcSvd(ComputeThinU | ComputeThinV).solve(strainAtGaussPt);
        MatrixXd stressAtNodes = shapeAtGaussPt.bdcSvd(ComputeThinU | ComputeThinV).solve(stressAtGaussPt);
        // Notes on LLS system:
        // Several options for solving a linear least squares system A x = b:
        // 1. SVD decomposition
        // VectorXd x = A.bdcSvd(ComputeThinU | ComputeThinV).solve(b);
        // 2. QR decomposition
        // VectorXd x = A.colPivHouseholderQr().solve(b);
        // from fast to slow, unstable to stable: householderQr()-->colPivHouseholderQr()-->fullPivHouseholderQr()
        // 3. Normal equations (A^T*A)*x = A^T*b
        // VectorXd x = (A.transpose() * A).ldlt().solve(A.transpose() * b)

        // Set calculated strain and stress value to every node (to be accumulated and averaged node-wise later)
        for (int n = 0; n < numNodes; n++)
            mesh.nodeArray()[nodeList(n)]->setStrainAndStress(strainAtNodes.row(n), stressAtNodes.row(n));
    }

}

void Analysis::averageStrainAndStress()
{
    // Average nodal strain and stress, and write the displacement information at the same time!
    for (int i = 0; i < mesh.nodeCount(); i++) {
        mesh.nodeArray()[i]->setDisp(nodalDisp(3 * i), nodalDisp(3 * i + 1), nodalDisp(3 * i + 2));
        nodalStrain.row(i) = mesh.nodeArray()[i]->averageStrain().transpose();
        nodalStress.row(i) = mesh.nodeArray()[i]->averageStress().transpose();
    }
}

void Analysis::printDisp() const
{
    std::cout << "Nodal Displacement: ";
    std::cout << std::endl;
    for (int i = 0; i < mesh.nodeCount(); i++) {
      std::cout << "Node " << i << " : " << mesh.nodeArray()[i]->getDisp().transpose() << std::endl;
    }
    std::cout << std::endl;
}

void Analysis::printStrain() const
{
    std::cout << "Averaged nodal strain: ";
    std::cout << std::endl;
    for (int i = 0; i < mesh.nodeCount(); i++) {
        std::cout << "Node " << i << " : " << nodalStrain.row(i) << std::endl;
    }
    std::cout << std::endl;
}

void Analysis::printStress() const
{
    std::cout << "Averaged nodal stress: ";
    std::cout << std::endl;
    for (int i = 0; i < mesh.nodeCount(); i++) {
        std::cout << "Node " << i << " : " << nodalStress.row(i) << std::endl;
    }
    std::cout << std::endl;
}

void Analysis::writeToVTK(std::string const & fileName) const
{
    // Format: https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
    std::ofstream file(fileName);

    file << "# vtk DataFile Version 2.0\n";
    file << "./vtk\n";
    file << "ASCII\n";
    file << "\n";
    file << "DATASET UNSTRUCTURED_GRID\n";
    file << "POINTS " << mesh.nodeCount() << " double\n";
    for (int i = 0; i < mesh.nodeCount(); i++){
        file << (double)(mesh.nodeArray()[i]->getGlobalCoord())(0) << " " << (double)(mesh.nodeArray()[i]->getGlobalCoord())(1) << " " << (double)(mesh.nodeArray()[i]->getGlobalCoord())(2) << "\n";
    }
    file << "\n";

    // CELLS n size
    // numPoints, i, j, k, ...
    // ...
    // n: number of elements
    // size: number of integers in each row. This includes the starting element type, so "+1"
    file << "CELLS " << mesh.elementCount() << " " << (20 + 1) * mesh.elementCount() << "\n";
    for (int i = 0; i < mesh.elementCount(); i++){
        file << 20 << " " << mesh.elementArray()[i]->getNodeList().transpose() << "\n"; // node index conforms with ParaView
    }
    file << "\n";

    // CELL_TYPES n
    // type0
    // type1
    // ...
    // n: as before
    file << "CELL_TYPES " << mesh.elementCount() << "\n";
    for (int i = 0; i < mesh.elementCount(); i++){
        file << 25 << "\n"; // 25: quadratic hexahedron
    }
    file << "\n";

    file << "POINT_DATA " << mesh.nodeCount() << "\n";

    file << "SCALARS " << "Displacement_X " << "double " << "1" << "\n";
    file << "LOOKUP_TABLE " << "default" << "\n";
    for (int i = 0; i < mesh.nodeCount(); i++){
        file << nodalDisp(3 * i) << "\n";
    }

    file << "SCALARS " << "Displacement_Y " << "double " << "1" << "\n";
    file << "LOOKUP_TABLE " << "default" << "\n";
    for (int i = 0; i < mesh.nodeCount(); i++){
        file << nodalDisp(3 * i + 1) << "\n";
    }

    file << "SCALARS " << "Displacement_Z " << "double " << "1" << "\n";
    file << "LOOKUP_TABLE " << "default" << "\n";
    for (int i = 0; i < mesh.nodeCount(); i++){
        file << nodalDisp(3 * i + 2) << "\n";
    }

    file << "SCALARS " << "Normal_X " << "double " << "1" << "\n";
    file << "LOOKUP_TABLE " << "default" << "\n";
    for (int i = 0; i < mesh.nodeCount(); i++)
        file << nodalStress(i, 0) << "\n";

    file << "SCALARS " << "Normal_Y " << "double " << "1" << "\n";
    file << "LOOKUP_TABLE " << "default" << "\n";
    for (int i = 0; i < mesh.nodeCount(); i++)
        file << nodalStress(i, 1) << "\n";

    file << "SCALARS " << "Normal_Z " << "double " << "1" << "\n";
    file << "LOOKUP_TABLE " << "default" << "\n";
    for (int i = 0; i < mesh.nodeCount(); i++)
        file << nodalStress(i, 2) << "\n";

    file << "SCALARS " << "Shear_XY " << "double " << "1" << "\n";
    file << "LOOKUP_TABLE " << "default" << "\n";
    for (int i = 0; i < mesh.nodeCount(); i++)
        file << nodalStress(i, 3) << "\n";

    file << "SCALARS " << "Shear_YZ " << "double " << "1" << "\n";
    file << "LOOKUP_TABLE " << "default" << "\n";
    for (int i = 0; i < mesh.nodeCount(); i++)
        file << nodalStress(i, 4) << "\n";

    file << "SCALARS " << "Shear_ZX " << "double " << "1" << "\n";
    file << "LOOKUP_TABLE " << "default" << "\n";
    for (int i = 0; i < mesh.nodeCount(); i++)
        file << nodalStress(i, 5) << "\n";

    file << "SCALARS " << "Distance_X " << "double " << "1" << "\n";
    file << "LOOKUP_TABLE " << "default" << "\n";
    for (int i = 0; i < mesh.nodeCount(); i++)
        file << (mesh.nodeArray()[i]->getGlobalCoord())(0)<< "\n";

    file << "SCALARS " << "Distance_Y " << "double " << "1" << "\n";
    file << "LOOKUP_TABLE " << "default" << "\n";
    for (int i = 0; i < mesh.nodeCount(); i++)
        file << (mesh.nodeArray()[i]->getGlobalCoord())(1)<< "\n";

    file << "SCALARS " << "Distance_Z " << "double " << "1" << "\n";
    file << "LOOKUP_TABLE " << "default" << "\n";
    for (int i = 0; i < mesh.nodeCount(); i++)
        file << (mesh.nodeArray()[i]->getGlobalCoord())(2)<< "\n";

    file.close();
}
