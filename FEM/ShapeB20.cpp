/**
 * @file ShapeB20.cpp
 * Implementation of ShapeB20 class.
 *
 * @author Haohang Huang
 * @date May 19, 2019
 */

#include "ShapeB20.h"
#include <cmath>

ShapeB20::ShapeB20(const int & nodes, const int & gaussians, const int & edges, const int & edgeNodes, const int & edgeGaussians, const int & faces, const int & faceNodes, const int & faceGaussians) :
    Shape(nodes, gaussians, edges, edgeNodes, edgeGaussians, faces, faceNodes, faceGaussians) // Call the constructor of base class
{
    /*
    * The sketch and index of the B20 element conforms with ParaView format.
    *
    * z = +1 level:
    * 4 -- 15 -- 7
    * |          |
    * 12         14
    * |          |
    * 5 -- 13 -- 6
    *
    * z = 0 level:
    * 16 ------- 19
    * |          |
    * |          |
    * |          |
    * 17 ------- 18
    *
    * z = -1 level:
    * 0 -- 11 -- 3
    * |          |
    * 8         10
    * |          |
    * 1 -- 9 --  2
    */
    // Local xi-eta-zeta coordinates of nodes
    // Corner nodes
    nodeCoord_[0] << -1, -1, -1;
    nodeCoord_[1] << 1, -1, -1;
    nodeCoord_[2] << 1, 1, -1;
    nodeCoord_[3] << -1, 1, -1;
    nodeCoord_[4] << -1, -1, 1;
    nodeCoord_[5] << 1, -1, 1;
    nodeCoord_[6] << 1, 1, 1;
    nodeCoord_[7] << -1, 1, 1;
    // Mid-side nodes
    nodeCoord_[8] << 0, -1, -1;
    nodeCoord_[9] << 1, 0, -1;
    nodeCoord_[10] << 0, 1, -1;
    nodeCoord_[11] << -1, 0, -1;
    nodeCoord_[12] << 0, -1, 1;
    nodeCoord_[13] << 1, 0, 1;
    nodeCoord_[14] << 0, 1, 1;
    nodeCoord_[15] << -1, 0, 1;
    nodeCoord_[16] << -1, -1, 0;
    nodeCoord_[17] << 1, -1, 0;
    nodeCoord_[18] << 1, 1, 0;
    nodeCoord_[19] << -1, 1, 0;

    // Local xi-eta-zeta coordinates of gaussian points
    // Ordering of Gaussian points doesn't matter that much, just ensure point & weight match each other
    double temp = std::sqrt(0.6);
    gaussianPt_[0] << -temp, -temp, -temp;
    gaussianPt_[1] << -temp, -temp, 0;
    gaussianPt_[2] << -temp, -temp, temp;
    gaussianPt_[3] << -temp, 0, -temp;
    gaussianPt_[4] << -temp, 0, 0;
    gaussianPt_[5] << -temp, 0, temp;
    gaussianPt_[6] << -temp, temp, -temp;
    gaussianPt_[7] << -temp, temp, 0;
    gaussianPt_[8] << -temp, temp, temp;
    gaussianPt_[9] << 0, -temp, -temp;
    gaussianPt_[10] << 0, -temp, 0;
    gaussianPt_[11] << 0, -temp, temp;
    gaussianPt_[12] << 0, 0, -temp;
    gaussianPt_[13] << 0, 0, 0;
    gaussianPt_[14] << 0, 0, temp;
    gaussianPt_[15] << 0, temp, -temp;
    gaussianPt_[16] << 0, temp, 0;
    gaussianPt_[17] << 0, temp, temp;
    gaussianPt_[18] << temp, -temp, -temp;
    gaussianPt_[19] << temp, -temp, 0;
    gaussianPt_[20] << temp, -temp, temp;
    gaussianPt_[21] << temp, 0, -temp;
    gaussianPt_[22] << temp, 0, 0;
    gaussianPt_[23] << temp, 0, temp;
    gaussianPt_[24] << temp, temp, -temp;
    gaussianPt_[25] << temp, temp, 0;
    gaussianPt_[26] << temp, temp, temp;

    // Gaussian weights
    // Note: corner/side here means there are 3 Gaussian points along each x/y/z (or xi/eta/zeta) direction, like this
    // x --------- x --------- x
    // corner --- side --- corner
    // "side" actually means "middle", so middel Gaussian points has larger weight (8/9)
    // "corner" means "end", so end Gaussian points has smaller weight (5/9)
    // The weight at a 3D Gaussian point is just multiply its relative side/corner location in EACH direction
    // e.g, the Gaussian on an edge will be "side" for one direction and "corner" for the other two directions, therefore side * corner * corner
    double corner = 5.0 / 9.0; // do not use 5/9! that is 0 due to rounding down!
    double side = 8.0 / 9.0;
    gaussianWt_[0] = corner * corner * corner;
    gaussianWt_[1] = corner * corner * side;
    gaussianWt_[2] = corner * corner * corner;
    gaussianWt_[3] = corner * corner * side;
    gaussianWt_[4] = corner * side * side;
    gaussianWt_[5] = corner * corner * side;
    gaussianWt_[6] = corner * corner * corner;
    gaussianWt_[7] = corner * corner * side;
    gaussianWt_[8] = corner * corner * corner;
    gaussianWt_[9] = corner * corner * side;
    gaussianWt_[10] = corner * side * side;
    gaussianWt_[11] = corner * corner * side;
    gaussianWt_[12] = corner * side * side;
    gaussianWt_[13] = side * side * side; // center
    gaussianWt_[14] = corner * side * side;
    gaussianWt_[15] = corner * corner * side;
    gaussianWt_[16] = corner * side * side;
    gaussianWt_[17] = corner * corner * side;
    gaussianWt_[18] = corner * corner * corner;
    gaussianWt_[19] = corner * corner * side;
    gaussianWt_[20] = corner * corner * corner;
    gaussianWt_[21] = corner * corner * side;
    gaussianWt_[22] = corner * side * side;
    gaussianWt_[23] = corner * corner * side;
    gaussianWt_[24] = corner * corner * corner;
    gaussianWt_[25] = corner * corner * side;
    gaussianWt_[26] = corner * corner * corner;

    // Edge Gaussian points (for applying edge load)
    // 0 ---- 1 ----- 2
    edgeGaussianPt_[0] = -temp;
    edgeGaussianPt_[1] = 0;
    edgeGaussianPt_[2] = temp;

    // Edge Gaussian weights
    edgeGaussianWt_[0] = corner;
    edgeGaussianWt_[1] = side; // center
    edgeGaussianWt_[2] = corner;

    // Edge list (see ShapeB20.h, Numbering from bottom up, counter-clockwise, node index should be consistent with the Gaussian weights)
    // Horizontal edges
    std::vector<int> edge0{0,8,1};
    std::vector<int> edge1{1,9,2};
    std::vector<int> edge2{3,10,2};
    std::vector<int> edge3{0,11,3};
    std::vector<int> edge4{4,12,5};
    std::vector<int> edge5{5,13,6};
    std::vector<int> edge6{7,14,6};
    std::vector<int> edge7{4,15,7};
    // Vertical edges
    std::vector<int> edge8{4,16,0};
    std::vector<int> edge9{5,17,1};
    std::vector<int> edge10{6,18,2};
    std::vector<int> edge11{7,19,3};
    edgeList_[0] = edge0;
    edgeList_[1] = edge1;
    edgeList_[2] = edge2;
    edgeList_[3] = edge3;
    edgeList_[4] = edge4;
    edgeList_[5] = edge5;
    edgeList_[6] = edge6;
    edgeList_[7] = edge7;
    edgeList_[8] = edge8;
    edgeList_[9] = edge9;
    edgeList_[10] = edge10;
    edgeList_[11] = edge11;

    // Face Gaussian points (for applying face load)
    // 6 -- 7 -- 8
    // |    |    |
    // 3 -- 4 -- 5
    // |    |    |
    // 0 -- 1 -- 2
    faceGaussianPt_[0] << -temp, -temp;
    faceGaussianPt_[1] << 0, -temp;
    faceGaussianPt_[2] << temp, -temp;
    faceGaussianPt_[3] << -temp, 0;
    faceGaussianPt_[4] << 0, 0;
    faceGaussianPt_[5] << temp, 0;
    faceGaussianPt_[6] << -temp, temp;
    faceGaussianPt_[7] << 0, temp;
    faceGaussianPt_[8] << temp, temp;

    // Face Gaussian weights
    faceGaussianWt_[0] = corner * corner;
    faceGaussianWt_[2] = corner * corner;
    faceGaussianWt_[6] = corner * corner;
    faceGaussianWt_[8] = corner * corner;
    faceGaussianWt_[1] = side * corner;
    faceGaussianWt_[3] = side * corner;
    faceGaussianWt_[5] = side * corner;
    faceGaussianWt_[7] = side * corner;
    faceGaussianWt_[4] = side * side; // center

    // Face list
    // 3 -- 6 -- 2
    // |         |
    // 7         5
    // |         |
    // 0 -- 4 -- 1
    std::vector<int> face0{1,2,3,0,9,10,11,8};
    std::vector<int> face1{5,6,7,4,13,14,15,12};
    std::vector<int> face2{0,1,5,4,8,17,12,16};
    std::vector<int> face3{1,2,6,5,9,18,13,17};
    std::vector<int> face4{2,3,7,6,10,19,14,18};
    std::vector<int> face5{0,3,7,4,11,19,15,16};
    faceList_[0] = face0; // bottom
    faceList_[1] = face1; // top
    faceList_[2] = face2; // left
    faceList_[3] = face3; // front
    faceList_[4] = face4; // right
    faceList_[5] = face5; // back

    // After setting up the above parameters, cache the shape function by pre-computing them
    _cacheShape();

}

ShapeB20::~ShapeB20()
{
}

VectorXd ShapeB20::functionVec(const Vector3d & point) const
{ // 20x1 Vector
    VectorXd result(numNodes_);
    double Ni = 0;
    for (int i = 0; i < numNodes_; i++) {
        if (i < 8) { // 8 corner nodes
          // Ni = (1+xi_i*xi)(1+eta_i*eta)(1+zeta_i*zeta)(xi_i*xi+eta_i*eta+zeta_i*zeta-2)/8
          Ni = (1 + nodeCoord_[i](0) * point(0)) * (1 + nodeCoord_[i](1) * point(1)) * (1 + nodeCoord_[i](2) * point(2)) * (nodeCoord_[i](0) * point(0) + nodeCoord_[i](1) * point(1) + nodeCoord_[i](2) * point(2) - 2) / 8;
        }
        if (i == 8 || i == 10 || i == 12 || i == 14) { // xi = 0 mid-side nodes
          // Ni = (1-xi^2)(1+eta_i*eta)(1+zeta_i*zeta)/4
          Ni = (1 - point(0) * point(0)) * (1 + nodeCoord_[i](1) * point(1)) * (1 + nodeCoord_[i](2) * point(2)) / 4;
        }
        if (i == 9 || i == 11 || i == 13 || i == 15) { // eta = 0 mid-side nodes
          // Ni = (1+xi_i*xi)(1-eta^2)(1+zeta_i*zeta)/4
          Ni = (1 + nodeCoord_[i](0) * point(0)) * (1 - point(1) * point(1)) * (1 + nodeCoord_[i](2) * point(2)) / 4;
        }
        if (i == 16 || i == 17 || i == 18 || i == 19) { // zeta = 0 mid-side nodes
          // Ni = (1+xi_i*xi)(1+eta_i*eta)(1-zeta^2)/4
          Ni = (1 + nodeCoord_[i](0) * point(0)) * (1 + nodeCoord_[i](1) * point(1)) * (1 - point(2) * point(2)) / 4;
        }
        result(i) = Ni;
    }
    return result;
}

MatrixXd ShapeB20::functionMat(const Vector3d & point) const
{ // 3x60 matrix
    MatrixXd result = MatrixXd::Zero(3, 3 * numNodes_);
    double Ni = 0;
    for (int i = 0; i < numNodes_; i++) {
        if (i < 8) { // 8 corner nodes
          // Ni = (1+xi_i*xi)(1+eta_i*eta)(1+zeta_i*zeta)(xi_i*xi+eta_i*eta+zeta_i*zeta-2)/8
          Ni = (1 + nodeCoord_[i](0) * point(0)) * (1 + nodeCoord_[i](1) * point(1)) * (1 + nodeCoord_[i](2) * point(2)) * (nodeCoord_[i](0) * point(0) + nodeCoord_[i](1) * point(1) + nodeCoord_[i](2) * point(2) - 2) / 8;
        }
        if (i == 8 || i == 10 || i == 12 || i == 14) { // xi = 0 mid-side nodes
          // Ni = (1-xi^2)(1+eta_i*eta)(1+zeta_i*zeta)/4
          Ni = (1 - point(0) * point(0)) * (1 + nodeCoord_[i](1) * point(1)) * (1 + nodeCoord_[i](2) * point(2)) / 4;
        }
        if (i == 9 || i == 11 || i == 13 || i == 15) { // eta = 0 mid-side nodes
          // Ni = (1+xi_i*xi)(1-eta^2)(1+zeta_i*zeta)/4
          Ni = (1 + nodeCoord_[i](0) * point(0)) * (1 - point(1) * point(1)) * (1 + nodeCoord_[i](2) * point(2)) / 4;
        }
        if (i == 16 || i == 17 || i == 18 || i == 19) { // zeta = 0 mid-side nodes
          // Ni = (1+xi_i*xi)(1+eta_i*eta)(1-zeta^2)/4
          Ni = (1 + nodeCoord_[i](0) * point(0)) * (1 + nodeCoord_[i](1) * point(1)) * (1 - point(2) * point(2)) / 4;
        }
        result(0, 3 * i) = Ni;
        result(1, 3 * i + 1) = Ni;
        result(2, 3 * i + 2) = Ni;
    }
    return result;
}

MatrixXd ShapeB20::functionDeriv(const Vector3d & point) const
{ // 3x20 matrix
    MatrixXd result(3, numNodes_);
    for (int i = 0; i < numNodes_; i++) {
      if (i < 8) { // 8 corner nodes
        // d(Ni)/d(xi)=xi_i*(1+eta_i*eta)*(1+zeta_i*zeta)*(2*xi_i*xi+eta_i*eta+zeta_i*zeta-1)/8
        // d(Ni)/d(eta), d(Ni)/d(zeta) symmetric
        result(0, i) = nodeCoord_[i](0) * (1 + nodeCoord_[i](1) * point(1)) * (1 + nodeCoord_[i](2) * point(2)) * (2 * nodeCoord_[i](0) * point(0) + nodeCoord_[i](1) * point(1) + nodeCoord_[i](2) * point(2) - 1)/ 8;
        result(1, i) = nodeCoord_[i](1) * (1 + nodeCoord_[i](0) * point(0)) * (1 + nodeCoord_[i](2) * point(2)) * (2 * nodeCoord_[i](1) * point(1) + nodeCoord_[i](0) * point(0) + nodeCoord_[i](2) * point(2) - 1)/ 8;
        result(2, i) = nodeCoord_[i](2) * (1 + nodeCoord_[i](0) * point(0)) * (1 + nodeCoord_[i](1) * point(1)) * (2 * nodeCoord_[i](2) * point(2) + nodeCoord_[i](0) * point(0) + nodeCoord_[i](1) * point(1) - 1)/ 8;
      }
      if (i == 8 || i == 10 || i == 12 || i == 14) { // xi = 0 mid-side nodes
        // d(Ni)/d(xi)=-xi*(1+eta_i*eta)*(1+zeta_i*zeta)/2
        result(0, i) = - point(0) * (1 + nodeCoord_[i](1) * point(1)) * (1 + nodeCoord_[i](2) * point(2)) / 2;
        // d(Ni)/d(eta)=(1-xi^2)*eta_i*(1+zeta_i*zeta)/4
        result(1, i) = (1 - point(0) * point(0)) * nodeCoord_[i](1) * (1 + nodeCoord_[i](2) * point(2)) / 4;
        // d(Ni)/d(zeta)=(1-xi^2)*(1+eta_i*eta)*zeta_i/4
        result(2, i) = (1 - point(0) * point(0)) * (1 + nodeCoord_[i](1) * point(1)) * nodeCoord_[i](2) / 4;
      }
      if (i == 9 || i == 11 || i == 13 || i == 15) { // eta = 0 mid-side nodes
        // d(Ni)/d(xi)=xi_i*(1-eta^2)*(1+zeta_i*zeta)/4
        result(0, i) = nodeCoord_[i](0) * (1 - point(1) * point(1)) * (1 + nodeCoord_[i](2) * point(2)) / 4;
        // d(Ni)/d(eta)=-(1+xi_i*xi)*eta*(1+zeta_i*zeta)/2
        result(1, i) = - (1 + nodeCoord_[i](0) * point(0)) * point(1) * (1 + nodeCoord_[i](2) * point(2)) / 2;
        // d(Ni)/d(zeta)=(1+xi_i*xi)*(1-eta^2)*zeta_i/4
        result(2, i) = (1 + nodeCoord_[i](0) * point(0)) * (1 - point(1) * point(1)) * nodeCoord_[i](2) / 4;
      }
      if (i == 16 || i == 17 || i == 18 || i == 19) { // zeta = 0 mid-side nodes
        // d(Ni)/d(xi)=xi_i*(1+eta_i*eta)*(1-zeta^2)/4
        result(0, i) = nodeCoord_[i](0) * (1 + nodeCoord_[i](1) * point(1)) * (1 - point(2) * point(2)) / 4;
        // d(Ni)/d(eta)=(1+xi_i*xi)*eta_i*(1-zeta^2)/4
        result(1, i) = (1 + nodeCoord_[i](0) * point(0)) * nodeCoord_[i](1) * (1 - point(2) * point(2)) / 4;
        // d(Ni)/d(zeta)=-(1+xi_i*xi)*(1+eta_i*eta)*zeta/2
        result(2, i) = - (1 + nodeCoord_[i](0) * point(0)) * (1 + nodeCoord_[i](1) * point(1)) * point(2) / 2;
      }
    }
    return result;
}

VectorXd ShapeB20::edgeFunctionVec(const double & point) const
{ // 3x1 vector
    VectorXd result(numEdgeNodes_);
    double Ni = 0;
    for (int i = 0; i < numEdgeNodes_; i++) {
        if (i == 0) // left/bottom point
          Ni = - point * (1 - point) / 2; // Ni = -x(1-x)/2
        if (i == 2) // right/top point
          Ni = point * (1 + point) / 2; // Ni = x(1+x)/2
        if (i == 1) // mid-side point
          Ni = 1 - point * point; // Ni = (1-x^2)
        result(i) = Ni;
    }
    return result;
}

MatrixXd ShapeB20::edgeFunctionMat(const double & point) const
{ // 3x9 matrix, 3 gaussian point and 3 node at each edge
    MatrixXd result = MatrixXd::Zero(3, 3 * numEdgeNodes_);
    double Ni = 0;
    for (int i = 0; i < numEdgeNodes_; i++) {
        if (i == 0) // left/bottom point
          Ni = - point * (1 - point) / 2; // Ni = -x(1-x)/2
        if (i == 2) // right/top point
          Ni = point * (1 + point) / 2; // Ni = x(1+x)/2
        if (i == 1) // mid-side point
          Ni = 1 - point * point; // Ni = (1-x^2)
        result(0, 3 * i) = Ni;
        result(1, 3 * i + 1) = Ni;
        result(2, 3 * i + 2) = Ni;
    }
    return result;
}

VectorXd ShapeB20::edgeFunctionDeriv(const double & point) const
{ // 3x1 vector
    VectorXd result(numEdgeNodes_);
    double Ni = 0;
    for (int i = 0; i < numEdgeNodes_; i++) {
        if (i == 0) // left/bottom point
          Ni = (2 * point - 1) / 2; // dNi/dx = (2x-1)/2
        if (i == 2) // right/top point
          Ni = (2 * point + 1) / 2; // dNi/dx = (2x+1)/2
        if (i == 1) // mid-side point
          Ni = - 2 * point; // dNi/dx = -2x
        result(i) = Ni;
    }
    return result;
}

VectorXd ShapeB20::faceFunctionVec(const Vector2d & point) const
{ // 8x1 vector. This is exactly the same in ShapeQ8 b/c a face is 2D

    // Local xi-eta coordinates of nodes
    std::vector<Vector2d> nodeCoord(numFaceNodes_);
    // Corner nodes
    nodeCoord[0] << -1, -1;
    nodeCoord[1] << 1, -1;
    nodeCoord[2] << 1, 1;
    nodeCoord[3] << -1, 1;
    // Mid-side nodes
    nodeCoord[4] << 0, -1;
    nodeCoord[5] << 1, 0;
    nodeCoord[6] << 0, 1;
    nodeCoord[7] << -1, 0;

    VectorXd result(numFaceNodes_);
    double Ni = 0;
    for (int i = 0; i < numFaceNodes_; i++) {
        if (i < 4) { // 4 corner nodes
          // Ni = (1+xi_i*xi)(1+eta_i*eta)(xi_i*xi+eta_i*eta-1)/4
          Ni = (1 + nodeCoord[i](0) * point(0)) * (1 + nodeCoord[i](1) * point(1)) * (nodeCoord[i](0) * point(0) + nodeCoord[i](1) * point(1) - 1) / 4;
        }
        if (i == 4 || i == 6) { // xi = 0 mid-side nodes
          // Ni = (1-xi^2)(1+eta_i*eta)/2
          Ni = (1 - point(0) * point(0)) * (1 + nodeCoord[i](1) * point(1)) / 2;
        }
        if (i== 5 || i == 7) { // eta = 0 mid-side nodes
          // Ni = (1+xi_i*xi)(1-eta^2)/2
          Ni = (1 + nodeCoord[i](0) * point(0)) * (1 - point(1) * point(1)) / 2;
        }
        result(i) = Ni;
    }
    return result;
}

MatrixXd ShapeB20::faceFunctionMat(const Vector2d & point) const
{ // 3x24 matrix

    // Local xi-eta coordinates of nodes
    std::vector<Vector2d> nodeCoord(numFaceNodes_);
    // Corner nodes
    nodeCoord[0] << -1, -1;
    nodeCoord[1] << 1, -1;
    nodeCoord[2] << 1, 1;
    nodeCoord[3] << -1, 1;
    // Mid-side nodes
    nodeCoord[4] << 0, -1;
    nodeCoord[5] << 1, 0;
    nodeCoord[6] << 0, 1;
    nodeCoord[7] << -1, 0;

    MatrixXd result = MatrixXd::Zero(3, 3 * numFaceNodes_);
    double Ni = 0;
    for (int i = 0; i < numFaceNodes_; i++) {
        if (i < 4) { // 4 corner nodes
          // Ni = (1+xi_i*xi)(1+eta_i*eta)(xi_i*xi+eta_i*eta-1)/4
          Ni = (1 + nodeCoord[i](0) * point(0)) * (1 + nodeCoord[i](1) * point(1)) * (nodeCoord[i](0) * point(0) + nodeCoord[i](1) * point(1) - 1) / 4;
        }
        if (i == 4 || i == 6) { // xi = 0 mid-side nodes
          // Ni = (1-xi^2)(1+eta_i*eta)/2
          Ni = (1 - point(0) * point(0)) * (1 + nodeCoord[i](1) * point(1)) / 2;
        }
        if (i == 5 || i == 7) { // eta = 0 mid-side nodes
          // Ni = (1+xi_i*xi)(1-eta^2)/2
          Ni = (1 + nodeCoord[i](0) * point(0)) * (1 - point(1) * point(1)) / 2;
        }
        result(0, 3 * i) = Ni;
        result(1, 3 * i + 1) = Ni;
        result(2, 3 * i + 2) = Ni;
    }
    return result;
}

MatrixXd ShapeB20::faceFunctionDeriv(const Vector2d & point) const
{ // 8x2 matrix

    // Local xi-eta coordinates of nodes
    std::vector<Vector2d> nodeCoord(numFaceNodes_);
    // Corner nodes
    nodeCoord[0] << -1, -1;
    nodeCoord[1] << 1, -1;
    nodeCoord[2] << 1, 1;
    nodeCoord[3] << -1, 1;
    // Mid-side nodes
    nodeCoord[4] << 0, -1;
    nodeCoord[5] << 1, 0;
    nodeCoord[6] << 0, 1;
    nodeCoord[7] << -1, 0;

    MatrixXd result(numFaceNodes_, 2);
    for (int i = 0; i < numFaceNodes_; i++) {
      if (i < 4) { // 4 corner nodes
        // d(Ni)/d(xi)=xi_i*(1+eta_i*eta)*(2*xi_i*xi+eta_i*eta)/4
        result(i, 0) = nodeCoord[i](0) * (1 + nodeCoord[i](1) * point(1)) * (2 * nodeCoord[i](0) * point(0) + nodeCoord[i](1) * point(1))/ 4;
        // d(Ni)/d(eta)=eta_i*(1+xi_i*xi)*(2*eta_i*eta+xi_i*xi)/4
        result(i, 1) = nodeCoord[i](1) * (1 + nodeCoord[i](0) * point(0)) * (2 * nodeCoord[i](1) * point(1) + nodeCoord[i](0) * point(0))/ 4;
      }
      if (i == 4 || i == 6) { // xi = 0 mid-side nodes
        // d(Ni)/d(xi)= -xi*(1+eta_i*eta)
        result(i, 0) = - point(0) * (1 + nodeCoord[i](1) * point(1));
        // d(Ni)/d(eta)= (1-xi^2)*eta_i/2
        result(i, 1) = (1 - point(0) * point(0)) * nodeCoord[i](1) / 2;
      }
      if (i== 5 || i == 7) { // eta = 0 mid-side nodes
        // d(Ni)/d(xi)= (1-eta^2)*xi_i/2
        result(i, 0) = (1 - point(1) * point(1)) * nodeCoord[i](0) / 2;
        // d(Ni)/d(eta)= -eta*(1+xi_i*xi)
        result(i, 1) = - point(1) * (1 + nodeCoord[i](0) * point(0));
      }
    }
    return result;
}
