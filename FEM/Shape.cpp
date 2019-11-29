/**
 * @file Shape.cpp
 * Implementation of Shape class.
 *
 * @author Haohang Huang
 * @date May 19, 2019
 */

#include "Shape.h"
#include <cmath>

Shape::Shape(const int & nodes, const int & gaussians, const int & edges, const int & edgeNodes, const int & edgeGaussians, const int & faces, const int & faceNodes, const int & faceGaussians) :
// use initializer list to define the size of vector!
    numNodes_(nodes), nodeCoord_(nodes),
    numGaussianPts_(gaussians), gaussianPt_(gaussians), gaussianWt_(gaussians),
    shapeVec_(gaussians), shapeMat_(gaussians), shapeDeriv_(gaussians),
    numEdges_(edges), numEdgeNodes_(edgeNodes), edgeList_(edges),
    numEdgeGaussianPts_(edgeGaussians), edgeGaussianPt_(edgeGaussians), edgeGaussianWt_(edgeGaussians),
    edgeShapeVec_(edgeGaussians), edgeShapeMat_(edgeGaussians), edgeShapeDeriv_(edgeGaussians),
    numFaces_(faces), numFaceNodes_(faceNodes), faceList_(faces),
    numFaceGaussianPts_(faceGaussians), faceGaussianPt_(faceGaussians), faceGaussianWt_(faceGaussians),
    faceShapeVec_(faceGaussians), faceShapeMat_(faceGaussians), faceShapeDeriv_(faceGaussians)
{
    // Do nothing, only specify the size of each member variable. Computation
    // are done in the constructor of each derived class.
}

const VectorXd & Shape::functionVec(const int & i) const
{
    return shapeVec_[i];
}

const MatrixXd & Shape::functionMat(const int & i) const
{
    return shapeMat_[i];
}

const MatrixXd & Shape::functionDeriv(const int & i) const
{
    return shapeDeriv_[i];
}

const VectorXd & Shape::edgeFunctionVec(const int & i) const
{
    return edgeShapeVec_[i];
}

const MatrixXd & Shape::edgeFunctionMat(const int & i) const
{
    return edgeShapeMat_[i];
}

const VectorXd & Shape::edgeFunctionDeriv(const int & i) const
{
    return edgeShapeDeriv_[i];
}

const VectorXd & Shape::faceFunctionVec(const int & i) const
{
    return faceShapeVec_[i];
}

const MatrixXd & Shape::faceFunctionMat(const int & i) const
{
    return faceShapeMat_[i];
}

const MatrixXd & Shape::faceFunctionDeriv(const int & i) const
{
    return faceShapeDeriv_[i];
}

const std::vector<Vector3d> & Shape::gaussianPt() const
{
    return gaussianPt_;
}

const Vector3d & Shape::gaussianPt(const int & i) const
{
    return gaussianPt_[i];
}

const std::vector<double> & Shape::gaussianWt() const
{
    return gaussianWt_;
}

const double & Shape::gaussianWt(const int & i) const
{
    return gaussianWt_[i];
}

const std::vector<double> & Shape::edgeGaussianPt() const
{
    return edgeGaussianPt_;
}

const std::vector<double> & Shape::edgeGaussianWt() const
{
    return edgeGaussianWt_;
}

const std::vector<int> & Shape::edge(const int & i) const
{
    return edgeList_[i];
}

const std::vector<Vector2d> & Shape::faceGaussianPt() const
{
    return faceGaussianPt_;
}

const std::vector<double> & Shape::faceGaussianWt() const
{
    return faceGaussianWt_;
}

const std::vector<int> & Shape::face(const int & i) const
{
    return faceList_[i];
}

void Shape::_cacheShape()
{
    for (int g = 0; g < numGaussianPts_; g++) {
        shapeVec_[g] = functionVec(gaussianPt_[g]);
        shapeMat_[g] = functionMat(gaussianPt_[g]);
        shapeDeriv_[g] = functionDeriv(gaussianPt_[g]);
    }

    for (int n = 0; n < numEdgeGaussianPts_; n++) {
        edgeShapeVec_[n] = edgeFunctionVec(edgeGaussianPt_[n]);
        edgeShapeMat_[n] = edgeFunctionMat(edgeGaussianPt_[n]);
        edgeShapeDeriv_[n] = edgeFunctionDeriv(edgeGaussianPt_[n]);
    }

    for (int n = 0; n < numFaceGaussianPts_; n++) {
        faceShapeVec_[n] = faceFunctionVec(faceGaussianPt_[n]);
        faceShapeMat_[n] = faceFunctionMat(faceGaussianPt_[n]);
        faceShapeDeriv_[n] = faceFunctionDeriv(faceGaussianPt_[n]);
    }
}
