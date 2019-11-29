/**
 * @file Element.cpp
 * Implementation of Element class.
 *
 * @author Haohang Huang
 * @date May 16, 2019
 */

// #define _USE_MATH_DEFINES // for use of M_PI. 3D doesn't need PI
// #include <cmath>
#include "Element.h"

Element::Element() {}

Element::Element(const int & index, const std::vector<int> & nodeList, Node** const meshNode, Material* const material) // previously here I pass in vector<Node> which is very expensive, now I pass in vector<int> & and a pointer to the pool of nodes create in mesh
  : index_(index), size_(static_cast<int>(nodeList.size())),
    nodeList_(size_), nodeCoord_(size_, 3),
    localStiffness_(MatrixXd::Zero(3 * size_, 3 * size_)),
    nodalForce_(VectorXd::Zero(3 * size_)),
    material_(material) // assign material
{
    for (int i = 0; i < size_; i++) {
      nodeList_(i) = nodeList[i];
      nodeCoord_.row(i) = meshNode[nodeList[i]]->getGlobalCoord();
    }

}

Element::Element(Element const & other)
{
    copy_(other);
}

Element const & Element::operator=(Element const & other)
{
    if (this != &other) {
        clear_();
        copy_(other);
    }
    return *this;
}

Element::~Element()
{
    clear_();
}

Material* Element::material() const
{
    return material_;
}

const MatrixXd & Element::localStiffness() const
{
    return localStiffness_;
}

const VectorXd & Element::nodalForce() const
{
    return nodalForce_;
}

MatrixXd Element::EMatrix(const VectorXd & modulus) const
{
    if (!material_->nonlinearity)
        return material_->EMatrix();
    else
        return material_->EMatrix(modulus);
}

const Vector3d & Element::bodyForce() const
{
    return material_->bodyForce();
}

const VectorXd & Element::thermalStrain() const
{
    return material_->thermalStrain();
}

MatrixXd Element::jacobian(const Vector3d & point) const
{
    return shape()->functionDeriv(point) * nodeCoord_;
}

MatrixXd Element::BMatrix(const Vector3d & point) const
{
    MatrixXd B = MatrixXd::Zero(6, 3 * size_);
    MatrixXd globalDeriv = (shape()->functionDeriv(point) * nodeCoord_).inverse() * shape()->functionDeriv(point);

    // Place into B matrix
    for (int n = 0; n < size_; n++) {
        double Nx = globalDeriv(0, n); // dNi/dx
        double Ny = globalDeriv(1, n); // dNi/dy
        double Nz = globalDeriv(2, n); // dNi/dz
        B(0, 3 * n) =  Nx;
        B(1, 3 * n + 1) = Ny;
        B(2, 3 * n + 2) = Nz;
        B(3, 3 * n) = Ny;
        B(3, 3 * n + 1) = Nx;
        B(4, 3 * n + 1) = Nz;
        B(4, 3 * n + 2) = Ny;
        B(5, 3 * n) = Nz;
        B(5, 3 * n + 2) = Nx;
    }

    return B;
}

void Element::computeStiffnessAndForce()
{
    // @BUG (solved) the same issue as the applyForce() in Analysis class,
    // Initialization!!! Here the localStiffness_ and nodalForce_ are member variables
    // of Element class, so when we iterate on each element, the value will accumulate if
    // we don't propertly initialize it. Initialization in ctor only works for linear elastic
    // but not for nonlinear analysis which requires resetting these variables every time.
    localStiffness_ = MatrixXd::Zero(3 * size_, 3 * size_);
    nodalForce_ = VectorXd::Zero(3 * size_);

    // Loop each Gaussian point and integration by quadrature
    for (int i = 0; i < shape()->gaussianPt().size(); i++) {
        // Local stiffness matrix
        // sum B^T * E * B * |J| * W(i) at all Gaussian points
        localStiffness_ += _BMatrix(i).transpose() * EMatrix(modulusAtGaussPt.row(i)) * _BMatrix(i) * _jacobianDet(i) * shape()->gaussianWt(i);

        // Body force
        // sum N^T * F * |J| * W(i) at all Gaussian points
        nodalForce_ += shape()->functionMat(i).transpose() * bodyForce() * _jacobianDet(i) * shape()->gaussianWt(i);

        // Temperature load
        // sum B^T * E * e0 * |J| * W(i) at all Gaussian points
        nodalForce_ += _BMatrix(i).transpose() * EMatrix(modulusAtGaussPt.row(i)) * thermalStrain() * _jacobianDet(i) * shape()->gaussianWt(i);
    }

}

void Element::computerForce()
{
    // Same as above
    nodalForce_ = VectorXd::Zero(3 * size_);
    for (int i = 0; i < shape()->gaussianPt().size(); i++) {
      // Body force
      // sum N^T * F * |J| * W(i) at all Gaussian points
      nodalForce_ += shape()->functionMat(i).transpose() * bodyForce() * _jacobianDet(i) * shape()->gaussianWt(i);

      // Temperature load
      // sum B^T * E * e0 * |J| * W(i) at all Gaussian points
      nodalForce_ += _BMatrix(i).transpose() * EMatrix(modulusAtGaussPt.row(i)) * thermalStrain() * _jacobianDet(i) * shape()->gaussianWt(i);
    }
}

const int & Element::getIndex() const
{
    return index_;
}

const int & Element::getSize() const
{
    return size_;
}

const VectorXi & Element::getNodeList() const
{
    return nodeList_;
}

const MatrixXd & Element::getNodeCoord() const
{
    return nodeCoord_;
}

MatrixXd Element::_BMatrix(const int & i) const
{
    // Example dimensions are given for ElementQ8 type
    MatrixXd B = MatrixXd::Zero(6, 3 * size_); // 6x24, B matrix
    // [dN1/dx  0   0   | dN2/dx  0   0   | ... ]
    // [0   dN1/dy  0   | 0   dN2/dy  0   | ... ]
    // [0   0   dN1/dz  | 0   0   dN2/dz  | ... ]
    // [dN1/dy dN1/dx 0 | dN2/dy dN2/dx 0 | ... ]
    // [0 dN1/dz dN1/dy | 0 dN2/dz dN2/dy | ... ]
    // [dN1/dz 0 dN1/dx | dN2/dz 0 dN2/dx | ... ]
    // where dN/dx = J^-1 * dN/dxi, dN/dy = J^-1 * dN/deta, dN/dz = J^-1 * dN/dzeta
    // local coordinates (or parent coordinates) are called xi, eta, zeta.

    // Shape function Ni (N1...N8) is actually Ni(xi, eta, zeta). Different Gaussian
    // points have different (xi, eta, zeta) coords, so their Ni's will be all different! As such their B's are different

    // 3x8 global derivatives [dNi/dx; dNi/dy; dNi/dz] = 3x3 inversed Jacobian [J^-1] * 3x8 local derivatives [dNi/dxi; dNi/deta; dNi/dzeta]
    // where 3x3 Jacobian = 3x8 local derivatives [dNi/dxi; dNi/deta; dNi/dzeta] * 8x3 node coordinates [xi yi zi]
    MatrixXd globalDeriv = (shape()->functionDeriv(i) * nodeCoord_).inverse() * shape()->functionDeriv(i);

    // Place into B matrix
    for (int n = 0; n < size_; n++) {
        double Nx = globalDeriv(0, n); // dNi/dx
        double Ny = globalDeriv(1, n); // dNi/dy
        double Nz = globalDeriv(2, n); // dNi/dz
        B(0, 3 * n) =  Nx;
        B(1, 3 * n + 1) = Ny;
        B(2, 3 * n + 2) = Nz;
        B(3, 3 * n) = Ny;
        B(3, 3 * n + 1) = Nx;
        B(4, 3 * n + 1) = Nz;
        B(4, 3 * n + 2) = Ny;
        B(5, 3 * n) = Nz;
        B(5, 3 * n + 2) = Nx;
    }

    return B;
}

double Element::_jacobianDet(const int & i) const
{
    return (shape()->functionDeriv(i) * nodeCoord_).determinant();
}

void Element::clear_()
{ // No dynamically allocated memory within Element class, so do nothing
}

void Element::copy_(Element const & other)
{
    index_ = other.index_;
    size_ = other.size_;
    nodeList_ = other.nodeList_;
    nodeCoord_ = other.nodeCoord_;
    material_ = other.material_;
    localStiffness_ = other.localStiffness_;
    nodalForce_ = other.nodalForce_;

}
