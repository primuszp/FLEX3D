/**
 * @file ElementB20.cpp
 * Implementation of ElementB20 class.
 *
 * @author Haohang Huang
 * @date May 19, 2019
 */

#include "ElementB20.h"

// static member should be initialized outside the class body because it doesn't
// depend on any instance of that class and is initialized before any instance is
// created. Initialization here does not require the static member to be public (protected is ok).
// See constructor details in Element.h
ElementB20::staticMembers ElementB20::statics(20,27,12,3,3,6,8,9);

ElementB20::ElementB20()
{
}

ElementB20::ElementB20(const int & index, const std::vector<int> & nodeList, Node** const meshNode, Material* const material)
    : Element(index, nodeList, meshNode, material) // call the constructor of base class in the initializer list!
{
    if (!material->anisotropy) {
        modulusAtGaussPt = VectorXd::Constant(statics.shape->gaussianPt().size(), material_->modulus()); // 27 x 1 vector, the length depends on element type, so can only be initialized in derived class
    } else {
        modulusAtGaussPt = MatrixXd(statics.shape->gaussianPt().size(), 4); // 27 x 4 Matrix. Each column: Mx, My, Mz, shear modulus G
        modulusAtGaussPt.col(0) = VectorXd::Constant(statics.shape->gaussianPt().size(), material_->modulusX());
        modulusAtGaussPt.col(1) = VectorXd::Constant(statics.shape->gaussianPt().size(), material_->modulusY());
        modulusAtGaussPt.col(2) = VectorXd::Constant(statics.shape->gaussianPt().size(), material_->modulusZ());
        modulusAtGaussPt.col(3) = VectorXd::Constant(statics.shape->gaussianPt().size(), material_->modulusG());
    }
}

ElementB20::~ElementB20()
{
}

Shape* ElementB20::shape() const
{
    return statics.shape;
}
