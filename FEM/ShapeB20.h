/**
 * @file ShapeB20.h
 * Derived class from Shape for the shape information of isoparametric B20 element.
 *
 * @author Haohang Huang
 * @date May 19, 2019
 */

#ifndef ShapeB20_h
#define ShapeB20_h

#include "Shape.h"

/* Derived class for the shape of isoparametric B20 element. See shape declaration in ElementB20.cpp.
 */
class ShapeB20 : public Shape
{
    public:
        /* See the documentation of base class Shape. */
        ShapeB20(const int & nodes, const int & gaussians, const int & edges, const int & edgeNodes, const int & edgeGaussians, const int & faces, const int & faceNodes, const int & faceGaussians);
        ~ShapeB20();

        VectorXd functionVec(const Vector3d & point) const;
        MatrixXd functionMat(const Vector3d & point) const;
        MatrixXd functionDeriv(const Vector3d & point) const;

        VectorXd edgeFunctionVec(const double & point) const;
        MatrixXd edgeFunctionMat(const double & point) const;
        VectorXd edgeFunctionDeriv(const double & point) const;

        VectorXd faceFunctionVec(const Vector2d & point) const;
        MatrixXd faceFunctionMat(const Vector2d & point) const;
        MatrixXd faceFunctionDeriv(const Vector2d & point) const;

};

#endif /* ShapeB20_h */
