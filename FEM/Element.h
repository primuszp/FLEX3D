/**
 * @file Element.h
 * Abstract base class for various isoparametric elements.
 *
 * @author Haohang Huang
 * @date May 16, 2019
 */

#ifndef Element_h
#define Element_h

#include "Node.h"
#include "Shape.h"
#include "Material.h"
#include "ShapeB20.h"

/* Abstract base Element class with shared public methods and pure virtual methods.
 *
 * For polymorphism of a combination of different element types (e.g., T3, Q6, Q8), we
 * create a base class with virtual methods and enable different implementation
 * using C++ inheritance feature. Methods can be accessed using a generic
 * pointer Element*.
 */
class Element
{
    public:
        /**
         * Default constructor for Element.
         */
        Element();

        /**
         * Custom constructor to create an element with given node indices.
         *
         * @param index The index number of current element.
         * @param nodeList The list of node indices this element is consist of.
         * @param meshNode A const pointer to the node pool of the mesh.
         * @param material A const pointer to the material class.
         *
         * @note Old version of this function pass in all the nodes, which is expensive.
         * This optimized version only pass in a pointer to access the node pool of the mesh.
         */
        Element(const int & index, const std::vector<int> & nodeList, Node** const meshNode, Material* const material);

        /**
         * Copy constructor.
         */
        Element(Element const & other);

        /**
         * Assignment operator.
         */
        Element const & operator=(Element const & other);

        /**
         * Destructor.
         *
         * @note As an abstract class, the destructor must be virtual.
         */
        virtual ~Element();

        /**
         * Get the shape of this element.
         *
         * @return A pointer to the shape object of this element.
         *
         * @note For different element, we would have different shape property,
         * so we return a generic polymorphed Shape class pointer to achieve
         * the standardized behavior of our program.
         * It should also be noted that this is the only pure virtual method
         * in Element class, because the Shape property of each element type is stored
         * as a static member of that derived class and cannot be generalized
         * into a inheritable method.
         * This is really an elegant solution to pave our path towards the hybrid
         * element case. Note that we access using Shape pointer without storing
         * a Shape instance or even a Shape pointer in the Element class, which
         * is also memory efficient for large-scale problems.
         * And returning a const pointer secures the usage of Shape object.
         */
        virtual Shape* shape() const = 0; // return Shape * & is insecure

        /**
         * Get the material information of the element.
         *
         * @return A pointer to the material.
         */
        Material* material() const;

        /**
         * Compute the local element stiffness matrix formulated by
         * differentiation of element strain energy and numerical integration
         * (Gaussian quadrature).
         *
         * @return The 3n-by-3n element stiffness matrix where n is number of
         * nodes belong to this element.
         */
        const MatrixXd & localStiffness() const;

        /**
         * Get the distributed nodal force from body force and temperature loads.
         *
         * @return The nodal force as a 3n-by-1 vector.
         */
        const VectorXd & nodalForce() const;

        /**
         * Get the stress-strain constitutive matrix of the element.
         *
         * @param modulus The modulus value in E matrix. For linear elastic problem,
         * it will just use the constant E; for nonlinear elastic, it will compute
         * E matrix based on specific modulus.
         * @return The 6-by-6 E matrix.
         */
        MatrixXd EMatrix(const VectorXd & modulus) const;

        /**
         * Get the body force (unit weight) of the element.
         *
         * @return The 3-by-1 body force vector.
         */
        const Vector3d & bodyForce() const;

        /**
         * Get the thermal strain to be used in the stress computation.
         *
         * @return The thermal strain as a 6-by-1 vector for 3D problem.
         */
        const VectorXd & thermalStrain() const;

        /**
         * Compute the Jacobian coordinates transformation matrix at an arbitrary
         * point.
         *
         * @param point The point where the Jacobian matrix to be evaluated.
         * @return The 3-by-3 Jacobian matrix.
         *
         * @note There are two sets of jacobian, radius, BMatrix methods in the
         * class. These ones are for calculating at arbitrary point, therefore
         * they are public method; the other ones  calculates at Gaussian points
         * only, which are for the convenience of the formulation of stiffness
         * matrix and force vector by pre-computing them.
         */
        MatrixXd jacobian(const Vector3d & point) const;

        /**
         * Compute the strain-displacement constitutive matrix at an arbitrary point.
         *
         * e = D*u = D*(N*u) = (D*N)u = Bu, so B is the differential matrix
         * multiplying the shape function.
         *
         * @param point The point where the B matrix to be evaluated.
         * @return The 6-by-3n B matrix where n is the number of nodes belong to
         * this element.
         */
        MatrixXd BMatrix(const Vector3d & point) const;

        /** A g-by-1 vector (isotropic) or g-by-4 matrix (anisotropic), where g is the number of Gaussian points of this element */
        MatrixXd modulusAtGaussPt; // for nonlinear analysis, made public for easier access

        /**
         * Helper function for the computation of element stiffness matrix
         * and nodal force vector (body force and temperature load only,
         * since they just depend on the material assigned to the elements).
         */
        void computeStiffnessAndForce();

        /**
         * Helper function for the computation of nodal force vector (body force
         * and temperature load). The assebly of stiffness matrix and force vector
         * are separated due to the possible redundancy introduced in nonlinear
         * scheme.
         */
        void computerForce();

        /**
         * Get the index of this element.
         *
         * @return The index.
         */
        const int & getIndex() const;

        /**
         * Get the size/type of this element.
         *
         * @return The size/type.
         */
        const int & getSize() const;

        /**
         * Get the list of nodes belong to this element.
         *
         * @return The node list as a n-by-1 vector.
         */
        const VectorXi & getNodeList() const;

        /**
         * Get the stacked coordinates matrix of the nodes belong to this element.
         *
         * @return The node coordinates as a n-by-3 matrix.
         */
        const MatrixXd & getNodeCoord() const;

    protected: // make as protected for derived classes to access much easier!

        /* Private helper structure for the bullet-proof memory management of
         * the static members defined in the derived classes.
         *
         * @note This is the most elegant solution to pave our path for a powerful
         * hybrid element system. There are several improvements made throughout
         * the design and development of this module:
         *
         * 1. "Shape" is a general property of a given element type, which does
         * not depend on the instances of that class. In other words, if we
         * add Shape as a class member, we will have 100 Shape objects in memory when
         * we create 100 elements, which is very unwise. Thus, we use a static
         * member that constitues an intrinsic property of the whole class and can be
         * shared by all its instances (with memory just allocated once).
         *
         * 2. At first attempt we defined a static Shape object in the class as:
         *
         * in .h:   static ShapeQ8 shape;
         * in .cpp: ShapeQ8 ElementQ8::shape;
         *
         * However, when we want to use the shape in our hybrid system (for
         * apply load, strain and stress computation, etc), it is inconvenient
         * that we need to figure out which type the element is and then choose
         * the corresponding shape. This requires a if/else check and becomes
         * redundant when element types increase. Thus, a polymorphed Shape
         * pointer is used to provide the flexibility.
         *
         * 3. So far so good. But, in order to use a pointer we need to allocate
         * it on heap by Shape* shape = new ShapeQ8, then where can we delete it
         * to prevent memory leak? One option is to delete this static pointer
         * before we exit the main(), but this means this class is not
         * self-contained and requires the client to be careful enough to delete
         * it at the end, which no one can promise. Thus, we use an internal
         * structure to manage all the static members we have in the class, and
         * the constructor & destructor of this structure takes over the
         * create & destroy task for us. Before the program starts, the static
         * structure is initialized by its constructor; when the program exists,
         * the structure is destroyed and its destructor is called so that the
         * shape object on heap is deallocated. This is a perfect solution for
         * memory management.
         *
         * It is also good to note that a nested class (class/struct defined inside a
         * class) can be inherited as well.
         */
        struct staticMembers {
              /** A pointer to the shape instance of this element type. */
              Shape* shape; // we can have more static members, if needed

              /**
               * Constructor for the strucuture.
               *
               * @param nodes Number of nodes of this shape of element.
               * @param gaussians Number of Gaussian integration points of this
               * shape of element.
               * @param edges Number of edges of this shape of element.
               * @param edgeNodes Number of nodes each edge consists of.
               * @param edgeGaussians Number of Gaussian intergration points at each
               * edge.
               */
              staticMembers(int nodes, int gaussians, int edges, int edgeNodes, int edgeGaussians, int faces, int faceNodes, int faceGaussians) {
                  // Create instances of different types of shape
                  switch (nodes) {
                      case 20 :
                          shape = new ShapeB20(nodes, gaussians, edges, edgeNodes, edgeGaussians, faces, faceNodes, faceGaussians);
                          break;
                  }
              }

              /**
               * Destructor for the strucuture. Self-contained class/struct.
               */
              ~staticMembers() {
                  delete shape; shape = NULL;
              }
          };

        /** Zero-based index of the element */
        int index_;

        /** Size/Type of the element, 8 for Q8, 6 for Q6, etc */
        int size_;

        /** The list of the nodes belong to this element */
        VectorXi nodeList_;

        /** The global coordinates of the nodes belong to this element as a
         * n-by-3 matrix.
         */
        MatrixXd nodeCoord_;

        /** The 3n-by-3n local stiffness matrix where n is the number of nodes */
        MatrixXd localStiffness_; // member variable cannot have the same name as member function

        /** The 3n-by-1 nodal force vector where n is the number of nodes */
        VectorXd nodalForce_;

        /**
         * A pointer to the material class. Material should not be stored as static
         * member as Shape does, because Shape is a generic property of a Element
         * while the Material can vary among the same element type. In addition,
         * the list of all material is maintained at the Mesh level, so the
         * memory allocation and deallocation is not handled inside Element. */
        Material* material_;

        /**
         * Private helper function for computing the B matrix (the strain-displacement
         * transformation matrix, e = D * u = D * N * u = B * u) at ith Gaussian point.
         *
         * @param i The ith Gaussian point where the B matrix to be evaluated.
         * @return The B matrix.
         * @note For the BMatrix at arbitrary point, use the public method BMatrix().
         */
        MatrixXd _BMatrix(const int & i) const;

        /**
         * Private helper function for computing the determinant of Jacobian
         * matrix at ith Gaussian point.
         *
         * @param i The Gaussian point where the Jacobian matrix to be evaluated.
         * @return The determinant of Jacobian matrix.
         */
        double _jacobianDet(const int & i) const;

        /**
         * Private helper function for deleting the current element, used in destructor
         * and assignment operator.
         */
        void clear_();

        /**
         * Private helper function for deep-copying the element, used in copy
         * constructor and assignment operator.
         */
        void copy_(Element const & other);

};

#endif /* Element_h */
