/**
 * @file Node.h
 * Node class for the basic coordinates and mechanical properties.
 *
 * @author Haohang Huang
 * @date May 16, 2019
 */

#ifndef Node_h
#define Node_h

#include "Eigen/Eigen"

using namespace Eigen;

/* Node class for the basic coordinates and mechanical properties information.
 */
class Node
{
    public:
        /**
         * Default constructor for Node.
         */
        Node();

        /**
         * Custom constructor to create a node with given node coordinates.
         *
         * @param index The index number of current node.
         * @param x The global x-coordinate of the node.
         * @param y The global y-coordinate of the node.
         * @param z The global z-coordinate of the node.
         * @note All pass-by-ref.
         */
        Node(const int & index, const double & x, const double & y, const double & z);

        /**
         * Copy constructor.
         */
        Node(Node const & other);

        /**
         * Assignment operator.
         */
        Node const & operator=(Node const & other);

        /**
         * Destructor.
         */
        ~Node();

        /**
         * Assign the initial global coordinates (x,y,z) of the node.
         *
         * @param x The global x-coordinate of the node.
         * @param y The global y-coordinate of the node.
         * @param z The global z-coordinate of the node.
         */
        void setGlobalCoord(const double & x, const double & y, const double & z);

        /**
         * Assign the final solved displacements (u,v,w) at the node.
         *
         * @param u The x-displacement at the node.
         * @param v The y-displacement at the node.
         * @param w The z-displacement at the node.
         */
        void setDisp(const double & u, const double & v, const double & w);

        /**
         * Assign the final solved force (Fx, Fy, Fz) at the node.
         *
         * @param Fx The x-force applied the node.
         * @param Fy The y-force applied the node.
         * @param Fz The z-force applied the node.
         */
        void setForce(const double & Fx, const double & Fy, const double & Fz);

        /**
         * Cumulate the stress and strain values at this node by summing up
         * the values from each adjacent element.
         *
         * @param strain The strain vector to be added to this node. The form is:
         * [e(x), e(y), e(z), gamma(xy), gamma(yz), gamma(zx))].
         * @param stress The stress vector to be added to this node. The form is:
         * [s(x), s(y), s(z), t(xy), t(yz), t(zx)].
         *
         * @note This is an internal step for computing the averaged strain and
         * stress at each node. The cumulative value will be averaged later.
         */
        void setStrainAndStress(const VectorXd & strain, const VectorXd & stress);

        /**
         * Average the strain vector at this node.
         *
         * @return The average strain vector.
         *
         * @note Node objects are dynamically allocated on heap memory, so we
         * can return-by-ref.
         */
        const VectorXd & averageStrain();

        /**
         * Average the stress vector at this node.
         *
         * @return The average stress vector.
         */
        const VectorXd & averageStress();

        /**
         * Get the index of this node.
         *
         * @return The index.
         */
        const int & getIndex() const;

        /**
         * Get the initial global coordinates (x,y,z) of the node.
         *
         * @return The global coordinates.
         */
        const Vector3d & getGlobalCoord() const;

        /**
         * Get the solved displacements (u,v,w) at the node.
         *
         * @return The displacements.
         */
        const Vector3d & getDisp() const;

        /**
         * Get the solved force (Fx, Fy, Fz) at the node.
         *
         * @return The forces.
         */
        const Vector3d & getForce() const;

    private:
        /** Zero-based index of the node */
        int index_;

        /** Global coordinates (x,y,z) */
        Vector3d globalCoord_;

        /** Displacements (u,v,w) */
        Vector3d disp_;

        /** Forces (Fx,Fy,Fz) */
        Vector3d force_;

        /** Strain vector [e(x), e(y), e(z), gamma(xy), gamma(yz), gamma(zx))] */
        VectorXd strain_;

        /** Strain vector [s(x), s(y), s(z), t(xy), t(yz), t(zx)] */
        VectorXd stress_;

        /** Number of adjacent elements at this node */
        int averageCount_;

        /**
         * Private helper function for deleting the current node, used in destructor
         * and assignment operator.
         */
        void clear_();

        /**
         * Private helper function for deep-copying the node, used in copy
         * constructor and assignment operator.
         */
        void copy_(Node const & other);

};

#endif /* Node_h */
