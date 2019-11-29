/**
 * @file Mesh.h
 * Mesh class for the entire mesh information read from input file.
 *
 * @author Haohang Huang
 * @date May 21, 2019
 */

#ifndef Mesh_h
#define Mesh_h

#include "Node.h"
#include "Element.h"
#include "Material.h"
#include <vector>

/* Mesh class for storing the node, element, material, boundary, load
 * information read from input file.
 */
class Mesh
{
    public:
        /**
         * Default constructor for Mesh.
         */
        Mesh();

        /**
         * Custom constructor to create a mesh with given file name.
         *
         * @param fileName The path to the input file.
         */
        Mesh(std::string const & fileName);

        /**
         * Destructor for Mesh.
         */
        ~Mesh();

        /**
         * Read in all the information and build the mesh.
         *
         * @param fileName The path to the input file.
         */
        void readFromFile(std::string const & fileName);

        /**
         * Get the total number of the nodes.
         *
         * @return The total number of the nodes.
         */
        const int & nodeCount() const;

        /**
         * Get the total number of the elements.
         *
         * @return The total number of the elements.
         */
        const int & elementCount() const;

        /**
         * Get the node at given zero-based index.
         *
         * @param index The index for query.
         * @return A pointer to the node being queried.
         */
        Node* getNode(const int & index) const;

        /**
         * Get the element at given zero-based index.
         *
         * @param index The index for query.
         * @return A pointer to the element being queried.
         */
        Element* getElement(const int & index) const;

        /**
         * Get the node array.
         *
         * @return A pointer to the node pool.
         */
        Node** nodeArray() const;

        /**
         * Get the element array.
         *
         * @return A pointer to the element pool.
         *
         * @note There are two ways to query a node/element:
         * mesh.nodeArray()[i] OR mesh.getNode(i)
         */
        Element** elementArray() const;

        /** A list of layered materials */
        std::vector<Material*> materialList;

        /** Analysis type */
        bool nonlinear;

        /** Iteration parameters used in nonlinear scheme */
        std::vector<double> iterations;

        /** A DOF list (see above) of the node with applied point load */
        std::vector<int> loadNodeList;

        /** The point load values at the node */
        std::vector<double> loadValue;

        /** A list of the elements with applied edge load */
        std::vector<int> edgeLoadElementList;

        /** A list of loaded edges in each edge-loaded element above */
        std::vector<std::vector<int> > loadEdgeList;

        /** A list of the edge load values (Fx, Fy, Fz) for each edge-loaded element above */
        std::vector<std::vector<double> > edgeLoadValue;

        /** A list of the elements with applied face load */
        std::vector<int> faceLoadElementList;

        /** A list of loaded faces in each face-loaded element above */
        std::vector<std::vector<int> > loadFaceList;

        /** A list of the face load values (Fx, Fy, Fz) for each face-loaded element above */
        std::vector<std::vector<double> > faceLoadValue;

        /** A list of the degree of freedom that are at boundary
         *
         * @note Actually this is not NODE list, but degree of freedom list,
         * e.g., if node 1,2 are X&Y fixed, node 4 is fully fixed, then the DOF list
         * is (3,4 | 6,7 | 12,13,14) and the boundary value list below is (0,0,0,0,0,0,0).
         */
        std::vector<int> boundaryNodeList;

        /** The boundary values corresponding to the DOF list above */
        std::vector<double> boundaryValue;

    private:
        /** Total number of nodes */
        int nodeCount_;

        /** Total number of elements */
        int elementCount_;

        /** A pointer to the node pool
         *
         * @note In old version, meshNode_ is a arr<Node> instead of arr<Node*>,
         * so when I do meshNode_ = new Node[size]; meshNode_[i] = Node(...),
         * actually the Node ctor is first called and a copy assignment operator
         * is called, which is redundant.
         */
        Node** meshNode_;

        /** A pointer to the element pool */
        Element** meshElement_;

        /**
         * A templated private helper function for readFromFile().
         *
         * Parse a line in the input file and place the values into an array.
         *
         * @param readLine The line to be parsed.
         * @param parseLine The array to store the parsed result.
         */
        template<typename T>
        void parseLine(std::string const & readLine, std::vector<T> & parseLine) const;
};

#endif /* Mesh_h */
