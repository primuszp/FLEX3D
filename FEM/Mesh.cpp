/**
 * @file Mesh.cpp
 * Implementation of Mesh class.
 *
 * @author Haohang Huang
 * @date May 21, 2019
 */

#include "Mesh.h"
#include "ElementB20.h"
#include "LinearElastic.h"
#include "NonlinearElastic.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

Mesh::Mesh()
{
}

Mesh::Mesh(std::string const & fileName)
{
    readFromFile(fileName);
}

Mesh::~Mesh()
{
    // Delete all nodes
    for (int i = 0; i < nodeCount_; i++) {
        delete meshNode_[i]; meshNode_[i] = NULL;
    }
    delete[] meshNode_; meshNode_ = NULL;

    // Delete all elements
    for (int i = 0; i < elementCount_; i++) {
        delete meshElement_[i]; meshElement_[i] = NULL;
    }
    delete[] meshElement_; meshElement_ = NULL;

    // Delete all materials
    for (auto & m : materialList) {
        delete m; m = NULL;
    }

}

void Mesh::readFromFile(std::string const & fileName)
{
    std::ifstream file(fileName); // open file as read-only stream // same as std::ifstream file; file.open(fileName); // ifstream is read-only type
    std::string readLine; // for each single line
    if (!file) {
        std::cerr << "ERROR: Cannot open file... Please check your file name or file path. Aborting." << std::endl;
        exit(-1);
    } // error handling

    // -------------------------------------------------------------------------
    // ------------------------ General Information ----------------------------
    // -------------------------------------------------------------------------
    // Read the first line of the file as mesh summary info including:
    // [0] Total number of nodes
    // [1] Total number of elements
    // [2] Total number of different types of material property the elements should be assigned to
    // [3] Total number of point loads
    // [4] Total number of edge loads
    // [5] Total number of face loads
    // [6] Total number of displacement constraints in X-direction
    // [7] Total number of displacement constraints in Y-direction
    // [8] Total number of displacement constraints in Z-direction
    std::vector<int> meshSummary;
    std::getline(file, readLine);
    parseLine(readLine, meshSummary);
    nodeCount_ = meshSummary[0];
    elementCount_ = meshSummary[1];
    int elementProperties = meshSummary[2];
    int pointLoad = meshSummary[3];
    int edgeLoad = meshSummary[4];
    int faceLoad = meshSummary[5];
    int boundaryX = meshSummary[6];
    int boundaryY = meshSummary[7];
    int boundaryZ = meshSummary[8];
    std::vector<int>().swap(meshSummary); // enforce to free the vector memory

    // Create node & element array on HEAP
    meshNode_ = new Node*[nodeCount_];
    meshElement_ = new Element*[elementCount_];

    // -------------------------------------------------------------------------
    // ------------------------ Material Properties ----------------------------
    // -------------------------------------------------------------------------
    // Read element properties (modulus, Poisson's ratio, thermal parameters, etc)
    // Format:
    // Line 1: start & end index of element, material isotropy and linearity. e.g., 0 35 0 0 means element No.0~35 are isotropic and linear elastic material.
    // Line 2 and so forth: specific material properties
    std::map<int, int> layerMap; // use an ordered map to decide layer No. by range finding, e.g., N layers, store N start indices of the element as [0 N1) [N1 N2) [N2 N3)
    std::vector<double> elementProperty;
    materialList.reserve(elementProperties);
    nonlinear = false;
    for (int i = 0; i < elementProperties; i++) {
        std::getline(file, readLine);
        std::vector<int> range;
        parseLine(readLine, range); // range will store 0 19 0 0
        layerMap[range[0]] = i; // record the lower bound of the range
        std::getline(file, readLine);
        parseLine(readLine, elementProperty);
        if (range[3] == 0) // linear elastic
            materialList.push_back(new LinearElastic(range[2], range[3], elementProperty)); // dynamically allocated, remember to delete in destructor!
        else { // nonlinear elastic
            nonlinear = true;
            // for nonlinear layer it should read two more lines about the material model parameters
            std::getline(file, readLine);
            int model = std::stoi(readLine, NULL); // [K-theta:1 Uzan:2 UT-Austin:3 MEPDG:4]
            std::getline(file, readLine);
            std::vector<double> parameters;
            parseLine(readLine, parameters);
            materialList.push_back(new NonlinearElastic(range[2], range[3], elementProperty, model, parameters));
        }
        // 0 if isotropic, 1 if cross-anisotropic; 0 if linear elastic, 1 if nonlinear elastic
    }
    std::vector<double>().swap(elementProperty);

    // For nonlinear analysis it should read one additional line after all layers, for iteration parameters (No. of load increments, damping ratios) if there are at least one nonlinear material
    if (nonlinear) {
        std::getline(file, readLine);
        parseLine(readLine, iterations); // vector "iterations" declared as a member variable of Mesh class, to be used in ctor of Nonlinear class
    }

    // -------------------------------------------------------------------------
    // -------------------- Loading (Point & Edge & Face) ----------------------
    // -------------------------------------------------------------------------
    loadNodeList.clear();
    loadValue.clear(); // these two are the public member variables of Mesh class, to be used in Analysis->applyForce()
    // Read point loads
    if (pointLoad > 0) {
        std::vector<int> pointNodeList;
        std::vector<double> pointLoadValue;
        std::getline(file, readLine);
        parseLine(readLine, pointNodeList);
        std::getline(file, readLine);
        parseLine(readLine, pointLoadValue);
        for (int i = 0; i < pointLoad; i++) {
            // Convert node index to DOF index
            loadNodeList.push_back(3 * pointNodeList[i]);
            loadNodeList.push_back(3 * pointNodeList[i] + 1);
            loadNodeList.push_back(3 * pointNodeList[i] + 2);
            // Assign load value
            loadValue.push_back(pointLoadValue[3 * i]);
            loadValue.push_back(pointLoadValue[3 * i + 1]);
            loadValue.push_back(pointLoadValue[3 * i + 2]);
        }
        std::vector<int>().swap(pointNodeList);
        std::vector<double>().swap(pointLoadValue);
    }

    edgeLoadElementList.clear();
    loadEdgeList.clear();
    edgeLoadValue.clear();
    // Read edge loads
    if (edgeLoad > 0) {
        std::vector<int> line1;
        std::vector<int> edgeNo;
        std::vector<double> line2;
        for (int i = 0; i < edgeLoad; i++) {
            std::getline(file, readLine);
            parseLine(readLine, line1);
            edgeLoadElementList.emplace_back(line1.front()); // 1st entry: element index
            edgeNo.assign(line1.begin() + 1, line1.end()); // rest entries: edge list // Note: assign() will rewrite edgeNo. using erase() is also ok
            loadEdgeList.emplace_back(edgeNo);

            std::getline(file, readLine);
            parseLine(readLine, line2);
            edgeLoadValue.emplace_back(line2);
        }
    }

    faceLoadElementList.clear();
    loadFaceList.clear();
    faceLoadValue.clear();
    // Read face loads
    if (faceLoad > 0) {
        std::vector<int> line1;
        std::vector<int> faceNo;
        std::vector<double> line2;
        for (int i = 0; i < faceLoad; i++) {
            std::getline(file, readLine);
            parseLine(readLine, line1);
            faceLoadElementList.emplace_back(line1.front()); // 1st entry: element index
            faceNo.assign(line1.begin() + 1, line1.end()); // rest entries: face list // Note: assign() will rewrite edgeNo
            loadFaceList.emplace_back(faceNo);

            std::getline(file, readLine);
            parseLine(readLine, line2);
            faceLoadValue.emplace_back(line2);
        }
    }

    // -------------------------------------------------------------------------
    // ------------------------ Node Coordinates -------------------------------
    // -------------------------------------------------------------------------
    // Read node coordinates
    std::vector<double> nodeCoord(3); // 3D
    for (int i = 0; i < nodeCount_; i++) {
        std::getline(file, readLine);
        parseLine(readLine, nodeCoord);
        meshNode_[i] = new Node(i, nodeCoord[0], nodeCoord[1], nodeCoord[2]);
    }
    std::vector<double>().swap(nodeCoord);

    // -------------------------------------------------------------------------
    // ------------------------ Element Indices --------------------------------
    // -------------------------------------------------------------------------
    // Read element's node list and create the corresponding element type
    std::vector<int> elementNodeList;
    for (int i = 0; i < elementCount_; i++) {
        std::getline(file, readLine);
        // Option 2: denote element node number at the beginning of the line
        std::string::size_type j = readLine.find(' ', 0); // find the first space
        int size = std::stoi(readLine.substr(0, j));
        readLine.erase(0, j + 1);
        parseLine(readLine, elementNodeList);
        // Find the material type of the current element
        // std::map<int, int>::iterator it = layerMap.lower_bound(i);
        // Material* material = materialList[it->second];
        // (Solved) @BUG // lower_bound will give the included index, upper_bound is non-included...WRONG! lower_bound will give the first no-less-than element! Not the real "lower bound" as we assumed. Should use upper_bound() - 1
        // Fix:
        std::map<int, int>::iterator it = layerMap.upper_bound(i); // upper_bound gives the first key that will go AFTER i
        Material* material = materialList[(--it)->second];
        // Create instances of different types of element
        switch (size) {
            case 20 :
                meshElement_[i] = new ElementB20(i, elementNodeList, meshNode_, material);
                break;
        }
    }
    std::vector<int>().swap(elementNodeList);

    // -------------------------------------------------------------------------
    // ------------------------ Boundary Conditions ----------------------------
    // -------------------------------------------------------------------------
    boundaryNodeList.clear();
    boundaryValue.clear(); // these two are the public member variables of Mesh class, to be used in Analysis->assembleStiffness()
    // Read and set X-direction boundary condition (no need to set to Node at this time, the nodal displacement will be finalized after sovling Ku=F in Analysis->solveDisp())
    if (boundaryX > 0) {
        std::vector<int> boundaryXNodeList;
        std::vector<double> boundaryXValue;
        std::getline(file, readLine);
        parseLine(readLine, boundaryXNodeList);
        std::getline(file, readLine);
        parseLine(readLine, boundaryXValue);
        for (int i = 0; i < boundaryX; i++) {
            boundaryNodeList.push_back(3 * boundaryXNodeList[i]); // convert to DOF index
            boundaryValue.push_back(boundaryXValue[i]);
        }
        std::vector<int>().swap(boundaryXNodeList);
        std::vector<double>().swap(boundaryXValue);
    }

    // Read and set Y-direction boundary condition
    if (boundaryY > 0) {
        std::vector<int> boundaryYNodeList;
        std::vector<double> boundaryYValue;
        std::getline(file, readLine);
        parseLine(readLine, boundaryYNodeList);
        std::getline(file, readLine);
        parseLine(readLine, boundaryYValue);
        for (int i = 0; i < boundaryY; i++) {
            boundaryNodeList.push_back(3 * boundaryYNodeList[i] + 1);
            boundaryValue.push_back(boundaryYValue[i]);
        }
        std::vector<int>().swap(boundaryYNodeList);
        std::vector<double>().swap(boundaryYValue);
    }

    // Read and set Z-direction boundary condition
    if (boundaryY > 0) {
        std::vector<int> boundaryZNodeList;
        std::vector<double> boundaryZValue;
        std::getline(file, readLine);
        parseLine(readLine, boundaryZNodeList);
        std::getline(file, readLine);
        parseLine(readLine, boundaryZValue);
        for (int i = 0; i < boundaryZ; i++) {
            boundaryNodeList.push_back(3 * boundaryZNodeList[i] + 2);
            boundaryValue.push_back(boundaryZValue[i]);
        }
        std::vector<int>().swap(boundaryZNodeList);
        std::vector<double>().swap(boundaryZValue);
    }

    // Complete read-in
    file.close();
}

const int & Mesh::nodeCount() const
{
    return nodeCount_;
}

const int & Mesh::elementCount() const
{
    return elementCount_;
}

Node* Mesh::getNode(const int & index) const
{
    return meshNode_[index];
}

Element* Mesh::getElement(const int & index) const
{
    return meshElement_[index];
}

Node** Mesh::nodeArray() const
{
    return meshNode_;
}

Element** Mesh::elementArray() const
{
    return meshElement_;
}

template<typename T>
void Mesh::parseLine(std::string const & readLine, std::vector<T> & parseLine) const
{
    // Idea:
    // use stringstream to parse the string (separate with " " by default) into segments
    // or if we want to specify the delimiter, use string.find() as in option 2

    std::stringstream ss(readLine);
    T segment;
    parseLine.clear();
    while (ss >> segment)
        parseLine.push_back(segment);

    // Option 2: a more general string segmentation approach
    // void parse(const std::string & s, char delimiter, std::vector<std::string> & v) {
    //     // Note:
    //     // 1. std::string.find(delimiter, start) finds the given "delimiter" in a given range ["start", end]
    //     // 2. std::string.substr(start, length) gives the substring begins at "start" and spans "length". "length" can be oversize
    //     // 3. std::string::npos is a special sign for "not found", typically a very large number
    //     std::string::size_type i = 0; // at head
    //     std::string::size_type j = s.find(delimiter, 0); // find the location of first delimiter
    //     while (j != std::string::npos) { // find until the last delimiter
    //         v.push_back(s.substr(i, j - i)); // parse the substring segment before the delimiter
    //         i = ++j; // move start point to one-pass the found delimiter
    //         j = s.find(delimiter, i); // find next delimiter in the string left
    //     }
    //     if(j == std::string::npos) // if finally no delimiter (reach the end), read in the last segment (an oversized length is ok). Here we use s.length() b/c when there's no delimiter we can just read the whole string.
    //         v.push_back(s.substr(i, s.length()));
    // }
    // // Usage
    // std::string str = "Hello World!";
    // std::vector<std::string> v;
    // parse(str,' ',v); // use space as delimiter
    // for(auto & e : v)
    //     std::cout << e << " "; // output is exactly the same
}
