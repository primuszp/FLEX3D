/**
 * @file main.cpp
 * Main execution interface.
 *
 * @author Haohang Huang
 * @date May 22, 2019
 */

#include <iostream>
#include <fstream>
#include <string>
#include <chrono>

#include "Linear.h"
#include "Nonlinear.h"

// 2D-->3D notes:
// 1. Vector2d to Vector3d, VectorXd should be treated carefully in the initializer list
// 2. stress & strain tensor dimension 4-->6
// 3. More input parameters for constructor
// Edit Flow: Node --> Shape --> Element --> Analysis --> Material --> Mesh

int main(int argc, char const *argv[]) {
    //--------------------------------------------------------------------------
    //--------------------------------Main program------------------------------
    //--------------------------------------------------------------------------

auto start = std::chrono::high_resolution_clock::now();

bool xcode = false;
    if (!xcode) {
        // Command line version
        for (int i = 1; i < argc; i++) {
            std::string inFile(argv[i]);
            std::string inFileName = inFile + ".txt";
            std::string outVTKName = inFile + ".vtk";

            Mesh mesh = Mesh(inFileName); // on stack, make sure lifetime of 'mesh' is longer than Analysis case

            Analysis* caseType; // 'case' is a reserved keyword for switch(), so can't declare a variable with name "case"
			if (mesh.nonlinear) {
				std::cout << "> Nonlinear analysis scheme" << std::endl;
				caseType = new Nonlinear(mesh);
			}
			else {
				std::cout << "> Linear analysis scheme" << std::endl;
				caseType = new Linear(mesh);
			}

            caseType->solve();
            // caseType->printDisp();
            // caseType->printStrain();
            // caseType->printStress();
            caseType->writeToVTK(outVTKName);
            delete caseType; caseType = NULL;
        }
    } else {
        // Xcode version
        // Designate file name
        std::string inFileName = "../test_input.txt";
        std::string outVTKName = "../test_output.vtk";
        //
        // std::cout << "Enter the input file name: ";
        // std::string inFileName;
        // std::getline(std::cin, inFileName);
        // inFileName += ".txt";
        //
        // std::cout << "Enter the output VTK file name: ";
        // std::string outFileName;
        // std::getline(std::cin, outFileName);
        // std::string outVTKName = outFileName + ".vtk";
        // outFileName += ".txt";

        Mesh mesh = Mesh(inFileName);
        Analysis* caseType;
        if (mesh.nonlinear)
            caseType = new Nonlinear(mesh);
        else
            caseType = new Linear(mesh);

        caseType->solve();
        // caseType->printDisp();
        // caseType->printStrain();
        // caseType->printStress();
        caseType->writeToVTK(outVTKName);
        delete caseType; caseType = NULL;
    }

// auto finish = std::chrono::high_resolution_clock::now();
// auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
// std::chrono::duration<double> elapsed = finish - start; // in second
// std::cout << "Elapsed time total: " << elapsed.count() << " ms" << std::endl;

// std::cout << "-----------------------------------------------------";
// std::cout << std::endl;
// std::cout << "Mission Completed! Thanks for using our program!" << std::endl;

    return 0;
}
