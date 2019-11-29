# FLEX3D
3D Finite Element Methods (FEM) Program for Flexible Pavement Analysis.

The program is developed by Haohang Huang @ University of Illinois at Urbana-Champaign (UIUC).

The program name represents Flexible pavements ("**FLEX**").

The code is in C++.

## Installation Guide
Clone the repository:
`git clone https://github.com/symphonylyh/FLEX3D`

Repo strcuture:
```bash
FLEX3D
├── CMakeLists.txt
├── FEM
│   └── Eigen
├── LICENSE
├── README.md
└── compile.sh
```

Note: `CMakeLists.txt` and `compile.sh` are only used on MacOS or Linux platform.

### Windows (Tested on Windows 10, Visual Studio 2017)
For Windows platform, it is recommended to compile the program by Microsoft Visual Studio.
1. File--New--Project--Visual C++--Windows Desktop--Windows console application--Name & Location. Set name to "build" and location to the path of this folder you just cloned. Leave "Create directory for solution" ticked. Remove the starter files VS automatically generated (they are usually `pch.cpp`, `pch.h`, etc.)

2. Tool bar--Change "Solution Configuration" to "Release", "Solution platform" to "x64". "Release" program is faster than "Debug"'s.

3. Tool bar--Project--Add existing items--Choose all files in folder "FEM"

4. Left Pane--Right click on project name "build"--Properties--Right click on the project name-Properties-C/C++-Precompiled headers-Precompiles header--Not using precompile header

5. Right click project--Properties--General--Target Name--Change to the name you like, e.g. `main3d`

6. Right click project--Build. You'll see some complaints from the compiler, but it doesn't matter.

The executable `main3d` will be under `./build/x64/Release/`

### MacOS and Linux
For Unix systems, it is much easier to compile. CMake and Xcode command line tools should be installed. Then simply run the bash script
`sh compile.sh`
and the executable will be under directory `./build/FEM/main`
