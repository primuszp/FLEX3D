# C-FLEX2D
2D Axisymmetric Finite Element Methods (FEM) Program for Flexible Pavement Analysis.

The program is developed by Haohang Huang and Jiayi Luo @ University of Illinois at Urbana-Champaign (UIUC). The advisor is [Professor Erol Tutumluer](https://cee.illinois.edu/directory/profile/tutumlue).

The development of this software program is sponsored by [Engineer Research and Development Center (ERDC)](https://www.erdc.usace.army.mil), United States Army Corps of Engineers (USACE). The ERDC project manager is [Jeb S. Tingle](https://www.erdc.usace.army.mil/Media/Images/igphoto/2002117275/).

The program name represents Corps of Engineers ("**C**") and Flexible pavements ("**FLEX**").

The code is in C++ and Python. 2D FEM engine is written in C++. Other related modules, extensions, and Graphical User Interface (GUI) are written in Python.

## Installation Guide
Clone the repository:
`git clone https://github.com/symphonylyh/C-FLEX2D`

Repo strcuture:
```bash
C-FLEX2D
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

5. Right click project--Properties--General--Target Name--Change to the name you like, e.g. `main2d`

6. Right click project--Build. You'll see some complaints from the compiler, but it doesn't matter.

The executable `main2d` will be under `./build/x64/Release/`

### MacOS and Linux
For Unix systems, it is much easier to compile. CMake and Xcode command line tools should be installed. Then simply run the bash script
`sh compile.sh`
and the executable will be under directory `./build/FEM/main`
