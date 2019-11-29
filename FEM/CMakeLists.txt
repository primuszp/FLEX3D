# This is the 2-level cmake file for each source code folder.

# Source code from current folder
aux_source_directory(. SRC_LIST)

# Executable
add_executable(main ${SRC_LIST})

# Source code from subfolders
# add_subdirectory(Eigen) // since Eigen is a header-only library, you don't need to do any compilation for it. The CMakeLists.txt under "Eigen" folder is useless

# Link libraries (the library object name should match the one of subfolder)
# target_link_libraries(main Eigen) // since Eigen is a header-only library, you don't need to do any compilation for it. TODO: how to link the compiled Eigen library image?

# Compiler flags
add_definitions(-O3 -Wall -Wextra -Wno-tautological-compare -Wno-sign-compare -Wmissing-variable-declarations -Wmissing-declarations -Wno-error=unused-const-variable -Wno-error=unused-parameter -std=c++11)
# -O3 for compiler optimization
