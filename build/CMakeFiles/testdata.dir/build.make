# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/augusto/aom

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/augusto/aom/build

# Utility rule file for testdata.

# Include the progress variables for this target.
include CMakeFiles/testdata.dir/progress.make

testdata: CMakeFiles/testdata.dir/build.make

.PHONY : testdata

# Rule to build all files generated by this target.
CMakeFiles/testdata.dir/build: testdata

.PHONY : CMakeFiles/testdata.dir/build

CMakeFiles/testdata.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/testdata.dir/cmake_clean.cmake
.PHONY : CMakeFiles/testdata.dir/clean

CMakeFiles/testdata.dir/depend:
	cd /home/augusto/aom/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/augusto/aom /home/augusto/aom /home/augusto/aom/build /home/augusto/aom/build /home/augusto/aom/build/CMakeFiles/testdata.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/testdata.dir/depend

