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

# Utility rule file for testdata_134.

# Include the progress variables for this target.
include CMakeFiles/testdata_134.dir/progress.make

CMakeFiles/testdata_134:
	/usr/bin/cmake -DAOM_CONFIG_DIR="/home/augusto/aom/build" -DAOM_ROOT="/home/augusto/aom" -DAOM_TEST_FILE="av1-1-b8-00-quantizer-45.ivf" -DAOM_TEST_CHECKSUM=0960bf018ada4224b8344519cf091850d50a57bd -P /home/augusto/aom/test/test_data_download_worker.cmake

testdata_134: CMakeFiles/testdata_134
testdata_134: CMakeFiles/testdata_134.dir/build.make

.PHONY : testdata_134

# Rule to build all files generated by this target.
CMakeFiles/testdata_134.dir/build: testdata_134

.PHONY : CMakeFiles/testdata_134.dir/build

CMakeFiles/testdata_134.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/testdata_134.dir/cmake_clean.cmake
.PHONY : CMakeFiles/testdata_134.dir/clean

CMakeFiles/testdata_134.dir/depend:
	cd /home/augusto/aom/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/augusto/aom /home/augusto/aom /home/augusto/aom/build /home/augusto/aom/build /home/augusto/aom/build/CMakeFiles/testdata_134.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/testdata_134.dir/depend

