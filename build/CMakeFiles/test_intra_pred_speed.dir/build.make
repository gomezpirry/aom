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

# Include any dependencies generated for this target.
include CMakeFiles/test_intra_pred_speed.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/test_intra_pred_speed.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/test_intra_pred_speed.dir/flags.make

CMakeFiles/test_intra_pred_speed.dir/gen_src/usage_exit.c.o: CMakeFiles/test_intra_pred_speed.dir/flags.make
CMakeFiles/test_intra_pred_speed.dir/gen_src/usage_exit.c.o: gen_src/usage_exit.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/augusto/aom/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/test_intra_pred_speed.dir/gen_src/usage_exit.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/test_intra_pred_speed.dir/gen_src/usage_exit.c.o   -c /home/augusto/aom/build/gen_src/usage_exit.c

CMakeFiles/test_intra_pred_speed.dir/gen_src/usage_exit.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/test_intra_pred_speed.dir/gen_src/usage_exit.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/augusto/aom/build/gen_src/usage_exit.c > CMakeFiles/test_intra_pred_speed.dir/gen_src/usage_exit.c.i

CMakeFiles/test_intra_pred_speed.dir/gen_src/usage_exit.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/test_intra_pred_speed.dir/gen_src/usage_exit.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/augusto/aom/build/gen_src/usage_exit.c -o CMakeFiles/test_intra_pred_speed.dir/gen_src/usage_exit.c.s

CMakeFiles/test_intra_pred_speed.dir/gen_src/usage_exit.c.o.requires:

.PHONY : CMakeFiles/test_intra_pred_speed.dir/gen_src/usage_exit.c.o.requires

CMakeFiles/test_intra_pred_speed.dir/gen_src/usage_exit.c.o.provides: CMakeFiles/test_intra_pred_speed.dir/gen_src/usage_exit.c.o.requires
	$(MAKE) -f CMakeFiles/test_intra_pred_speed.dir/build.make CMakeFiles/test_intra_pred_speed.dir/gen_src/usage_exit.c.o.provides.build
.PHONY : CMakeFiles/test_intra_pred_speed.dir/gen_src/usage_exit.c.o.provides

CMakeFiles/test_intra_pred_speed.dir/gen_src/usage_exit.c.o.provides.build: CMakeFiles/test_intra_pred_speed.dir/gen_src/usage_exit.c.o


CMakeFiles/test_intra_pred_speed.dir/test/test_intra_pred_speed.cc.o: CMakeFiles/test_intra_pred_speed.dir/flags.make
CMakeFiles/test_intra_pred_speed.dir/test/test_intra_pred_speed.cc.o: ../test/test_intra_pred_speed.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/augusto/aom/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/test_intra_pred_speed.dir/test/test_intra_pred_speed.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test_intra_pred_speed.dir/test/test_intra_pred_speed.cc.o -c /home/augusto/aom/test/test_intra_pred_speed.cc

CMakeFiles/test_intra_pred_speed.dir/test/test_intra_pred_speed.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_intra_pred_speed.dir/test/test_intra_pred_speed.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/augusto/aom/test/test_intra_pred_speed.cc > CMakeFiles/test_intra_pred_speed.dir/test/test_intra_pred_speed.cc.i

CMakeFiles/test_intra_pred_speed.dir/test/test_intra_pred_speed.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_intra_pred_speed.dir/test/test_intra_pred_speed.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/augusto/aom/test/test_intra_pred_speed.cc -o CMakeFiles/test_intra_pred_speed.dir/test/test_intra_pred_speed.cc.s

CMakeFiles/test_intra_pred_speed.dir/test/test_intra_pred_speed.cc.o.requires:

.PHONY : CMakeFiles/test_intra_pred_speed.dir/test/test_intra_pred_speed.cc.o.requires

CMakeFiles/test_intra_pred_speed.dir/test/test_intra_pred_speed.cc.o.provides: CMakeFiles/test_intra_pred_speed.dir/test/test_intra_pred_speed.cc.o.requires
	$(MAKE) -f CMakeFiles/test_intra_pred_speed.dir/build.make CMakeFiles/test_intra_pred_speed.dir/test/test_intra_pred_speed.cc.o.provides.build
.PHONY : CMakeFiles/test_intra_pred_speed.dir/test/test_intra_pred_speed.cc.o.provides

CMakeFiles/test_intra_pred_speed.dir/test/test_intra_pred_speed.cc.o.provides.build: CMakeFiles/test_intra_pred_speed.dir/test/test_intra_pred_speed.cc.o


# Object files for target test_intra_pred_speed
test_intra_pred_speed_OBJECTS = \
"CMakeFiles/test_intra_pred_speed.dir/gen_src/usage_exit.c.o" \
"CMakeFiles/test_intra_pred_speed.dir/test/test_intra_pred_speed.cc.o"

# External object files for target test_intra_pred_speed
test_intra_pred_speed_EXTERNAL_OBJECTS = \
"/home/augusto/aom/build/CMakeFiles/aom_common_app_util.dir/common/args.c.o" \
"/home/augusto/aom/build/CMakeFiles/aom_common_app_util.dir/common/av1_config.c.o" \
"/home/augusto/aom/build/CMakeFiles/aom_common_app_util.dir/common/md5_utils.c.o" \
"/home/augusto/aom/build/CMakeFiles/aom_common_app_util.dir/common/tools_common.c.o" \
"/home/augusto/aom/build/CMakeFiles/aom_common_app_util.dir/common/rawenc.c.o" \
"/home/augusto/aom/build/CMakeFiles/aom_common_app_util.dir/common/y4menc.c.o"

test_intra_pred_speed: CMakeFiles/test_intra_pred_speed.dir/gen_src/usage_exit.c.o
test_intra_pred_speed: CMakeFiles/test_intra_pred_speed.dir/test/test_intra_pred_speed.cc.o
test_intra_pred_speed: CMakeFiles/aom_common_app_util.dir/common/args.c.o
test_intra_pred_speed: CMakeFiles/aom_common_app_util.dir/common/av1_config.c.o
test_intra_pred_speed: CMakeFiles/aom_common_app_util.dir/common/md5_utils.c.o
test_intra_pred_speed: CMakeFiles/aom_common_app_util.dir/common/tools_common.c.o
test_intra_pred_speed: CMakeFiles/aom_common_app_util.dir/common/rawenc.c.o
test_intra_pred_speed: CMakeFiles/aom_common_app_util.dir/common/y4menc.c.o
test_intra_pred_speed: CMakeFiles/test_intra_pred_speed.dir/build.make
test_intra_pred_speed: libaom.a
test_intra_pred_speed: third_party/googletest/src/googletest/libgtest.a
test_intra_pred_speed: CMakeFiles/test_intra_pred_speed.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/augusto/aom/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable test_intra_pred_speed"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_intra_pred_speed.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/test_intra_pred_speed.dir/build: test_intra_pred_speed

.PHONY : CMakeFiles/test_intra_pred_speed.dir/build

CMakeFiles/test_intra_pred_speed.dir/requires: CMakeFiles/test_intra_pred_speed.dir/gen_src/usage_exit.c.o.requires
CMakeFiles/test_intra_pred_speed.dir/requires: CMakeFiles/test_intra_pred_speed.dir/test/test_intra_pred_speed.cc.o.requires

.PHONY : CMakeFiles/test_intra_pred_speed.dir/requires

CMakeFiles/test_intra_pred_speed.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/test_intra_pred_speed.dir/cmake_clean.cmake
.PHONY : CMakeFiles/test_intra_pred_speed.dir/clean

CMakeFiles/test_intra_pred_speed.dir/depend:
	cd /home/augusto/aom/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/augusto/aom /home/augusto/aom /home/augusto/aom/build /home/augusto/aom/build /home/augusto/aom/build/CMakeFiles/test_intra_pred_speed.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/test_intra_pred_speed.dir/depend

