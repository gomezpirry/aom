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
include CMakeFiles/aom_scale.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/aom_scale.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/aom_scale.dir/flags.make

CMakeFiles/aom_scale.dir/aom_scale/generic/aom_scale.c.o: CMakeFiles/aom_scale.dir/flags.make
CMakeFiles/aom_scale.dir/aom_scale/generic/aom_scale.c.o: ../aom_scale/generic/aom_scale.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/augusto/aom/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/aom_scale.dir/aom_scale/generic/aom_scale.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/aom_scale.dir/aom_scale/generic/aom_scale.c.o   -c /home/augusto/aom/aom_scale/generic/aom_scale.c

CMakeFiles/aom_scale.dir/aom_scale/generic/aom_scale.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/aom_scale.dir/aom_scale/generic/aom_scale.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/augusto/aom/aom_scale/generic/aom_scale.c > CMakeFiles/aom_scale.dir/aom_scale/generic/aom_scale.c.i

CMakeFiles/aom_scale.dir/aom_scale/generic/aom_scale.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/aom_scale.dir/aom_scale/generic/aom_scale.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/augusto/aom/aom_scale/generic/aom_scale.c -o CMakeFiles/aom_scale.dir/aom_scale/generic/aom_scale.c.s

CMakeFiles/aom_scale.dir/aom_scale/generic/aom_scale.c.o.requires:

.PHONY : CMakeFiles/aom_scale.dir/aom_scale/generic/aom_scale.c.o.requires

CMakeFiles/aom_scale.dir/aom_scale/generic/aom_scale.c.o.provides: CMakeFiles/aom_scale.dir/aom_scale/generic/aom_scale.c.o.requires
	$(MAKE) -f CMakeFiles/aom_scale.dir/build.make CMakeFiles/aom_scale.dir/aom_scale/generic/aom_scale.c.o.provides.build
.PHONY : CMakeFiles/aom_scale.dir/aom_scale/generic/aom_scale.c.o.provides

CMakeFiles/aom_scale.dir/aom_scale/generic/aom_scale.c.o.provides.build: CMakeFiles/aom_scale.dir/aom_scale/generic/aom_scale.c.o


CMakeFiles/aom_scale.dir/aom_scale/generic/gen_scalers.c.o: CMakeFiles/aom_scale.dir/flags.make
CMakeFiles/aom_scale.dir/aom_scale/generic/gen_scalers.c.o: ../aom_scale/generic/gen_scalers.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/augusto/aom/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/aom_scale.dir/aom_scale/generic/gen_scalers.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/aom_scale.dir/aom_scale/generic/gen_scalers.c.o   -c /home/augusto/aom/aom_scale/generic/gen_scalers.c

CMakeFiles/aom_scale.dir/aom_scale/generic/gen_scalers.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/aom_scale.dir/aom_scale/generic/gen_scalers.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/augusto/aom/aom_scale/generic/gen_scalers.c > CMakeFiles/aom_scale.dir/aom_scale/generic/gen_scalers.c.i

CMakeFiles/aom_scale.dir/aom_scale/generic/gen_scalers.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/aom_scale.dir/aom_scale/generic/gen_scalers.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/augusto/aom/aom_scale/generic/gen_scalers.c -o CMakeFiles/aom_scale.dir/aom_scale/generic/gen_scalers.c.s

CMakeFiles/aom_scale.dir/aom_scale/generic/gen_scalers.c.o.requires:

.PHONY : CMakeFiles/aom_scale.dir/aom_scale/generic/gen_scalers.c.o.requires

CMakeFiles/aom_scale.dir/aom_scale/generic/gen_scalers.c.o.provides: CMakeFiles/aom_scale.dir/aom_scale/generic/gen_scalers.c.o.requires
	$(MAKE) -f CMakeFiles/aom_scale.dir/build.make CMakeFiles/aom_scale.dir/aom_scale/generic/gen_scalers.c.o.provides.build
.PHONY : CMakeFiles/aom_scale.dir/aom_scale/generic/gen_scalers.c.o.provides

CMakeFiles/aom_scale.dir/aom_scale/generic/gen_scalers.c.o.provides.build: CMakeFiles/aom_scale.dir/aom_scale/generic/gen_scalers.c.o


CMakeFiles/aom_scale.dir/aom_scale/generic/yv12config.c.o: CMakeFiles/aom_scale.dir/flags.make
CMakeFiles/aom_scale.dir/aom_scale/generic/yv12config.c.o: ../aom_scale/generic/yv12config.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/augusto/aom/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object CMakeFiles/aom_scale.dir/aom_scale/generic/yv12config.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/aom_scale.dir/aom_scale/generic/yv12config.c.o   -c /home/augusto/aom/aom_scale/generic/yv12config.c

CMakeFiles/aom_scale.dir/aom_scale/generic/yv12config.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/aom_scale.dir/aom_scale/generic/yv12config.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/augusto/aom/aom_scale/generic/yv12config.c > CMakeFiles/aom_scale.dir/aom_scale/generic/yv12config.c.i

CMakeFiles/aom_scale.dir/aom_scale/generic/yv12config.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/aom_scale.dir/aom_scale/generic/yv12config.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/augusto/aom/aom_scale/generic/yv12config.c -o CMakeFiles/aom_scale.dir/aom_scale/generic/yv12config.c.s

CMakeFiles/aom_scale.dir/aom_scale/generic/yv12config.c.o.requires:

.PHONY : CMakeFiles/aom_scale.dir/aom_scale/generic/yv12config.c.o.requires

CMakeFiles/aom_scale.dir/aom_scale/generic/yv12config.c.o.provides: CMakeFiles/aom_scale.dir/aom_scale/generic/yv12config.c.o.requires
	$(MAKE) -f CMakeFiles/aom_scale.dir/build.make CMakeFiles/aom_scale.dir/aom_scale/generic/yv12config.c.o.provides.build
.PHONY : CMakeFiles/aom_scale.dir/aom_scale/generic/yv12config.c.o.provides

CMakeFiles/aom_scale.dir/aom_scale/generic/yv12config.c.o.provides.build: CMakeFiles/aom_scale.dir/aom_scale/generic/yv12config.c.o


CMakeFiles/aom_scale.dir/aom_scale/generic/yv12extend.c.o: CMakeFiles/aom_scale.dir/flags.make
CMakeFiles/aom_scale.dir/aom_scale/generic/yv12extend.c.o: ../aom_scale/generic/yv12extend.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/augusto/aom/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object CMakeFiles/aom_scale.dir/aom_scale/generic/yv12extend.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/aom_scale.dir/aom_scale/generic/yv12extend.c.o   -c /home/augusto/aom/aom_scale/generic/yv12extend.c

CMakeFiles/aom_scale.dir/aom_scale/generic/yv12extend.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/aom_scale.dir/aom_scale/generic/yv12extend.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/augusto/aom/aom_scale/generic/yv12extend.c > CMakeFiles/aom_scale.dir/aom_scale/generic/yv12extend.c.i

CMakeFiles/aom_scale.dir/aom_scale/generic/yv12extend.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/aom_scale.dir/aom_scale/generic/yv12extend.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/augusto/aom/aom_scale/generic/yv12extend.c -o CMakeFiles/aom_scale.dir/aom_scale/generic/yv12extend.c.s

CMakeFiles/aom_scale.dir/aom_scale/generic/yv12extend.c.o.requires:

.PHONY : CMakeFiles/aom_scale.dir/aom_scale/generic/yv12extend.c.o.requires

CMakeFiles/aom_scale.dir/aom_scale/generic/yv12extend.c.o.provides: CMakeFiles/aom_scale.dir/aom_scale/generic/yv12extend.c.o.requires
	$(MAKE) -f CMakeFiles/aom_scale.dir/build.make CMakeFiles/aom_scale.dir/aom_scale/generic/yv12extend.c.o.provides.build
.PHONY : CMakeFiles/aom_scale.dir/aom_scale/generic/yv12extend.c.o.provides

CMakeFiles/aom_scale.dir/aom_scale/generic/yv12extend.c.o.provides.build: CMakeFiles/aom_scale.dir/aom_scale/generic/yv12extend.c.o


aom_scale: CMakeFiles/aom_scale.dir/aom_scale/generic/aom_scale.c.o
aom_scale: CMakeFiles/aom_scale.dir/aom_scale/generic/gen_scalers.c.o
aom_scale: CMakeFiles/aom_scale.dir/aom_scale/generic/yv12config.c.o
aom_scale: CMakeFiles/aom_scale.dir/aom_scale/generic/yv12extend.c.o
aom_scale: CMakeFiles/aom_scale.dir/build.make

.PHONY : aom_scale

# Rule to build all files generated by this target.
CMakeFiles/aom_scale.dir/build: aom_scale

.PHONY : CMakeFiles/aom_scale.dir/build

CMakeFiles/aom_scale.dir/requires: CMakeFiles/aom_scale.dir/aom_scale/generic/aom_scale.c.o.requires
CMakeFiles/aom_scale.dir/requires: CMakeFiles/aom_scale.dir/aom_scale/generic/gen_scalers.c.o.requires
CMakeFiles/aom_scale.dir/requires: CMakeFiles/aom_scale.dir/aom_scale/generic/yv12config.c.o.requires
CMakeFiles/aom_scale.dir/requires: CMakeFiles/aom_scale.dir/aom_scale/generic/yv12extend.c.o.requires

.PHONY : CMakeFiles/aom_scale.dir/requires

CMakeFiles/aom_scale.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/aom_scale.dir/cmake_clean.cmake
.PHONY : CMakeFiles/aom_scale.dir/clean

CMakeFiles/aom_scale.dir/depend:
	cd /home/augusto/aom/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/augusto/aom /home/augusto/aom /home/augusto/aom/build /home/augusto/aom/build /home/augusto/aom/build/CMakeFiles/aom_scale.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/aom_scale.dir/depend

