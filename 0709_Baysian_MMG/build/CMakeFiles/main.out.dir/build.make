# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /home/black/Desktop/sakanashi/0709_Baysian_MMG

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/black/Desktop/sakanashi/0709_Baysian_MMG/build

# Include any dependencies generated for this target.
include CMakeFiles/main.out.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/main.out.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/main.out.dir/flags.make

CMakeFiles/main.out.dir/cmaes.f90.o: CMakeFiles/main.out.dir/flags.make
CMakeFiles/main.out.dir/cmaes.f90.o: ../cmaes.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/black/Desktop/sakanashi/0709_Baysian_MMG/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object CMakeFiles/main.out.dir/cmaes.f90.o"
	/usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/black/Desktop/sakanashi/0709_Baysian_MMG/cmaes.f90 -o CMakeFiles/main.out.dir/cmaes.f90.o

CMakeFiles/main.out.dir/cmaes.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/main.out.dir/cmaes.f90.i"
	/usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/black/Desktop/sakanashi/0709_Baysian_MMG/cmaes.f90 > CMakeFiles/main.out.dir/cmaes.f90.i

CMakeFiles/main.out.dir/cmaes.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/main.out.dir/cmaes.f90.s"
	/usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/black/Desktop/sakanashi/0709_Baysian_MMG/cmaes.f90 -o CMakeFiles/main.out.dir/cmaes.f90.s

CMakeFiles/main.out.dir/eigen.f90.o: CMakeFiles/main.out.dir/flags.make
CMakeFiles/main.out.dir/eigen.f90.o: ../eigen.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/black/Desktop/sakanashi/0709_Baysian_MMG/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building Fortran object CMakeFiles/main.out.dir/eigen.f90.o"
	/usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/black/Desktop/sakanashi/0709_Baysian_MMG/eigen.f90 -o CMakeFiles/main.out.dir/eigen.f90.o

CMakeFiles/main.out.dir/eigen.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/main.out.dir/eigen.f90.i"
	/usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/black/Desktop/sakanashi/0709_Baysian_MMG/eigen.f90 > CMakeFiles/main.out.dir/eigen.f90.i

CMakeFiles/main.out.dir/eigen.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/main.out.dir/eigen.f90.s"
	/usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/black/Desktop/sakanashi/0709_Baysian_MMG/eigen.f90 -o CMakeFiles/main.out.dir/eigen.f90.s

CMakeFiles/main.out.dir/random.f90.o: CMakeFiles/main.out.dir/flags.make
CMakeFiles/main.out.dir/random.f90.o: ../random.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/black/Desktop/sakanashi/0709_Baysian_MMG/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building Fortran object CMakeFiles/main.out.dir/random.f90.o"
	/usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/black/Desktop/sakanashi/0709_Baysian_MMG/random.f90 -o CMakeFiles/main.out.dir/random.f90.o

CMakeFiles/main.out.dir/random.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/main.out.dir/random.f90.i"
	/usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/black/Desktop/sakanashi/0709_Baysian_MMG/random.f90 > CMakeFiles/main.out.dir/random.f90.i

CMakeFiles/main.out.dir/random.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/main.out.dir/random.f90.s"
	/usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/black/Desktop/sakanashi/0709_Baysian_MMG/random.f90 -o CMakeFiles/main.out.dir/random.f90.s

CMakeFiles/main.out.dir/problem_SI.f90.o: CMakeFiles/main.out.dir/flags.make
CMakeFiles/main.out.dir/problem_SI.f90.o: ../problem_SI.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/black/Desktop/sakanashi/0709_Baysian_MMG/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building Fortran object CMakeFiles/main.out.dir/problem_SI.f90.o"
	/usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/black/Desktop/sakanashi/0709_Baysian_MMG/problem_SI.f90 -o CMakeFiles/main.out.dir/problem_SI.f90.o

CMakeFiles/main.out.dir/problem_SI.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/main.out.dir/problem_SI.f90.i"
	/usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/black/Desktop/sakanashi/0709_Baysian_MMG/problem_SI.f90 > CMakeFiles/main.out.dir/problem_SI.f90.i

CMakeFiles/main.out.dir/problem_SI.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/main.out.dir/problem_SI.f90.s"
	/usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/black/Desktop/sakanashi/0709_Baysian_MMG/problem_SI.f90 -o CMakeFiles/main.out.dir/problem_SI.f90.s

CMakeFiles/main.out.dir/MMG_model.f90.o: CMakeFiles/main.out.dir/flags.make
CMakeFiles/main.out.dir/MMG_model.f90.o: ../MMG_model.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/black/Desktop/sakanashi/0709_Baysian_MMG/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building Fortran object CMakeFiles/main.out.dir/MMG_model.f90.o"
	/usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/black/Desktop/sakanashi/0709_Baysian_MMG/MMG_model.f90 -o CMakeFiles/main.out.dir/MMG_model.f90.o

CMakeFiles/main.out.dir/MMG_model.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/main.out.dir/MMG_model.f90.i"
	/usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/black/Desktop/sakanashi/0709_Baysian_MMG/MMG_model.f90 > CMakeFiles/main.out.dir/MMG_model.f90.i

CMakeFiles/main.out.dir/MMG_model.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/main.out.dir/MMG_model.f90.s"
	/usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/black/Desktop/sakanashi/0709_Baysian_MMG/MMG_model.f90 -o CMakeFiles/main.out.dir/MMG_model.f90.s

CMakeFiles/main.out.dir/main.f90.o: CMakeFiles/main.out.dir/flags.make
CMakeFiles/main.out.dir/main.f90.o: ../main.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/black/Desktop/sakanashi/0709_Baysian_MMG/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building Fortran object CMakeFiles/main.out.dir/main.f90.o"
	/usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/black/Desktop/sakanashi/0709_Baysian_MMG/main.f90 -o CMakeFiles/main.out.dir/main.f90.o

CMakeFiles/main.out.dir/main.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/main.out.dir/main.f90.i"
	/usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/black/Desktop/sakanashi/0709_Baysian_MMG/main.f90 > CMakeFiles/main.out.dir/main.f90.i

CMakeFiles/main.out.dir/main.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/main.out.dir/main.f90.s"
	/usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/black/Desktop/sakanashi/0709_Baysian_MMG/main.f90 -o CMakeFiles/main.out.dir/main.f90.s

# Object files for target main.out
main_out_OBJECTS = \
"CMakeFiles/main.out.dir/cmaes.f90.o" \
"CMakeFiles/main.out.dir/eigen.f90.o" \
"CMakeFiles/main.out.dir/random.f90.o" \
"CMakeFiles/main.out.dir/problem_SI.f90.o" \
"CMakeFiles/main.out.dir/MMG_model.f90.o" \
"CMakeFiles/main.out.dir/main.f90.o"

# External object files for target main.out
main_out_EXTERNAL_OBJECTS =

main.out: CMakeFiles/main.out.dir/cmaes.f90.o
main.out: CMakeFiles/main.out.dir/eigen.f90.o
main.out: CMakeFiles/main.out.dir/random.f90.o
main.out: CMakeFiles/main.out.dir/problem_SI.f90.o
main.out: CMakeFiles/main.out.dir/MMG_model.f90.o
main.out: CMakeFiles/main.out.dir/main.f90.o
main.out: CMakeFiles/main.out.dir/build.make
main.out: /usr/lib/gcc/x86_64-linux-gnu/9/libgomp.so
main.out: /usr/lib/x86_64-linux-gnu/libpthread.so
main.out: CMakeFiles/main.out.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/black/Desktop/sakanashi/0709_Baysian_MMG/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Linking Fortran executable main.out"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/main.out.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/main.out.dir/build: main.out

.PHONY : CMakeFiles/main.out.dir/build

CMakeFiles/main.out.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/main.out.dir/cmake_clean.cmake
.PHONY : CMakeFiles/main.out.dir/clean

CMakeFiles/main.out.dir/depend:
	cd /home/black/Desktop/sakanashi/0709_Baysian_MMG/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/black/Desktop/sakanashi/0709_Baysian_MMG /home/black/Desktop/sakanashi/0709_Baysian_MMG /home/black/Desktop/sakanashi/0709_Baysian_MMG/build /home/black/Desktop/sakanashi/0709_Baysian_MMG/build /home/black/Desktop/sakanashi/0709_Baysian_MMG/build/CMakeFiles/main.out.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/main.out.dir/depend

