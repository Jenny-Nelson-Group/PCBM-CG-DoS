# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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
CMAKE_COMMAND = "/Applications/CMake 2.8-12.app/Contents/bin/cmake"

# The command to remove a file.
RM = "/Applications/CMake 2.8-12.app/Contents/bin/cmake" -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = "/Applications/CMake 2.8-12.app/Contents/bin/ccmake"

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build

# Include any dependencies generated for this target.
include share/template/CMakeFiles/template.dir/depend.make

# Include the progress variables for this target.
include share/template/CMakeFiles/template.dir/progress.make

# Include the compile flags for this target's objects.
include share/template/CMakeFiles/template.dir/flags.make

share/template/CMakeFiles/template.dir/template.c.o: share/template/CMakeFiles/template.dir/flags.make
share/template/CMakeFiles/template.dir/template.c.o: ../share/template/template.c
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object share/template/CMakeFiles/template.dir/template.c.o"
	cd /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/share/template && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/template.dir/template.c.o   -c /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/share/template/template.c

share/template/CMakeFiles/template.dir/template.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/template.dir/template.c.i"
	cd /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/share/template && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/share/template/template.c > CMakeFiles/template.dir/template.c.i

share/template/CMakeFiles/template.dir/template.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/template.dir/template.c.s"
	cd /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/share/template && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/share/template/template.c -o CMakeFiles/template.dir/template.c.s

share/template/CMakeFiles/template.dir/template.c.o.requires:
.PHONY : share/template/CMakeFiles/template.dir/template.c.o.requires

share/template/CMakeFiles/template.dir/template.c.o.provides: share/template/CMakeFiles/template.dir/template.c.o.requires
	$(MAKE) -f share/template/CMakeFiles/template.dir/build.make share/template/CMakeFiles/template.dir/template.c.o.provides.build
.PHONY : share/template/CMakeFiles/template.dir/template.c.o.provides

share/template/CMakeFiles/template.dir/template.c.o.provides.build: share/template/CMakeFiles/template.dir/template.c.o

# Object files for target template
template_OBJECTS = \
"CMakeFiles/template.dir/template.c.o"

# External object files for target template
template_EXTERNAL_OBJECTS =

share/template/template: share/template/CMakeFiles/template.dir/template.c.o
share/template/template: share/template/CMakeFiles/template.dir/build.make
share/template/template: src/gmxlib/libgmx.8.dylib
share/template/template: //Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/fftwBuild-prefix/lib/libfftw3f.a
share/template/template: share/template/CMakeFiles/template.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable template"
	cd /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/share/template && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/template.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
share/template/CMakeFiles/template.dir/build: share/template/template
.PHONY : share/template/CMakeFiles/template.dir/build

share/template/CMakeFiles/template.dir/requires: share/template/CMakeFiles/template.dir/template.c.o.requires
.PHONY : share/template/CMakeFiles/template.dir/requires

share/template/CMakeFiles/template.dir/clean:
	cd /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/share/template && $(CMAKE_COMMAND) -P CMakeFiles/template.dir/cmake_clean.cmake
.PHONY : share/template/CMakeFiles/template.dir/clean

share/template/CMakeFiles/template.dir/depend:
	cd /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5 /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/share/template /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/share/template /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/share/template/CMakeFiles/template.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : share/template/CMakeFiles/template.dir/depend

