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

# Utility rule file for fftwBuild.

# Include the progress variables for this target.
include src/contrib/fftw/CMakeFiles/fftwBuild.dir/progress.make

src/contrib/fftw/CMakeFiles/fftwBuild: src/contrib/fftw/CMakeFiles/fftwBuild-complete

src/contrib/fftw/CMakeFiles/fftwBuild-complete: src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-install
src/contrib/fftw/CMakeFiles/fftwBuild-complete: src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-mkdir
src/contrib/fftw/CMakeFiles/fftwBuild-complete: src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-download
src/contrib/fftw/CMakeFiles/fftwBuild-complete: src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-update
src/contrib/fftw/CMakeFiles/fftwBuild-complete: src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-patch
src/contrib/fftw/CMakeFiles/fftwBuild-complete: src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-configure
src/contrib/fftw/CMakeFiles/fftwBuild-complete: src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-build
src/contrib/fftw/CMakeFiles/fftwBuild-complete: src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-install
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Completed 'fftwBuild'"
	cd /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw && "/Applications/CMake 2.8-12.app/Contents/bin/cmake" -E make_directory /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/CMakeFiles
	cd /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw && "/Applications/CMake 2.8-12.app/Contents/bin/cmake" -E touch /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/CMakeFiles/fftwBuild-complete
	cd /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw && "/Applications/CMake 2.8-12.app/Contents/bin/cmake" -E touch /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-done

src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-install: src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-build
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Performing install step for 'fftwBuild'"
	cd /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-build && $(MAKE) install
	cd /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-build && "/Applications/CMake 2.8-12.app/Contents/bin/cmake" -E touch /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-install

src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-mkdir:
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Creating directories for 'fftwBuild'"
	cd /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw && "/Applications/CMake 2.8-12.app/Contents/bin/cmake" -E make_directory /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/fftwBuild-prefix/src/fftwBuild
	cd /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw && "/Applications/CMake 2.8-12.app/Contents/bin/cmake" -E make_directory /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-build
	cd /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw && "/Applications/CMake 2.8-12.app/Contents/bin/cmake" -E make_directory /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/fftwBuild-prefix
	cd /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw && "/Applications/CMake 2.8-12.app/Contents/bin/cmake" -E make_directory /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/fftwBuild-prefix/tmp
	cd /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw && "/Applications/CMake 2.8-12.app/Contents/bin/cmake" -E make_directory /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp
	cd /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw && "/Applications/CMake 2.8-12.app/Contents/bin/cmake" -E make_directory /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/fftwBuild-prefix/src
	cd /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw && "/Applications/CMake 2.8-12.app/Contents/bin/cmake" -E touch /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-mkdir

src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-download: src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-urlinfo.txt
src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-download: src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-mkdir
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Performing download step (download, verify and extract) for 'fftwBuild'"
	cd /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/fftwBuild-prefix/src && "/Applications/CMake 2.8-12.app/Contents/bin/cmake" -P /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/download-fftwBuild.cmake
	cd /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/fftwBuild-prefix/src && "/Applications/CMake 2.8-12.app/Contents/bin/cmake" -P /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/verify-fftwBuild.cmake
	cd /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/fftwBuild-prefix/src && "/Applications/CMake 2.8-12.app/Contents/bin/cmake" -P /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/extract-fftwBuild.cmake
	cd /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/fftwBuild-prefix/src && "/Applications/CMake 2.8-12.app/Contents/bin/cmake" -E touch /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-download

src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-update: src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-download
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/CMakeFiles $(CMAKE_PROGRESS_5)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "No update step for 'fftwBuild'"
	cd /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw && "/Applications/CMake 2.8-12.app/Contents/bin/cmake" -E touch /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-update

src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-patch: src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-download
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/CMakeFiles $(CMAKE_PROGRESS_6)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "No patch step for 'fftwBuild'"
	cd /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw && "/Applications/CMake 2.8-12.app/Contents/bin/cmake" -E touch /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-patch

src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-configure: src/contrib/fftw/fftwBuild-prefix/tmp/fftwBuild-cfgcmd.txt
src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-configure: src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-update
src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-configure: src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-patch
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/CMakeFiles $(CMAKE_PROGRESS_7)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Performing configure step for 'fftwBuild'"
	cd /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-build && /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/fftwBuild-prefix/src/fftwBuild/configure --prefix=/Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/fftwBuild-prefix --libdir=/Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/fftwBuild-prefix/lib --disable-shared --enable-static --enable-sse2 --enable-float
	cd /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-build && "/Applications/CMake 2.8-12.app/Contents/bin/cmake" -E touch /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-configure

src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-build: src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-configure
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/CMakeFiles $(CMAKE_PROGRESS_8)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Performing build step for 'fftwBuild'"
	cd /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-build && $(MAKE)
	cd /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-build && "/Applications/CMake 2.8-12.app/Contents/bin/cmake" -E touch /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-build

fftwBuild: src/contrib/fftw/CMakeFiles/fftwBuild
fftwBuild: src/contrib/fftw/CMakeFiles/fftwBuild-complete
fftwBuild: src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-install
fftwBuild: src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-mkdir
fftwBuild: src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-download
fftwBuild: src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-update
fftwBuild: src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-patch
fftwBuild: src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-configure
fftwBuild: src/contrib/fftw/fftwBuild-prefix/src/fftwBuild-stamp/fftwBuild-build
fftwBuild: src/contrib/fftw/CMakeFiles/fftwBuild.dir/build.make
.PHONY : fftwBuild

# Rule to build all files generated by this target.
src/contrib/fftw/CMakeFiles/fftwBuild.dir/build: fftwBuild
.PHONY : src/contrib/fftw/CMakeFiles/fftwBuild.dir/build

src/contrib/fftw/CMakeFiles/fftwBuild.dir/clean:
	cd /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw && $(CMAKE_COMMAND) -P CMakeFiles/fftwBuild.dir/cmake_clean.cmake
.PHONY : src/contrib/fftw/CMakeFiles/fftwBuild.dir/clean

src/contrib/fftw/CMakeFiles/fftwBuild.dir/depend:
	cd /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5 /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/src/contrib/fftw /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw /Users/er113/2014-05_beth_pcbm-cg-dos/gromacs-4.6.5/build/src/contrib/fftw/CMakeFiles/fftwBuild.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/contrib/fftw/CMakeFiles/fftwBuild.dir/depend
