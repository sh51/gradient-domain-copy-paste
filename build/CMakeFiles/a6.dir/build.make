# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.18.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.18.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/ryantucker/cs73/gradient-domain-copy-paste

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/ryantucker/cs73/gradient-domain-copy-paste/build

# Include any dependencies generated for this target.
include CMakeFiles/a6.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/a6.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/a6.dir/flags.make

CMakeFiles/a6.dir/src/a6_main.cpp.o: CMakeFiles/a6.dir/flags.make
CMakeFiles/a6.dir/src/a6_main.cpp.o: ../src/a6_main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/ryantucker/cs73/gradient-domain-copy-paste/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/a6.dir/src/a6_main.cpp.o"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/a6.dir/src/a6_main.cpp.o -c /Users/ryantucker/cs73/gradient-domain-copy-paste/src/a6_main.cpp

CMakeFiles/a6.dir/src/a6_main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/a6.dir/src/a6_main.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/ryantucker/cs73/gradient-domain-copy-paste/src/a6_main.cpp > CMakeFiles/a6.dir/src/a6_main.cpp.i

CMakeFiles/a6.dir/src/a6_main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/a6.dir/src/a6_main.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/ryantucker/cs73/gradient-domain-copy-paste/src/a6_main.cpp -o CMakeFiles/a6.dir/src/a6_main.cpp.s

# Object files for target a6
a6_OBJECTS = \
"CMakeFiles/a6.dir/src/a6_main.cpp.o"

# External object files for target a6
a6_EXTERNAL_OBJECTS =

a6: CMakeFiles/a6.dir/src/a6_main.cpp.o
a6: CMakeFiles/a6.dir/build.make
a6: libcommon_lib.a
a6: CMakeFiles/a6.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/ryantucker/cs73/gradient-domain-copy-paste/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable a6"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/a6.dir/link.txt --verbose=$(VERBOSE)
	/usr/local/Cellar/cmake/3.18.2/bin/cmake -E make_directory /Users/ryantucker/cs73/gradient-domain-copy-paste/data/output

# Rule to build all files generated by this target.
CMakeFiles/a6.dir/build: a6

.PHONY : CMakeFiles/a6.dir/build

CMakeFiles/a6.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/a6.dir/cmake_clean.cmake
.PHONY : CMakeFiles/a6.dir/clean

CMakeFiles/a6.dir/depend:
	cd /Users/ryantucker/cs73/gradient-domain-copy-paste/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/ryantucker/cs73/gradient-domain-copy-paste /Users/ryantucker/cs73/gradient-domain-copy-paste /Users/ryantucker/cs73/gradient-domain-copy-paste/build /Users/ryantucker/cs73/gradient-domain-copy-paste/build /Users/ryantucker/cs73/gradient-domain-copy-paste/build/CMakeFiles/a6.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/a6.dir/depend

