# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.17.1/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.17.1/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/huang/Downloads/CS73-gdcp

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/huang/Downloads/CS73-gdcp/build

# Include any dependencies generated for this target.
include CMakeFiles/a4.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/a4.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/a4.dir/flags.make

CMakeFiles/a4.dir/src/a4_main.cpp.o: CMakeFiles/a4.dir/flags.make
CMakeFiles/a4.dir/src/a4_main.cpp.o: ../src/a4_main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/huang/Downloads/CS73-gdcp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/a4.dir/src/a4_main.cpp.o"
	/usr/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/a4.dir/src/a4_main.cpp.o -c /Users/huang/Downloads/CS73-gdcp/src/a4_main.cpp

CMakeFiles/a4.dir/src/a4_main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/a4.dir/src/a4_main.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/huang/Downloads/CS73-gdcp/src/a4_main.cpp > CMakeFiles/a4.dir/src/a4_main.cpp.i

CMakeFiles/a4.dir/src/a4_main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/a4.dir/src/a4_main.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/huang/Downloads/CS73-gdcp/src/a4_main.cpp -o CMakeFiles/a4.dir/src/a4_main.cpp.s

CMakeFiles/a4.dir/src/filtering.cpp.o: CMakeFiles/a4.dir/flags.make
CMakeFiles/a4.dir/src/filtering.cpp.o: ../src/filtering.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/huang/Downloads/CS73-gdcp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/a4.dir/src/filtering.cpp.o"
	/usr/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/a4.dir/src/filtering.cpp.o -c /Users/huang/Downloads/CS73-gdcp/src/filtering.cpp

CMakeFiles/a4.dir/src/filtering.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/a4.dir/src/filtering.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/huang/Downloads/CS73-gdcp/src/filtering.cpp > CMakeFiles/a4.dir/src/filtering.cpp.i

CMakeFiles/a4.dir/src/filtering.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/a4.dir/src/filtering.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/huang/Downloads/CS73-gdcp/src/filtering.cpp -o CMakeFiles/a4.dir/src/filtering.cpp.s

CMakeFiles/a4.dir/src/a2.cpp.o: CMakeFiles/a4.dir/flags.make
CMakeFiles/a4.dir/src/a2.cpp.o: ../src/a2.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/huang/Downloads/CS73-gdcp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/a4.dir/src/a2.cpp.o"
	/usr/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/a4.dir/src/a2.cpp.o -c /Users/huang/Downloads/CS73-gdcp/src/a2.cpp

CMakeFiles/a4.dir/src/a2.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/a4.dir/src/a2.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/huang/Downloads/CS73-gdcp/src/a2.cpp > CMakeFiles/a4.dir/src/a2.cpp.i

CMakeFiles/a4.dir/src/a2.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/a4.dir/src/a2.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/huang/Downloads/CS73-gdcp/src/a2.cpp -o CMakeFiles/a4.dir/src/a2.cpp.s

# Object files for target a4
a4_OBJECTS = \
"CMakeFiles/a4.dir/src/a4_main.cpp.o" \
"CMakeFiles/a4.dir/src/filtering.cpp.o" \
"CMakeFiles/a4.dir/src/a2.cpp.o"

# External object files for target a4
a4_EXTERNAL_OBJECTS =

a4: CMakeFiles/a4.dir/src/a4_main.cpp.o
a4: CMakeFiles/a4.dir/src/filtering.cpp.o
a4: CMakeFiles/a4.dir/src/a2.cpp.o
a4: CMakeFiles/a4.dir/build.make
a4: libcommon_lib.a
a4: CMakeFiles/a4.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/huang/Downloads/CS73-gdcp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable a4"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/a4.dir/link.txt --verbose=$(VERBOSE)
	/usr/local/Cellar/cmake/3.17.1/bin/cmake -E make_directory /Users/huang/Downloads/CS73-gdcp/data/output

# Rule to build all files generated by this target.
CMakeFiles/a4.dir/build: a4

.PHONY : CMakeFiles/a4.dir/build

CMakeFiles/a4.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/a4.dir/cmake_clean.cmake
.PHONY : CMakeFiles/a4.dir/clean

CMakeFiles/a4.dir/depend:
	cd /Users/huang/Downloads/CS73-gdcp/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/huang/Downloads/CS73-gdcp /Users/huang/Downloads/CS73-gdcp /Users/huang/Downloads/CS73-gdcp/build /Users/huang/Downloads/CS73-gdcp/build /Users/huang/Downloads/CS73-gdcp/build/CMakeFiles/a4.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/a4.dir/depend

