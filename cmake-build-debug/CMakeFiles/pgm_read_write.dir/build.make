# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.12

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "D:\CLion 2018.2\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "D:\CLion 2018.2\bin\cmake\win\bin\cmake.exe" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:\Users\Administrator\Desktop\sift_cpp

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\Users\Administrator\Desktop\sift_cpp\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/pgm_read_write.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/pgm_read_write.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/pgm_read_write.dir/flags.make

CMakeFiles/pgm_read_write.dir/main.cpp.obj: CMakeFiles/pgm_read_write.dir/flags.make
CMakeFiles/pgm_read_write.dir/main.cpp.obj: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\Administrator\Desktop\sift_cpp\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/pgm_read_write.dir/main.cpp.obj"
	D:\MinGW64\bin\g++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\pgm_read_write.dir\main.cpp.obj -c C:\Users\Administrator\Desktop\sift_cpp\main.cpp

CMakeFiles/pgm_read_write.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pgm_read_write.dir/main.cpp.i"
	D:\MinGW64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\Administrator\Desktop\sift_cpp\main.cpp > CMakeFiles\pgm_read_write.dir\main.cpp.i

CMakeFiles/pgm_read_write.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pgm_read_write.dir/main.cpp.s"
	D:\MinGW64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\Administrator\Desktop\sift_cpp\main.cpp -o CMakeFiles\pgm_read_write.dir\main.cpp.s

CMakeFiles/pgm_read_write.dir/utils.cpp.obj: CMakeFiles/pgm_read_write.dir/flags.make
CMakeFiles/pgm_read_write.dir/utils.cpp.obj: ../utils.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\Administrator\Desktop\sift_cpp\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/pgm_read_write.dir/utils.cpp.obj"
	D:\MinGW64\bin\g++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\pgm_read_write.dir\utils.cpp.obj -c C:\Users\Administrator\Desktop\sift_cpp\utils.cpp

CMakeFiles/pgm_read_write.dir/utils.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pgm_read_write.dir/utils.cpp.i"
	D:\MinGW64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\Administrator\Desktop\sift_cpp\utils.cpp > CMakeFiles\pgm_read_write.dir\utils.cpp.i

CMakeFiles/pgm_read_write.dir/utils.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pgm_read_write.dir/utils.cpp.s"
	D:\MinGW64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\Administrator\Desktop\sift_cpp\utils.cpp -o CMakeFiles\pgm_read_write.dir\utils.cpp.s

# Object files for target pgm_read_write
pgm_read_write_OBJECTS = \
"CMakeFiles/pgm_read_write.dir/main.cpp.obj" \
"CMakeFiles/pgm_read_write.dir/utils.cpp.obj"

# External object files for target pgm_read_write
pgm_read_write_EXTERNAL_OBJECTS =

pgm_read_write.exe: CMakeFiles/pgm_read_write.dir/main.cpp.obj
pgm_read_write.exe: CMakeFiles/pgm_read_write.dir/utils.cpp.obj
pgm_read_write.exe: CMakeFiles/pgm_read_write.dir/build.make
pgm_read_write.exe: C:\Users\Administrator\Desktop\pgm_read_write\img.a
pgm_read_write.exe: CMakeFiles/pgm_read_write.dir/linklibs.rsp
pgm_read_write.exe: CMakeFiles/pgm_read_write.dir/objects1.rsp
pgm_read_write.exe: CMakeFiles/pgm_read_write.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=C:\Users\Administrator\Desktop\sift_cpp\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable pgm_read_write.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\pgm_read_write.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/pgm_read_write.dir/build: pgm_read_write.exe

.PHONY : CMakeFiles/pgm_read_write.dir/build

CMakeFiles/pgm_read_write.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\pgm_read_write.dir\cmake_clean.cmake
.PHONY : CMakeFiles/pgm_read_write.dir/clean

CMakeFiles/pgm_read_write.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" C:\Users\Administrator\Desktop\sift_cpp C:\Users\Administrator\Desktop\sift_cpp C:\Users\Administrator\Desktop\sift_cpp\cmake-build-debug C:\Users\Administrator\Desktop\sift_cpp\cmake-build-debug C:\Users\Administrator\Desktop\sift_cpp\cmake-build-debug\CMakeFiles\pgm_read_write.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/pgm_read_write.dir/depend

