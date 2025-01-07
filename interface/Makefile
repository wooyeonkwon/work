# Makefile example for ROOT-based project

# Set compiler and flags
CXX = g++
CXXFLAGS = `root-config --cflags`

# Set the library output name
LIBRARY_NAME = libAnalysisClasses.so

# Set source files
SOURCES = $(wildcard *.cc)

# Set object files
OBJECTS = $(SOURCES:.cc=.o)

# Set the target for building the library
all: $(LIBRARY_NAME)

$(LIBRARY_NAME): $(OBJECTS)
	$(CXX) -shared -o $@ $^ $(shell root-config --libs)


# Target for cleaning up object files and libraries
clean:
		rm -f $(OBJECTS) $(LIBRARY_NAME)
