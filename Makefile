BUILD = g++ -std=c++11 -ligraph
BUILDFLAGS = -Wall -Wextra

# Source files
SOURCES = main.cpp headers/readDataset.cpp headers/loadingBar.cpp headers/connectedComponents.cpp headers/conn_comp_igraph.cpp headers/compareWsnSample.cpp

# Header files
HEADERS = headers/readDataset.h headers/loadingBar.h headers/connectedComponents.h headers/conn_comp_igraph.h headers/compareWsnSample.h

# Object files
OBJECTS = $(SOURCES:.cpp=.o)

# Executable name
EXECUTABLE = main

# Default target
all: $(EXECUTABLE)

# Rule to create the executable
$(EXECUTABLE): $(OBJECTS)
	$(BUILD) $(BUILDFLAGS) -o $@ $^

# Rule to compile source files
%.o: %.cpp $(HEADERS)
	$(BUILD) $(BUILDFLAGS) -c -o $@ $<

# Clean target
clean:
	@echo "Cleaning .o files"
	rm -f $(EXECUTABLE) $(OBJECTS)
