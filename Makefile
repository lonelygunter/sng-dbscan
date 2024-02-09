BUILD = g++ -std=c++11
BUILDFLAGS = -Wall -Wextra

# Source files
SOURCES = main.cpp headers/readDataset.cpp headers/connectedComponents.cpp headers/compareWsnSample.cpp headers/loadingBar.cpp

# Header files
HEADERS = headers/readDataset.h headers/connectedComponents.h headers/compareWsnSample.h headers/loadingBar.h

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
