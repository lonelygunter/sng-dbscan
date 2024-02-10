BUILD = g++ -std=c++11 -ligraph
BUILDFLAGS = -Wall -Wextra

# Source files
SOURCES = main.cpp headers/*.cpp

# Header files
HEADERS = headers/*.h

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
