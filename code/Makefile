# Variables
CXX := g++
CXXFLAGS := -Wall -g -Wextra -std=c++17
TARGET := main
SRC := main.cpp

# Rules
all: $(TARGET)

.PHONY: $(TARGET)
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

clean:
	rm -f $(TARGET)

# Phony targets
.PHONY: all clean
