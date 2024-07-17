# Variables
CXX := g++
CXXFLAGS := -Wall -Wextra -std=c++17
TARGET := graph-generator 
SRC := graph-generator.cpp
INPUT_SRC := ../input

# Rules
all: $(TARGET)

.PHONY: $(TARGET)
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

run: $(OUT)
	./$(TARGET) $(arg) 0.10 > $(INPUT_SRC)/$(arg)-random0
	./$(TARGET) $(arg) 0.20 > $(INPUT_SRC)/$(arg)-random1
	./$(TARGET) $(arg) 0.25 > $(INPUT_SRC)/$(arg)-random2
	./$(TARGET) $(arg) 0.50 > $(INPUT_SRC)/$(arg)-random3
	./$(TARGET) $(arg) 0.50 > $(INPUT_SRC)/$(arg)-random4
	./$(TARGET) $(arg) 0.50 > $(INPUT_SRC)/$(arg)-random5
	./$(TARGET) $(arg) 0.75 > $(INPUT_SRC)/$(arg)-random6
	./$(TARGET) $(arg) 0.75 > $(INPUT_SRC)/$(arg)-random7
	./$(TARGET) $(arg) 0.75 > $(INPUT_SRC)/$(arg)-random8
	./$(TARGET) $(arg) 1.00 > $(INPUT_SRC)/$(arg)-random9

clean:
	rm -f $(TARGET)

# Phony targets
.PHONY: all clean
