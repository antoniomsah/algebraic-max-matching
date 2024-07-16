# Variables
CXX := g++
CXXFLAGS := -Wall -Wextra -std=c++17
TARGET := matrix-generator 
SRC := matrix-generator.cpp
INPUT_SRC := ../input

# Rules
all: $(TARGET)

.PHONY: $(TARGET)
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

run: $(OUT)
	./$(TARGET) $(arg) 10 > $(INPUT_SRC)/$(arg)-random0
	./$(TARGET) $(arg) 100 > $(INPUT_SRC)/$(arg)-random1
	./$(TARGET) $(arg) 1000 > $(INPUT_SRC)/$(arg)-random2
	./$(TARGET) $(arg) 1000 > $(INPUT_SRC)/$(arg)-random3
	./$(TARGET) $(arg) 1000000 > $(INPUT_SRC)/$(arg)-random4
	./$(TARGET) $(arg) 1000000 > $(INPUT_SRC)/$(arg)-random5
	./$(TARGET) $(arg) 1000000 > $(INPUT_SRC)/$(arg)-random6
	./$(TARGET) $(arg) 1000000000 > $(INPUT_SRC)/$(arg)-random7
	./$(TARGET) $(arg) 1000000000 > $(INPUT_SRC)/$(arg)-random8
	./$(TARGET) $(arg) 1000000000 > $(INPUT_SRC)/$(arg)-random9

clean:
	rm -f $(TARGET)

# Phony targets
.PHONY: all clean
