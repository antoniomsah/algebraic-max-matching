# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -O3

# Source file
SRC = main.cpp

# Output executable
OUT = main

# Input and output directories
INPUT_DIR = input
OUTPUT_DIR = output
JUDGE_DIR = judge

# Files
OUTPUT_FILES = $(wildcard $(OUTPUT_DIR)/*)

# Default target
all: $(OUT)

# Build the executable
$(OUT): $(SRC)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Run the program with each file in the input directory
run: $(OUT)
	mkdir -p $(OUTPUT_DIR)
	for file in $(INPUT_DIR)/*; do \
		basefile=$$(basename $$file); \
		echo "starting $$basefile"; \
		./$(OUT) < $$file > $(OUTPUT_DIR)/$$basefile; \
		echo "concluded $$basefile"; \
	done

# Compares files from output with judge
compare: 
	@for output_file in $(OUTPUT_FILES); do \
		judge_file=$(JUDGE_DIR)/$$(basename $$output_file); \
		if [ -f $$judge_file ]; then \
			if ! diff $$output_file $$judge_file > /dev/null; then \
				echo "Failed $$judge_file"; \
				diff $$output_file $$judge_file; \
				exit 1; \
			fi; \
		else \
				echo "Judge file $$judge_file not found"; \
		fi \
	done

# Clean up the generated files
clean:
	rm -f $(OUT) $(OUTPUT_DIR)/*

.PHONY: compare
