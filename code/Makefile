MAXIMUM_MATCHING_DIR := tests/maximum-matching
PERFECT_MATCHING_DIR := tests/perfect-matching
MATRIX_INVERSE_DIR := tests/matrix-inverse
ALG ?= 2

.PHONY: all clear

all: test_matrix_inverse test_perfect_matching test_maximum_matching

test_matrix_inverse:
	@echo "Testing Matrix Inverse algorithms..."
	@cd $(MATRIX_INVERSE_DIR) && \
		$(MAKE) && \
		$(MAKE) run && \
		$(MAKE) compare && \
		$(MAKE) clean

test_perfect_matching:
	@echo "Testing Perfect Matching algorithms..."
	@cd $(PERFECT_MATCHING_DIR) && \
		$(MAKE) && \
		$(MAKE) run "arg=$(ALG)" && \
		$(MAKE) compare && \
		$(MAKE) clean

test_maximum_matching:
	@echo "Testing Maximum Matching algorithms..."
	@cd $(MAXIMUM_MATCHING_DIR) && \
		$(MAKE) && \
		$(MAKE) run "arg=$(ALG)" && \
		$(MAKE) compare && \
		$(MAKE) clean

clean:
	@cd $(MATRIX_INVERSE_DIR) && \
		$(MAKE) clean
	@cd $(PERFECT_MATCHING_DIR) && \
		$(MAKE) clean
	@cd $(MAXIMUM_MATCHING_DIR) && \
		$(MAKE) clean
