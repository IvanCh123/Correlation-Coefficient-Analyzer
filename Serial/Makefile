# levdist makefile v1.0


# Variables =====================================

# Compiler options
CC=cc
CXX=g++
FLAGS=-Wall -Wextra -pthread -std=gnu11 -fopenmp
CFLAGS=$(FLAGS)
CXXFLAGS=$(FLAGS)
LIBS=

# Configure flags according to the target
debug: FLAGS += -g
release: FLAGS += -O3 -DNDEBUG

# Directories
BIN_DIR=bin
BUILD_DIR=build
DOC_DIR=doc
SRC_DIR=src
TEST_A_DIR=test
TEST_U_DIR=test_u
TEST_A_OUT_DIR=$(BUILD_DIR)/$(TEST_A_DIR)
TEST_U_OUT_DIR=$(BUILD_DIR)/$(TEST_U_DIR)
INST_DIR=$(HOME)/bin

# Project variables
APP_NAME=corr
#EXECUTABLE=$(BIN_DIR)/$(shell basename `pwd`)
EXECUTABLE=$(BIN_DIR)/$(APP_NAME)
FILTER_MEMCHECK=| grep 'lost'

# Files and test cases
HEADERS=$(wildcard $(SRC_DIR)/*.h)
SOURCES=$(wildcard $(SRC_DIR)/*.c*)
OBJECTS=$(SOURCES:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)

# Test cases ASCII and Unicode
TESTS_A=$(wildcard $(TEST_A_DIR)/*)
TESTS_A_OUT=$(TESTS_A:$(TEST_A_DIR)/%=$(TEST_A_OUT_DIR)/%)
TESTS_A_MEMCHECK=$(TESTS_A:$(TEST_A_DIR)/%=$(TEST_A_OUT_DIR)/%_memcheck)
TESTS_U=$(wildcard $(TEST_U_DIR)/*)
TESTS_U_OUT=$(TESTS_U:$(TEST_U_DIR)/%=$(TEST_U_OUT_DIR)/%)
TESTS_U_MEMCHECK=$(TESTS_U:$(TEST_U_DIR)/%=$(TEST_U_OUT_DIR)/%_memcheck)


# Building ======================================

# Default rule is debug target
debug: $(EXECUTABLE)
release: $(EXECUTABLE)

# All targets
all: debug doc test

# Linker call
$(EXECUTABLE): $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(CFLAGS) -o $@ $^ $(LIBS)

# Compile C source into .o object file
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c $(HEADERS) | $(BUILD_DIR)
	$(CC) -c $(CFLAGS) $< -o $@

# Compile C++ source into .o object file
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS) | $(BUILD_DIR)
	$(CXX) -c $(CXXFLAGS) $< -o $@

# Documentation
.PHONY: doc | $(DOC_DIR)
doc:
	doxygen Doxyfile


# Testing =======================================

# Function to test a test case directory
# $(call test_dir,in_dir,out_dir)
test_dir = @mkdir -p $2 ;\
	cd $(1) ;\
	../../$(EXECUTABLE) $(shell cat $(1)/args.txt) < input.txt > ../../$(2)/output.txt 2> ../../$(2)/error.txt || true ;\
	cd ../.. ;\
	$(HOME)/bin/icdiff $(1)/output.txt $(2)/output.txt || true ;\
	$(HOME)/bin/icdiff $(1)/error.txt $(2)/error.txt || true ;\
	rm -rf $(2) ;\
	$(eval testcnt=$(shell echo $$(($(testcnt)+1))))

# Testing in ASCII
.PHONY: test
test: debug $(TESTS_A_OUT) | $(TEST_A_OUT_DIR)
	@echo "$(testcnt) ascii cases tested, for details call 'make test -n'"

.PHONY: $(TEST_A_OUT_DIR)/%
$(TEST_A_OUT_DIR)/%: $(TEST_A_DIR)/% | icdiff
	$(call test_dir,$<,$@)

# Testing in Unicode
.PHONY: test_u
test_u: debug $(TESTS_U_OUT) | $(TEST_U_OUT_DIR)
	@echo "$(testcnt) unicode cases tested, for details call 'make test_u -n'"

.PHONY: $(TEST_U_OUT_DIR)/%
$(TEST_U_OUT_DIR)/%: $(TEST_U_DIR)/% | icdiff
	$(call test_dir,$<,$@)


# Profiling =====================================

# Function to profile a test case directory
# $(call memcheck_dir,in_dir)
memcheck_dir = 	@cd $(1) ;\
	echo "cd $(1) && valgrind ../../$(EXECUTABLE) -Q $(shell cat $(1)/args.txt) < input.txt" ;\
	valgrind --tool=memcheck ../../$(EXECUTABLE) -Q $(shell cat $(1)/args.txt) < input.txt 2>&1 > /dev/null $(FILTER_MEMCHECK) || true ;\
	cd ../..

# Memcheck with ASCII test cases
.PHONY: memcheck
memcheck: debug $(TESTS_A_MEMCHECK)

.PHONY: $(TEST_A_OUT_DIR)/%_memcheck
$(TEST_A_OUT_DIR)/%_memcheck: $(TEST_A_DIR)/%
	$(call memcheck_dir,$<)

# Memcheck with Unicode test cases
.PHONY: memcheck_u
memcheck_u: debug $(TESTS_U_MEMCHECK)

.PHONY: $(TEST_U_OUT_DIR)/%_memcheck
$(TEST_U_OUT_DIR)/%_memcheck: $(TEST_U_DIR)/%
	$(call memcheck_dir,$<)


# Misc ==========================================

# Remove output directories
.PHONY: clean
clean:
	rm -rf $(BIN_DIR) $(BUILD_DIR) $(DOC_DIR)

# Install a copy of the executable in ~/bin
install: release | $(INST_DIR)
	cp -p $(EXECUTABLE) $(INST_DIR)

# Uninstall the copy of the executable in ~/bin
.PHONY: uninstall
uninstall:
	rm -f $(INST_DIR)/$(APP_NAME)


# Directories ===================================
# Rules for creating output directories

$(BIN_DIR):
	mkdir -p $@

$(BUILD_DIR):
	mkdir -p $@

$(TEST_A_OUT_DIR):
	mkdir -p $@

$(TEST_U_OUT_DIR):
	mkdir -p $@

$(INST_DIR):
	mkdir -p $@


# 3rd party =====================================

icdiff: $(HOME)/bin/icdiff

$(HOME)/bin/icdiff:
	mkdir -p $(HOME)/bin
	wget -q https://raw.githubusercontent.com/jeffkaufman/icdiff/release-1.9.2/icdiff -O $@
	chmod +x $@
