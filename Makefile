CXX := g++

BUILD ?= RELEASE # set to DEBUG for debug build

ifeq ($(BUILD),DEBUG)
    CXXFLAGS := -std=c++20 -Wall -Wextra -pedantic -g -O0 -Iinc
else
    CXXFLAGS := -std=c++20 -Wall -Wextra -pedantic -O3 -march=native -DNDEBUG -Iinc
endif

LDFLAGS :=

# Directories
SRC_DIR := src
INC_DIR := inc
TEST_DIR := test
BUILD_DIR := build
BIN_DIR := bin

# Source files
SRC := $(wildcard $(SRC_DIR)/*.cpp)
TESTS := $(wildcard $(TEST_DIR)/*.cpp)

# Object files for library
OBJ := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRC))

# Test executables
TEST_BINS := $(patsubst $(TEST_DIR)/%.cpp,$(BIN_DIR)/%.out,$(TESTS))

.PHONY: all clean tests

all: tests

tests: $(TEST_BINS)

# Compile library source files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Build each test executable
$(BIN_DIR)/%.out: $(TEST_DIR)/%.cpp $(OBJ)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< $(OBJ) -o $@ $(LDFLAGS)

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)