# Makefile for multipole simulation

# Compiler
FC = gfortran

# Compiler flags
FFLAGS = -O3 -march=native -ffast-math -funroll-loops -flto
FFLAGS += -Wall -Wextra -Wpedantic
# FFLAGS += -g -fcheck=all -fbacktrace  # Uncomment for debugging

# Source files
MODULE = multipole_module.f90
FMM3D_MODULE = fmm3d_module.f90
FMM_OPS = fmm_operators.f90
MAIN_SRC = main.f90
TEST_ACCURACY_SRC = test_accuracy.f90
TEST_PERFORMANCE_SRC = test_performance.f90
TEST_FMM_PRECISION_SRC = test_fmm_precision.f90

OBJECTS = multipole_module.o main.o
FMM_OBJECTS = fmm3d_module.o fmm_operators.o
TARGET = multipole_sim
TEST_ACCURACY = test_accuracy
TEST_PERFORMANCE = test_performance
TEST_FMM_PRECISION = test_fmm_precision

# Python script
PYTHON = python3
COEFF_SCRIPT = compute_coefficients.py

# Default target
all: $(TARGET) fmm

# Compile the original program
$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $^
	@echo ""
	@echo "Build complete: $(TARGET)"
	@echo "Run with: ./$(TARGET)"

# Compile FMM precision test
fmm: $(TEST_FMM_PRECISION)
	@echo ""
	@echo "FMM high-precision test built: $(TEST_FMM_PRECISION)"
	@echo "Run with: ./$(TEST_FMM_PRECISION)"

# Test programs
$(TEST_ACCURACY): multipole_module.o $(TEST_ACCURACY_SRC)
	$(FC) $(FFLAGS) -o $@ multipole_module.o $(TEST_ACCURACY_SRC)

$(TEST_PERFORMANCE): multipole_module.o $(TEST_PERFORMANCE_SRC)
	$(FC) $(FFLAGS) -o $@ multipole_module.o $(TEST_PERFORMANCE_SRC)

$(TEST_FMM_PRECISION): $(FMM_OBJECTS) $(TEST_FMM_PRECISION_SRC)
	$(FC) $(FFLAGS) -o $@ fmm3d_module.o fmm_operators.o $(TEST_FMM_PRECISION_SRC)

# Module dependencies
main.o: multipole_module.o
test_accuracy.o: multipole_module.o
test_performance.o: multipole_module.o
fmm_operators.o: fmm3d_module.o
test_fmm_precision.o: fmm3d_module.o fmm_operators.o

# Compile rules
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# Run tests
test: $(TEST_ACCURACY) $(TEST_PERFORMANCE)
	@echo "Running accuracy test..."
	./$(TEST_ACCURACY)
	@echo ""
	@echo "Running performance benchmark..."
	./$(TEST_PERFORMANCE)

# Run FMM precision test
test-fmm: $(TEST_FMM_PRECISION)
	@echo "Running FMM high-precision test..."
	./$(TEST_FMM_PRECISION)

# Generate coefficients (optional, for future use)
coefficients:
	@mkdir -p coefficients
	$(PYTHON) $(COEFF_SCRIPT) 20
	@echo "Coefficients generated for p_max=20"

# Run the simulation
run: $(TARGET)
	./$(TARGET)

# Visualize results
viz: run
	$(PYTHON) visualize.py

# Clean build files
clean:
	rm -f *.o *.mod $(TARGET) $(TEST_ACCURACY) $(TEST_PERFORMANCE) $(TEST_FMM_PRECISION)

# Clean everything including output
cleanall: clean
	rm -rf output coefficients

# Create output directory
output:
	mkdir -p output

# Help
help:
	@echo "Multipole Simulation Makefile"
	@echo ""
	@echo "Targets:"
	@echo "  make              - Build the main program and FMM test"
	@echo "  make test         - Build and run accuracy & performance tests"
	@echo "  make test-fmm     - Build and run FMM high-precision test (10^-12)"
	@echo "  make run          - Build and run simulation"
	@echo "  make viz          - Run simulation and visualize"
	@echo "  make coefficients - Generate coefficient tables"
	@echo "  make clean        - Remove build files"
	@echo "  make cleanall     - Remove build files and output"
	@echo "  make help         - Show this help"

.PHONY: all run viz test test-fmm fmm clean cleanall coefficients help output
