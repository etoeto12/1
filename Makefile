# Makefile for multipole simulation
# Uniform grid FMM with complex moments for 1/r potential in z=0 plane

# Compiler
FC = gfortran

# Compiler flags
FFLAGS = -O3 -march=native -ffast-math -funroll-loops
FFLAGS += -Wall -Wextra

# Debug flags (uncomment for debugging)
# FFLAGS = -g -fcheck=all -fbacktrace -Wall

# Source files
MODULE_SRC = multipole_accurate.f90
MAIN_SRC = main.f90

# Object files
MODULE_OBJ = multipole_accurate.o
MAIN_OBJ = main.o

# Executable
TARGET = multipole_sim

# Default target
all: $(TARGET)

# Build executable
$(TARGET): $(MODULE_OBJ) $(MAIN_OBJ)
	$(FC) $(FFLAGS) -o $@ $^
	@echo ""
	@echo "=========================================="
	@echo "Build complete!"
	@echo "Run with: ./$(TARGET)"
	@echo "=========================================="

# Compile module
$(MODULE_OBJ): $(MODULE_SRC)
	$(FC) $(FFLAGS) -c $<

# Compile main program (depends on module)
$(MAIN_OBJ): $(MAIN_SRC) $(MODULE_OBJ)
	$(FC) $(FFLAGS) -c $<

# Run simulation
run: $(TARGET)
	@echo "Running simulation..."
	@echo ""
	./$(TARGET)

# Clean build artifacts
clean:
	rm -f *.o *.mod $(TARGET)
	@echo "Cleaned build artifacts"

# Clean everything including output files
cleanall: clean
	rm -f particles_*.csv
	@echo "Cleaned all files"

# Help
help:
	@echo "Multipole Simulation - Uniform Grid FMM"
	@echo ""
	@echo "Targets:"
	@echo "  make        - Build the program"
	@echo "  make run    - Build and run simulation"
	@echo "  make clean  - Remove build artifacts (.o, .mod, executable)"
	@echo "  make cleanall - Remove all files including output CSV files"
	@echo "  make help   - Show this help"
	@echo ""
	@echo "Usage:"
	@echo "  1. Edit parameters in main.f90"
	@echo "  2. Run 'make' to compile"
	@echo "  3. Run './multipole_sim' to execute"

.PHONY: all run clean cleanall help
