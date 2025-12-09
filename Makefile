# Makefile for multipole simulation

# Compiler
FC = gfortran

# Compiler flags
FFLAGS = -O3 -march=native -ffast-math -funroll-loops -flto
FFLAGS += -Wall -Wextra -Wpedantic
# FFLAGS += -g -fcheck=all -fbacktrace  # Uncomment for debugging

# Source files
SOURCES = multipole_module.f90 main.f90
OBJECTS = $(SOURCES:.f90=.o)
TARGET = multipole_sim

# Python script
PYTHON = python3
COEFF_SCRIPT = compute_coefficients.py

# Default target
all: $(TARGET)

# Compile the program
$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $^
	@echo ""
	@echo "Build complete: $(TARGET)"
	@echo "Run with: ./$(TARGET)"

# Module dependencies
main.o: multipole_module.o

# Compile rules
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

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
	rm -f *.o *.mod $(TARGET)

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
	@echo "  make              - Build the program"
	@echo "  make run          - Build and run simulation"
	@echo "  make viz          - Run simulation and visualize"
	@echo "  make coefficients - Generate coefficient tables"
	@echo "  make clean        - Remove build files"
	@echo "  make cleanall     - Remove build files and output"
	@echo "  make help         - Show this help"

.PHONY: all run viz clean cleanall coefficients help output
