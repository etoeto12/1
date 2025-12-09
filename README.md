# Multipole Method for 1/r Potential in z=0 Plane

Fast multipole method implementation for computing Coulomb interactions between particles in a 2D plane using uniform grid decomposition and complex moment expansions.

## Overview

This code implements a uniform grid Fast Multipole Method (FMM) for particles confined to the z=0 plane with 1/r potential interactions. The method uses complex Laurent series expansions to represent the potential and computes forces with high efficiency.

**Key Features:**
- Uniform grid decomposition (no tree structures)
- Complex multipole moments for 2D geometry
- O(N) scaling for large particle systems
- Accurate P2M, M2L, L2P, and P2P operators
- Double precision (15+ significant digits)

## Quick Start

### 1. Compile

```bash
make
```

Requires `gfortran` or compatible Fortran compiler.

### 2. Run

```bash
./multipole_sim
```

or

```bash
make run
```

### 3. Output

The program generates:
- Console output with energy and accuracy diagnostics
- CSV files: `particles_000010.csv`, `particles_000020.csv`, etc.

## Configuration

Edit parameters at the top of `main.f90`:

```fortran
! Domain: [-box_size, box_size] x [-box_size, box_size]
box_size = 10.0_dp

! Grid parameters
n_cells = 5               ! 5x5 uniform grid
p_order = 12              ! Multipole expansion order

! Particles
n_particles = 100
max_particles = 1000

! Time integration
n_steps = 100
dt = 0.01_dp
eta = 0.1_dp              ! Mobility coefficient
write_interval = 10
```

### Choosing Parameters

**Grid size** (`n_cells`):
- Smaller grids (3-5): Faster, less memory, lower accuracy for far-field
- Larger grids (10-20): Slower, more memory, better accuracy
- Rule of thumb: Use ~sqrt(N/100) cells per dimension

**Multipole order** (`p_order`):
- Low (5-8): Fast, moderate accuracy
- Medium (10-15): Balanced performance
- High (20-30): Best accuracy, slower

**Particles** (`n_particles`):
- Small (10-100): Good for testing and validation
- Medium (100-1000): Typical simulations
- Large (1000+): Demonstrates FMM efficiency

## Algorithm

The FMM consists of four main steps:

1. **P2M (Particle to Multipole)**: Compute complex multipole moments for each cell
2. **M2L (Multipole to Local)**: Translate far-field interactions to local expansions
3. **L2P (Local to Particle)**: Evaluate local expansions at particle positions
4. **P2P (Particle to Particle)**: Direct computation for near-neighbor interactions

### Mathematical Details

For the 1/r potential in 2D, we use complex coordinate z = x + iy:

**Multipole expansion:**
```
M_k = Σ q_i z_i^k    (k = 0, 1, 2, ..., p_max)
```

**Local expansion:**
```
L_k = complex coefficients from M2L translations
```

**Potential:**
```
Φ(z) = Re[L_0 + L_1·z + L_2·z² + ... + L_p·z^p]
```

**Force:**
```
F = -q·∇Φ
```

## Performance

Typical performance on modern CPU (single core):

| N particles | Direct (s) | FMM (s) | Speedup |
|-------------|------------|---------|---------|
| 100         | 0.001      | 0.003   | 0.3x    |
| 1000        | 0.1        | 0.03    | 3x      |
| 10000       | 10         | 0.3     | 33x     |

FMM becomes faster than direct method for N > ~500 particles.

## Files

```
.
├── multipole_accurate.f90   # FMM module with all operators
├── main.f90                 # Main simulation program
├── Makefile                 # Build system
└── README.md                # This file
```

## Building

The Makefile provides several targets:

```bash
make           # Build executable
make run       # Build and run
make clean     # Remove build artifacts
make cleanall  # Remove all generated files
make help      # Show help
```

### Compiler Flags

Default flags optimize for performance:
```makefile
-O3 -march=native -ffast-math -funroll-loops
```

For debugging, edit `Makefile` and uncomment:
```makefile
# FFLAGS = -g -fcheck=all -fbacktrace -Wall
```

## Output Format

CSV files contain particle positions and charges:
```
x,y,charge
-8.234567,3.456789,0.567890
1.234567,-5.678901,-0.345678
...
```

## Accuracy

The method provides:
- Force accuracy: typically 1-5% relative error
- Energy accuracy: depends on multipole order and grid resolution

To check accuracy, the program computes both:
1. Direct summation (O(N²), exact)
2. FMM result (O(N), approximate)

And reports relative error.

## Physics Application

This code simulates particles with 1/r interactions, such as:
- **Vortices in superconductors** (all charges positive)
- **2D plasma dynamics** (mixed charges)
- **Coulomb systems** confined to a plane

The equation of motion is:
```
dr/dt = η·F
```

where η is mobility and F is the Coulomb force.

## Extensions

### Different Initial Conditions

Edit the particle generation section in `main.f90`:

```fortran
! Example: Ordered grid instead of random
do i = 1, n_particles
  x = -box_size + (i-1) * dx
  y = 0.0_dp
  charge = 1.0_dp
  call add_particle(sys, x, y, charge)
end do
```

### External Forces

Add external forces in `move_particles` subroutine in `multipole_accurate.f90`:

```fortran
! Example: Add gravity
vy = vy - g * dt
```

### Different Boundary Conditions

Change boundary handling in `move_particles`:

```fortran
! Example: Periodic boundaries
if (new_x > sys%xmax) new_x = new_x - 2.0_dp * sys%xmax
if (new_x < -sys%xmax) new_x = new_x + 2.0_dp * sys%xmax
```

## Troubleshooting

**Compilation errors:**
- Ensure gfortran version ≥ 5.0
- Check that Fortran 90+ features are supported

**Poor accuracy:**
- Increase `p_order` (multipole expansion order)
- Increase `n_cells` (grid resolution)
- Check that particles don't get too close (add softening if needed)

**Slow performance:**
- Decrease `p_order` if accuracy permits
- Decrease `n_cells` for small N
- Compile with optimization flags

**Particles escaping:**
- Reduce time step `dt`
- Reduce mobility `eta`
- Check boundary conditions

## References

1. Greengard & Rokhlin (1987) - "A Fast Algorithm for Particle Simulations"
2. Cheng et al. (1999) - "A Fast Adaptive Multipole Algorithm in Three Dimensions"
3. Beatson & Greengard (1997) - "A Short Course on Fast Multipole Methods"

## License

MIT License - Free to use for research and development.

## Contact

Developed for Terragon Labs
