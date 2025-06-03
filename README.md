# Cardiac Calcium Dynamics Simulation

A comprehensive 2D cardiac electrophysiology simulation modeling calcium handling, stochastic calcium sparks, and electrical propagation in cardiac tissue.

## Overview

This simulation models the complex interplay between electrical activity and calcium dynamics in cardiac myocytes using a spatially-resolved approach. The code implements a detailed cellular model with:

- **Electrical Activity**: Multiple ion channels (Na⁺, K⁺, Ca²⁺) and membrane voltage dynamics
- **Calcium Handling**: Two-compartment calcium model with boundary and interior regions  
- **Stochastic Processes**: Calcium spark generation and termination using binomial processes
- **Spatial Coupling**: 2D tissue with gap junction-mediated electrical propagation

## Code Architecture

The simulation consists of three main components that work together:

### 1. Main Simulation Engine (`ead1.f`)
- **Purpose**: Core simulation loop and spatial integration
- **Key Features**:
  - 2D grid of cardiac cells (5×5 default, configurable)
  - Time-stepping with adaptive integration for voltage dynamics
  - Spatial diffusion of membrane voltage between cells
  - Stimulus protocol management
  - Data output and analysis

### 2. Voltage Current Modules (`vcurrents.f`)
- **Purpose**: Implements detailed ionic current calculations for membrane voltage
- **Channels Included**:
  - `ina()`: Fast sodium current (action potential upstroke)
  - `ikr()`: Rapid delayed rectifier K⁺ current
  - `iks()`: Slow delayed rectifier K⁺ current  
  - `ik1()`: Inward rectifier K⁺ current
  - `ito()`: Transient outward K⁺ current (fast + slow)
  - `inak()`: Na⁺/K⁺ pump current

### 3. Calcium Current Modules (`cacurrents.f`)
- **Purpose**: Handles calcium transport and stochastic spark processes
- **Key Functions**:
  - `binom()`: Binomial random number generation
  - `binevol()`: Stochastic evolution of calcium spark clusters
  - `uptake()`: SERCA pump calcium uptake with Hill kinetics
  - `total()`/`xfree()`: Calcium buffering calculations
  - `inaca()`: Na⁺/Ca²⁺ exchanger current
  - `ica()`: L-type calcium current with Goldman-Hodgkin-Katz formulation
  - `markov()`: 10-state Markov model for L-type calcium channels

## How the Components Interact

The main simulation engine (`ead1.f`) orchestrates the entire simulation by:

1. **Calling Voltage Currents**: At each time step, `ead1.f` calls functions from `vcurrents.f` to calculate all membrane currents (Na⁺, K⁺, pump currents)

2. **Calling Calcium Currents**: Simultaneously calls functions from `cacurrents.f` to handle:
   - L-type calcium current calculation (`ica()`)
   - Na⁺/Ca²⁺ exchanger current (`inaca()`) 
   - Calcium buffering and SR uptake (`total()`, `xfree()`, `uptake()`)
   - Stochastic calcium spark dynamics (`binom()`, `binevol()`)
   - L-type channel Markov state evolution (`markov()`)

3. **Integrating Dynamics**: Combines all currents to update membrane voltage and calcium concentrations

4. **Spatial Coupling**: Applies diffusion operators to propagate electrical activity across the 2D tissue

### 1. Main Simulation Engine (`ead1.f`)
- **Purpose**: Core simulation loop and spatial integration
- **Key Features**:
  - 2D grid of cardiac cells (5×5 default, configurable)
  - Time-stepping with adaptive integration for voltage dynamics
  - Spatial diffusion of membrane voltage between cells
  - Stimulus protocol management
  - Data output and analysis

### 2. Voltage Current Modules (`vcurrents.f`)
- **Purpose**: Implements detailed ionic current calculations for membrane voltage
- **Channels Included**:
  - `ina()`: Fast sodium current (action potential upstroke)
  - `ikr()`: Rapid delayed rectifier K⁺ current
  - `iks()`: Slow delayed rectifier K⁺ current  
  - `ik1()`: Inward rectifier K⁺ current
  - `ito()`: Transient outward K⁺ current (fast + slow)
  - `inak()`: Na⁺/K⁺ pump current

### 3. Calcium Current Modules (`cacurrents.f`)
- **Purpose**: Handles calcium transport and stochastic spark processes
- **Key Functions**:
  - `binom()`: Binomial random number generation
  - `binevol()`: Stochastic evolution of calcium spark clusters
  - `uptake()`: SERCA pump calcium uptake with Hill kinetics
  - `total()`/`xfree()`: Calcium buffering calculations
  - `inaca()`: Na⁺/Ca²⁺ exchanger current
  - `ica()`: L-type calcium current with Goldman-Hodgkin-Katz formulation
  - `markov()`: 10-state Markov model for L-type calcium channels

## Physical Model

### Cellular Structure
The model represents cardiac myocytes with sparse t-tubules where:
- **Boundary Region**: Cell periphery with calcium influx from L-type channels
- **Interior Region**: Cell center with calcium diffusion from boundary
- **SR Compartments**: Separate sarcoplasmic reticulum stores for each region

### Calcium Spark Dynamics
- **Spark Clusters**: 3000 total boundary clusters, subset can be sparking
- **Stochastic Evolution**: Binomial processes for spark initiation/termination
- **L-type Channel Coupling**: Calcium channel opening triggers spark recruitment
- **RyR Release**: Calcium-induced calcium release from sarcoplasmic reticulum

### Electrical Coupling
- **Gap Junctions**: Ohmic coupling between neighboring cells
- **Diffusion**: Explicit finite difference scheme with no-flux boundaries
- **Stability**: Dual-step diffusion for numerical stability

## Compilation and Running

### Prerequisites
- Intel Fortran Compiler (`ifx`) or compatible Fortran compiler
- UNIX/Linux environment (for file I/O)

### Compilation
```bash
# Compile all modules together
ifx -o cardiac_sim ead1.f vcurrents.f cacurrents.f

# Alternative with optimization
ifx -O2 -o cardiac_sim ead1.f vcurrents.f cacurrents.f
```

### Execution
```bash
./cardiac_sim
```

### Output Files
- `v1x.dat`: Voltage time series from center cell (normal case)
- `v2x.dat`: Voltage time series (enhanced RyR2 leak case)

## Model Parameters

### Key Physiological Parameters
- **Temperature**: 35°C (308 K)
- **Ion Concentrations**: 
  - [Na⁺]ₒ = 136 mM, [Na⁺]ᵢ = ~14 mM (cycle length dependent)
  - [K⁺]ₒ = 5.4 mM, [K⁺]ᵢ = 140 mM
  - [Ca²⁺]ₒ = 1.8 mM
- **Pacing**: 500 ms cycle length (default)
- **Grid Size**: 5×5 cells (75 μm × 75 μm with 15 μm spacing)

### Model Variants
Set `mk` parameter in main code:
- `mk = 1`: Normal calcium handling
- `mk = 2`: Enhanced RyR2 leak (increased spark sensitivity)

## Scientific Applications

### Research Areas
- **Arrhythmogenesis**: How calcium sparks trigger arrhythmias
- **Excitation-Contraction Coupling**: Calcium dynamics during normal beats
- **Drug Effects**: Modifying channel conductances or calcium handling
- **Disease States**: Altered RyR2 function, calcium overload conditions

### Measurable Outputs
- Action potential morphology and duration
- Calcium transient amplitude and kinetics  
- Spark frequency and spatial distribution
- Propagation velocity and conduction patterns

## Code Structure Details

### Main Integration Loop
1. **Stimulus Application**: Brief current pulse at cycle start
2. **Ion Channel Updates**: Calculate all membrane currents
3. **Voltage Integration**: Adaptive time-stepping for fast dynamics
4. **Calcium Dynamics**: Update concentrations and spark states
5. **Spatial Diffusion**: Couple cells via gap junctions
6. **Data Collection**: Record voltage and calcium variables

### Numerical Methods
- **Time Integration**: Explicit Euler with adaptive sub-stepping
- **Spatial Discretization**: Second-order finite differences
- **Stochastic Processes**: Exact binomial sampling
- **Markov Chains**: Exponential integration for channel states

## Validation and Testing

### Expected Behavior
- **Normal Beat**: ~300 ms action potential duration, calcium transient peak ~1-2 μM
- **Spatial Propagation**: Uniform activation across 5×5 grid
- **Spark Activity**: Stochastic calcium sparks during calcium transient
- **Stability**: Steady-state reached after several beats

### Common Issues
- **Segmentation Faults**: Check array bounds and initialization
- **Long Line Errors**: Fortran compilers may have line length limits
- **Numerical Instability**: Reduce time step if voltages become unphysical

## References

This model builds on established frameworks for cardiac electrophysiology:
- Ion channel formulations from Grandi et al. and Shannon et al.
- Calcium handling based on Sobie et al. and Williams et al.
- Stochastic calcium sparks following Stern et al. and Rios et al.

## License

[Specify license here - typically MIT, GPL, or academic use]

## Contributing

[Guidelines for contributions, bug reports, and feature requests]

## Contact

[Contact information for questions and collaboration]
