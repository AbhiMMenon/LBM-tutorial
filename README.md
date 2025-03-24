# Lattice Boltzmann D2Q9: flow past cylinder

Fortran implementation of Philip Mocz's python tutorial. Single processor for now.
For grid size 400x100 tau=0.6 (relaxation time to equilibrium density
function),  tau=0.7 for 800x200 grid to preserve viscous behaviour (also
more stable)

## Requirements

- gfortran
- make
- gnuplot
- ffmpeg

# Compile and run 
  
    mkdir data obj
    make
    ./LBMSolver

Writes ~600Mb raw binary data to `data`. Visualize vorticity with 

    make movie

![](output.gif)




