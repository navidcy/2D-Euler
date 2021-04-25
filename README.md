# 2D Euler

Study of the behavior of Euler equations at long times. The code is written in [Julia](https://julialang.org) and utilizes the Julia package [GeophysicalFlows.jl](http://github.com/FourierFlows/GeophysicalFlows.jl).

## Installation

First you need to [install Julia](https://julialang.org/downloads/). Then clone the repository, e.g.,

```
git clone https://github.com/navidcy/2D-Euler/.git
```

Enter the directory you've cloned the repository, e.g., 

```
cd 2D-Euler/
```

Then, after you edit `setup_and_run_simulation.jl` file with your parameter values, and while still inside the repository's main directory, you may run a simulation via

```
julia --project setup_and_run_simulation.jl
```

To make an animation of the simulation, first edit the `visualize_simulation.jl` or `visualize_movie.jl` script to point to the correct `.jld2` output files from the simulation you want to visualize, and then run, e.g.,

```
julia --project visualize_simulation.jl
```

## Cite

You can cite the [GeophysicalFlows.jl](http://github.com/FourierFlows/GeophysicalFlows.jl) package as:

> Constantinou et al., (2021). GeophysicalFlows.jl: Solvers for geophysical fluid dynamics problems in periodic domains on CPUs & GPUs. _Journal of Open Source Software_, **6(60)**, 3053, doi:[10.21105/joss.03053](https://doi.org/10.21105/joss.03053)
