# 2D Euler

Study of the behavior of Euler equations at long times. The code is written in [Julia](https://julialang.org) and utilizes the Julia package [GeophysicalFlows.jl](http://github.com/FourierFlows/GeophysicalFlows.jl).

## Installation

First you need to [install Julia](https://julialang.org/downloads/). We suggest using Julia version 1.6 or later.

Then clone the repository, e.g.,

```
git clone https://github.com/navidcy/2D-Euler.git
```

Enter the directory you've cloned the repository into, e.g., 

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

The bibtex entry is:

```bibtex
@article{GeophysicalFlowsJOSS,
  doi = {10.21105/joss.03053},
  url = {https://doi.org/10.21105/joss.03053},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {60},
  pages = {3053},
  author = {Navid C. Constantinou and Gregory LeClaire Wagner and Lia Siegelman and Brodie C. Pearson and André Palóczy},
  title = {GeophysicalFlows.jl: Solvers for geophysical fluid dynamics problems in periodic domains on CPUs \& GPUs},
  journal = {Journal of Open Source Software}
}
```
