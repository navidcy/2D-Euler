# # 2D decaying turbulence
#
# A simulation of decaying two-dimensional turbulence.

using Pkg; Pkg.instantiate()

using FourierFlows, Printf, Random, JLD2
 
using Random: seed!
using FFTW: rfft, irfft
using LinearAlgebra: ldiv!
using FourierFlows: parsevalsum

import GeophysicalFlows.TwoDNavierStokes
import GeophysicalFlows.TwoDNavierStokes: energy, enstrophy
import GeophysicalFlows: peakedisotropicspectrum


# ## Choosing a device: CPU or GPU

dev = CPU()     # Device (CPU/GPU)

# ## Numerical, domain, and simulation parameters
#
# First, we pick some numerical and physical parameters for our model.

α = 1                      # aspect ratio parameter α = Ly / Lx
ny, Ly  = 2*2*256, 2π            # grid resolution and domain length
nx, Lx  = Int(2*2*256 / α), Ly / α


## Then we pick the time-stepper parameters
	 dt =0.5*1e-3  # timestep
nsubs = 1000    # number of steps between each plot
 
 ν = 0.2*1e-4
 
initial_data_wavenumber = 2.0
initial_data_bandwidth = 0.5

k₀ = initial_data_wavenumber

 tfinal = 600.0 #/ (ν * k₀^2)
 nsteps = Int(round(tfinal / dt))

 grid = TwoDGrid(dev, nx, Lx, ny, Ly)
 x, y = gridpoints(grid)
 
 # No Forcing
 calcF!(args...) = nothing

# Initial condition

# seed!(1234)

E0 = 0.1


function random_initial_condition(grid, initial_data_wavenumber, initial_data_bandwidth; energy=0.1)
	K = @. sqrt(grid.Krsq)
	forcing_spectrum = @. exp(-(K - initial_data_wavenumber)^2 / (2*initial_data_bandwidth^2))
	
	energy_before_normalization = parsevalsum(forcing_spectrum .* grid.invKrsq / 2, grid) / (grid.Lx * grid.Ly)	
	
	forcing_spectrum *= energy / energy_before_normalization
	
	ξ = exp.(2π * im * rand(eltype(grid), size(K)))
	
	qh = @. sqrt(forcing_spectrum) * ξ
	q  = irfft(qh, grid.nx)
	
	return q
end


ζ₀ = random_initial_condition(grid, initial_data_wavenumber, initial_data_bandwidth; energy=E0)

# confirm that energy is E0
# prob = TwoDNavierStokes.Problem(dev; nx=nx, Lx=Lx, ny=ny, Ly=Ly, ν=ν, dt=dt, calcF=calcF!, stepper="ETDRK4")
# TwoDNavierStokes.set_ζ!(prob, ζ₀)
# dx, dy = Lx/nx, Ly/ny
# println(sum(@. (prob.vars.u^2 + prob.vars.v^2)/2 * dx * dy / (Lx * Ly)))
 
# plot initial condition
# using Plots
# p =heatmap(grid.x, grid.y, ζ₀', aspectratio=1)
# display(p)

include("simulation.jl")