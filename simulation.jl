# ## Problem setup
# We initialize a `Problem` by providing a set of keyword arguments. The
# `stepper` keyword defines the time-stepper to be used.
prob = TwoDNavierStokes.Problem(dev; nx=nx, Lx=Lx, ny=ny, Ly=Ly, ν=ν, dt=dt, calcF=calcF!, stepper="ETDRK4")

# Next we define some shortcuts for convenience.
sol, cl, vs, gr = prob.sol, prob.clock, prob.vars, prob.grid
x, y = gr.x, gr.y

TwoDNavierStokes.set_ζ!(prob, ζ₀)

energy_initial = energy(prob)

U = sqrt(2 * energy_initial)

Re = U * Ly / ν


# ## Diagnostics

"""
    vorticityL4(prob)
Returns the domain-averaged enstrophy, ∫ ζ⁴ dxdy / (Lx Ly), for the solution in `sol`.
"""
@inline function vorticityL4(prob)
  sol, vars, grid = prob.sol, prob.vars, prob.grid
  @. vars.ζh = sol
  ldiv!(vars.ζ, grid.rfftplan, vars.ζh)
  return sum(@. vars.ζ^4) * grid.dx * grid.dy / (grid.Lx * grid.Ly)
end

"""
    palinstrophy(prob)
Returns the domain-averaged palinstrophy, ∫ |∇ζ|² dxdy / (Lx Ly), for the solution in `sol`.
"""
@inline function palinstrophy(prob)
  sol, vars, grid = prob.sol, prob.vars, prob.grid
  palinstrophyh = vars.uh # use vars.uh as scratch variable

  @. palinstrophyh = grid.Krsq * abs2(sol)
  return 1 / (grid.Lx * grid.Ly) * parsevalsum(palinstrophyh, grid)
end


# Create Diagnostics -- `energy` and `enstrophy` functions are imported at the top.
E = Diagnostic(energy, prob; nsteps=nsteps)
Z2 = Diagnostic(enstrophy, prob; nsteps=nsteps)
Z4 = Diagnostic(vorticityL4, prob; nsteps=nsteps)
P = Diagnostic(palinstrophy, prob; nsteps=nsteps)
diags = [E, Z2, Z4, P] # A list of Diagnostics types passed to "stepforward!" will  be updated every timestep.


# ## Output

# We choose folder for outputing `.jld2` files and snapshots (`.png` files).
filepath = "./data/"
filename = joinpath(filepath, "euler_longtime.jld2")
filename_diags = joinpath(filepath, "euler_longtime_diags.jld2")

filename = FourierFlows.uniquepath(filename)
@info "Output will be saved at $filename."

filename_diags = FourierFlows.uniquepath(filename_diags)
@info "Diagnostics will be saved at $filename_diags."

# Do some basic file management
if isfile(filename); rm(filename); end
if isfile(filename_diags); rm(filename_diags); end

# And then create Output
get_sol(prob) = sol # extracts the Fourier-transformed solution
out = Output(prob, filename, (:zetah, get_sol))
saveproblem(out)


# ## Time-stepping the `Problem` forward

# We time-step the `Problem` forward in time.

saveoutput(out) # save initial condition

file = jldopen(filename, "a+")
FourierFlows.savefield(file, "params/nsubs", nsubs)
FourierFlows.savefield(file, "params/k₀", k₀)
close(file)

@info "Starting simulation..."

startwalltime = time()

for j=0:Int(round(nsteps/nsubs))-1
    
  if j%(1000/nsubs)==0
    TwoDNavierStokes.updatevars!(prob)
    
    cfl = cl.dt * maximum([maximum(vs.u) / gr.dx, maximum(vs.v) / gr.dy])
    
    estimated_remaining_walltime = (time()-startwalltime)/60 / cl.step * (nsteps-cl.step)
    log = @sprintf("step: %04d, t: %d, tν : %.4f, E: %.3e, Zₗ₂: %.3e, Zₗ₄: %.3e, P: %.3e, cfl: %.4f, walltime: %.2f min, estimated remaining walltime: %.2f min",
      cl.step, cl.t, ν*k₀^2*cl.t, E.data[E.i], Z2.data[Z2.i], Z4.data[Z4.i], P.data[P.i], cfl, (time()-startwalltime)/60, estimated_remaining_walltime)
    println(log)
  end  

  stepforward!(prob, diags, nsubs)
  dealias!(sol, gr)

  saveoutput(out)
end

@info @sprintf("Simulation finished after %.2f minutes.", (time()-startwalltime)/60)

savediagnostic(E, "energy", filename_diags)
savediagnostic(Z2, "enstrophyL2", filename_diags)
savediagnostic(Z4, "enstrophyL4", filename_diags)
savediagnostic(P, "palinstrophy", filename_diags)

@info "Run visualize_simulation.jl after prescribing the two filenames used for saving output at the top the visualize_simulation.jl script."

@info "Output for this simulation was saved at $filename and $filename_diags."
