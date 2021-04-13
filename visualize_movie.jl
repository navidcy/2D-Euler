ENV["GKSwstype"] = "100"

using Pkg; Pkg.instantiate()

using JLD2
using FourierFlows
using Plots
using Printf
using FFTW: irfft

filename = "./data/euler_longtime_3.jld2"

withoutgif(path) = (length(path)>3 && path[end-3:end] == ".gif") ? path[1:end-4] : path

"""
    uniquepath(path)

Returns `path` with a number appended if `isfile(path)`, incremented until `path` does not exist.
"""
function uniquegifpath(path)
  n = 1
  if isfile(path)
    path = withoutgif(path) * "_$n.gif"
  end
  while isfile(path)
    n += 1
    path = withoutgif(path)[1:end-length("_$(n-1)")] * "_$n.gif"
  end
  
  return path
end

moviegif_filename = "./movies/2deuler.gif"
moviegif_filename = uniquegifpath(moviegif_filename)
moviemp4_filename = moviegif_filename[1:end-4]*".mp4"

# ## Visualizing the simulation

function plot_output(x, y, ζ, ψ, t, k₀, ν, t_final)
    tν = ν * k₀^2 * t
    
    p_ζ = heatmap(x, y, ζ',
             aspectratio = 1,
                       c = :balance,
                    clim = (-2, 2),
                   xlims = (x[1], x[end]),
                   ylims = (y[1], y[end]),
                  xlabel = "x",
                  ylabel = "y",
                   title = "vorticity, tν = "*@sprintf("%.2f", tν),
              framestyle = :box)

    p_ψ = contourf(x, y, ψ',
               aspectratio = 1,
                    levels = 11,
                     xlims = (x[1], x[end]),
                     ylims = (y[1], y[end]),
                    xlabel = "x",
                    ylabel = "y",
                     title = "streamfunction",
                framestyle = :box)
                   
    l = @layout Plots.grid(1, 2)
             
    p = plot(p_ζ, p_ψ, layout = l, size = (900, 500))

    return p
end


# now let's make a movie of the flow fields

file = jldopen(filename)

nx, ny = file["grid/nx"], file["grid/ny"]
Lx, Ly = file["grid/Lx"], file["grid/Ly"]

global ν = file["params/ν"]
global k₀ = file["params/k₀"]
global nsubs = file["params/nsubs"]

grid = TwoDGrid(nx, Lx, ny, Ly)

x, y = grid.x, grid.y

iterations = parse.(Int, keys(file["snapshots/t"]))
final_iteration = iterations[end]

ζh_initial = file["snapshots/zetah/0"]
t_final = file["snapshots/t/$final_iteration"]

ζh = zeros(Complex{Float64}, (grid.nkr, grid.nl))
ψh = zeros(Complex{Float64}, (grid.nkr, grid.nl))

@. ζh = ζh_initial
ζ = irfft(deepcopy(ζh), grid.nx)
@. ψh = @. -grid.invKrsq * ζh
ψ = irfft(deepcopy(ψh), grid.nx)

global p = plot_output(x, y, ζ, ψ, 0.0, k₀, ν, t_final)

global total_iterations = length(iterations)

startwalltime = time()

anim = @animate for (i, iteration) in enumerate(iterations[1:end-1])
  
  if i%25== 0
    estimated_remaining_walltime = (time()-startwalltime)/60 / i * (total_iterations - i)
    log = @sprintf("estimated remaining walltime: %.2f min", estimated_remaining_walltime)
    println(log)
  end
  
  local t = file["snapshots/t/$iteration"]
  local tν = ν * k₀^2 * t
  
  local ζh = file["snapshots/zetah/$iteration"] 
  
  @. ψh = -grid.invKrsq * ζh
  
  local ζ = irfft(ζh, grid.nx)
  local ψ = irfft(ψh, grid.nx)

  p[1][1][:z] = ζ
  p[1][:title] = "vorticity, tν = " * @sprintf("%.4f", tν)
  p[2][1][:z] = ψ
  p[2][:title] = "streamfunction"
end

gif(anim, moviegif_filename, fps=14)
mp4(anim, moviemp4_filename, fps=14)
