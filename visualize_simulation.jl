ENV["GKSwstype"] = "100"

using Pkg; Pkg.instantiate()

using JLD2
using FourierFlows
using Plots
using Printf
using FFTW: irfft

filename = "./data/euler_longtime.jld2"
filename_diags = "./data/euler_longtime_diags.jld2"

normalize_diags_with_initial = false

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

moviegif_filename = "./movies/euler_longtime.gif"
moviegif_filename = uniquegifpath(moviegif_filename)
moviemp4_filename = moviegif_filename[1:end-4]*".mp4"

# ## Visualizing the simulation

function plot_output(x, y, ζ, ψ, t, k₀, ν, t_final)
    tν = ν * k₀^2 * t
    
    p_ζ = heatmap(x, y, ζ',
             aspectratio = 1,
                       c = :balance,
                    #clim = (-2, 2),
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

    p_diags1 = plot(4, # this means "a plot with two series"
                   label = ["|u|₂ / |u|₂(t=0)"  "|ζ|₂ / |ζ|₂(t=0)" "|ζ|₄ / |ζ|₄(t=0)" "|∇ζ|₂ / |∇ζ|₂(t=0)"],
                  legend = :topright,
               linewidth = 2,
                   alpha = 0.7, 
                  xlabel = "ν k₀² t",
                   xlims = (0, 1.01 * ν * k₀^2 * t_final),
                   # ylims = (0, 3),
                   # yscale = :log10
                   )
                   
    p_diags2 = plot(2, # this means "a plot with two series"
                   label = ["∂²ψ/∂x∂y(0, 0)" "(∂²/∂x²-∂²/∂y²)ψ(0, 0)"],
                  legend = :topright,
               linewidth = 2,
                   alpha = 0.7, 
                  xlabel = "ν k₀² t",
                   xlims = (0, 1.01 * ν * k₀^2 * t_final),
                   # yscale = :log10,
                   # ylims = (0, 3),
                   )

    p_diags3 = plot(3, # this means "a plot with two series"
                   label = ["max(|∂²ψ/∂x²|)" "max(|∂²ψ/∂y²|)" "max(|∂²ψ/∂x∂y|)"],
                  legend = :topright,
               linewidth = 2,
                   alpha = 0.7, 
                  xlabel = "ν k₀² t",
                   xlims = (0, 1.01 * ν * k₀^2 * t_final),
                   # yscale = :log10,
                   # ylims = (0, 3),
                   )
                   
    l = @layout [ Plots.grid(1, 2)
                  c{0.22h}
                  d{0.22h}
                  e{0.22h} ]
             
    p = plot(p_ζ, p_ψ, p_diags1, p_diags2, p_diags3, layout = l, size = (900, 1000))

    return p
end

diags = jldopen(filename_diags)


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

energyh = @. 1 / 2 * grid.invKrsq * abs2(ζh_initial)
energy_initial = 1 / (grid.Lx * grid.Ly) * FourierFlows.parsevalsum(energyh, grid)
U = sqrt(2*energy_initial)

global Re = U * grid.Lx / ν

global p = plot_output(x, y, ζ, ψ, 0.0, k₀, ν, t_final)

global total_iterations = length(iterations)

startwalltime = time()

global ∂²ψ∂x² = zeros(grid.nx, grid.ny)
global ∂²ψ∂y² = zeros(grid.nx, grid.ny)
global ∂²ψ∂x∂y = zeros(grid.nx, grid.ny)

global ∂²ψ∂x²h = zeros(Complex{Float64}, grid.nkr, grid.nl)
global ∂²ψ∂y²h = zeros(Complex{Float64}, grid.nkr, grid.nl)
global ∂²ψ∂x∂yh = zeros(Complex{Float64}, grid.nkr, grid.nl)

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
    
  hypeh = @. (grid.kr^2 - grid.l^2) * ψh
  hype = irfft(deepcopy(hypeh), grid.nx)
  hype_00 = hype[Int(grid.nx/2), Int(grid.ny/2)]

  @. ∂²ψ∂x²h = -grid.kr^2 * ψh
  @. ∂²ψ∂y²h = -grid.l^2 * ψh
  @. ∂²ψ∂x∂yh = -grid.kr * grid.l * ψh
  
  local ∂²ψ∂x² = irfft(∂²ψ∂x²h, grid.nx)  
  local ∂²ψ∂y² = irfft(∂²ψ∂y²h, grid.nx)  
  local ∂²ψ∂x∂y = irfft(∂²ψ∂x∂yh, grid.nx)
    
  ∂²ψ∂x∂y_00 = ∂²ψ∂x∂y[Int(grid.nx/2), Int(grid.ny/2)]
  
  local ζ = irfft(deepcopy(ζh), grid.nx)
  local ψ = irfft(deepcopy(ψh), grid.nx)

  p[1][1][:z] = ζ
  p[1][:title] = "vorticity, tν = " * @sprintf("%.4f", tν)
  p[2][1][:z] = ψ
  p[2][:title] = "streamfunction"
  p[3][:title] = "Re = " * @sprintf("%.2f", Re)
  
  t_diags = diags["diags/energy/t"][(i-1)*nsubs+1:i*nsubs]
  t_diags = diags["diags/energy/t"][(i-1)*nsubs+1:i*nsubs]
  tν_diags = ν * k₀^2 * t_diags

  ΔE  = (diags["diags/energy/data"][(i-1)*nsubs+1:i*nsubs]diags["diags/energy/data"][1]).^(1/2)
  ΔΖ₂ = (diags["diags/enstrophyL2/data"][(i-1)*nsubs+1:i*nsubs]).^(1/2)
  ΔΖ₄ = (diags["diags/enstrophyL4/data"][(i-1)*nsubs+1:i*nsubs]).^(1/4)
  ΔP  = (diags["diags/palinstrophy/data"][(i-1)*nsubs+1:i*nsubs]).^(1/2)
  
  if normalize_diags_with_initial == true
    ΔE  /= (diags["diags/palinstrophy/data"][1])^(1/2)
    ΔΖ₂ /= (diags["diags/enstrophyL2/data"][1])^(1/2)
    ΔΖ₄ /= (diags["diags/enstrophyL4/data"][1])^(1/4)
    ΔP  /= (diags["diags/palinstrophy/data"][1])^(1/2)
  end
  
  push!(p[3][1], tν_diags, log10.(ΔE))
  push!(p[3][2], tν_diags, log10.(ΔΖ₂))
  push!(p[3][3], tν_diags, log10.(ΔΖ₄))
  push!(p[3][4], tν_diags, log10.(ΔP))
  push!(p[4][1], tν, ∂²ψ∂x∂y_00)
  push!(p[4][2], tν, hype_00)
  push!(p[5][1], tν, maximum(abs.(∂²ψ∂x²)))
  push!(p[5][2], tν, maximum(abs.(∂²ψ∂y²)))
  push!(p[5][3], tν, maximum(abs.(∂²ψ∂x∂y)))
end

gif(anim, moviegif_filename, fps=14)
mp4(anim, moviemp4_filename, fps=14)
