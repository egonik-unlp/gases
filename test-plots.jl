include("sampling.jl")
using .MaxwellBoltzmann, Unitful, Plots, UnicodePlots
const m = 6.6335209e-20
gr()

temperatures = [120, 298, 320, 350]

p = Plots.plot()
p2 = Plots.plot()
for temperature ∈ temperatures
    rvs(temperature, 20000) |> ustrip |> v -> histogram!(p, v, label="$temperature K", alpha = .3)
    rvs(temperature, 20000) |> ustrip |> v -> 1//2*m*v.^2 |> Eₖ -> histogram!(p2, Eₖ, label="$temperature K", alpha = .3)
end

p
p2