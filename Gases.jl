include("./sampling.jl")

using .MaxwellBoltzmann
using Unitful
using Distances


const T = 298u"K"
	
mutable struct Gas 
	N :: Int64
	pos :: Matrix
	vel :: Matrix
	radius :: Float64
	mass :: Float64
	T :: Int64
	function Gas(;N, radius, mass, T )
		θ = 2π.*rand(N)
		V₀ = rvs(T,N) |> ustrip
		vel = map(V₀,θ) do velocity,angle
			velocity.* [cos(angle) , sin(angle)]
		end |> v -> hcat(v...)'
		pos = rand(N,2)
		new(N,pos,vel,radius,mass, T)
	end
end





mutable struct Simulation
	steps :: Int
	step :: Int
	gases :: Union{Gas,Vector{Gas}}
	function Simulation(steps ,gases...)
		new(steps,1,[gases...])
	end
end


argon = Gas(;N = 100,
	 radius = 1,
	 mass = 10^-7,
	 T = 298)

argon2 = Gas(;N = 100,
	 radius = 1,
	 mass = 10^-9,
	 T = 298)


s = Simulation(1000,argon,argon2)


function colisiones!(s:: Simulation)
	#la distancia que elegi es random
	#primero construimos una matriz que incluya las masas de las particulas
	data = [vcat(gas.pos', fill(gas.mass, size(gas.pos)[1])', gas.vel') for gas in s.gases] |> p -> hcat(p...)

	dist = pairwise(Euclidean(), data[begin:2,:], dims=2) |> d -> findall(d .< .1) 
	colisiones = [index for index in dist if index[1] < index[2] ] |> c -> getindex.(c, [1 2])
	map(eachrow(colisiones)) do (px,py)
		println("velocidad en colision p1 = $(data[4,px]), p2 =$(data[4,py])")
		println("masa en colision p1 = $(data[3,px]), p2 =$(data[3,py])")
	end
	#colisiones
end

function solve_colisiones(g::Gas)
	# Se tiene que tener en cuenta que podria haber colisiones entre distintas particulas que no 
	# tienen por que ser del mismo gas. Esto implicaria pensar que habria que concatenar de alguna 	      # las posiciones de todas las partículas.

	nothing
end




function advance_time!(s::Simulation)
	nothing	
end

