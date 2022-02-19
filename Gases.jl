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
	#gases :: Union{Gas,Vector{Gas}}
	vel :: Array{Float64}
	pos :: Array{Float64}
	mass :: Array{Float64}

	function Simulation(steps ,gases...)
		vel = vcat([g.vel for g in gases]...) 
		pos = vcat([g.pos for g in gases]...) 
		mass = vcat([fill(g.mass, g.N) for g in gases]...)
		new(steps,1, vel, pos, mass)
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


function colisiones!(s:: Simulation, distancia = .1)
	#la distancia que elegi como prueba de colisiones es random
	data = [vcat(gas.pos', fill(gas.mass, size(gas.pos)[1])', gas.vel') for gas in s.gases] |> p -> hcat(p...)
	# En este array la estructura es :
	# 1 - pos_x
	# 2 - pos_y
	# 3 - masa
	# 4 - vel_x
	# 5 - vel_y
	dist = pairwise(Euclidean(), data[begin:2,:], dims=2) |> d -> findall(d .< distancia) 
	colisiones = [index for index in dist if index[1] < index[2] ]# |> c -> getindex.(c, [1 2])
	map(cls) do cl
           i,j = (cl[1], cl[2])
           pos_i, vel_i = data[1:2,i], data[4:5,i] 
           pos_j, vel_j = data[1:2,j], data[4:5,j]
           rel_pos, rel_vel = pos_i - pos_j, vel_i - vel_j
           r_rel = dot(rel_pos, rel_pos)
           v_rel = dot(rel_vel, rel_pos)
           v_rel = (2 .* rel_pos .*v_rel)./(r_rel .- rel_vel)
           v_cm = (vel_i .+ vel_j)./ 2
           v_rel, v_cm
	end
end


function colisiones_(s::Simulation)
	data = [vcat(gas.pos', fill(gas.mass, size(gas.pos)[1])', gas.vel') for gas in s.gases] |> p -> hcat(p...)

	dist = pairwise(Euclidean(), data[begin:2,:], dims=2) |> d -> findall(d .< .1) 
	[index for index in dist if index[1] < index[2] ]
end

function build_data(s::Simulation)
	[vcat(gas.pos', fill(gas.mass, size(gas.pos)[1])', gas.vel') for gas in s.gases] |> p -> hcat(p...)
end





function solve_colisiones(g::Gas)
	# Se tiene que tener en cuenta que podria haber colisiones entre distintas particulas que no 
	# tienen por que ser del mismo gas. Esto implicaria pensar que habria que concatenar de alguna 	      # las posiciones de todas las partículas.

	nothing
end




function advance_time!(s::Simulation)
	nothing	
end

