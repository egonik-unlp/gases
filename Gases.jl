include("./sampling.jl")

using .MaxwellBoltzmann, Unitful, Distances, Plots


gr()


const T = 298u"K"
const dt = 1
	
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
	dist = pairwise(Euclidean(), s.pos', dims=2) |> d -> findall(d .< distancia) 
	colisiones = [index for index in dist if index[1] < index[2] ]# |> c -> getindex.(c, [1 2])
	for cl ∈ colisiones
		i,j = (cl[1], cl[2])
		pos_i, vel_i = s.pos[i], s.vel[i] 
		pos_j, vel_j = s.vel[j], s.vel[j]
		rel_pos, rel_vel = pos_i - pos_j, vel_i - vel_j
		r_rel = dot(rel_pos, rel_pos)
		v_rel = dot(rel_vel, rel_pos)
		v_rel = (2 .* rel_pos .*v_rel)./(r_rel .- rel_vel)
		v_cm = (vel_i .+ vel_j)./ 2
		s.vel[i] = v_cm - v_rel/2
		s.vel[j] = v_cm + v_rel/2
	end

end

function advance_time!(s::Simulation)
	s.pos += s.vel * dt
	colisiones!(s)
	s.step += dt
end

p = Plots.plot()

@gif for step in 1:s.steps
	scatter!(p, s.pos[:,1], s.pos[:,2])
	advance_time!(s)
end



Plots.savefig(p, "grafico.gif")