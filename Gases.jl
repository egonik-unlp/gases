include("./sampling.jl")

using .MaxwellBoltzmann, Unitful, Distances, Plots, LinearAlgebra, Distributions
const X, Y = 1,2
const rscale = 1e-6 # Escala
const tscale = 1e9 #

gr()


const T = 298u"K"
const dt = 1/15
	
struct Gas 
	N :: Int64
	pos :: Matrix
	vel :: Matrix
	radius :: Float64
	mass :: Float64
	T :: Int64
	function Gas(;N, radius, mass, T , lo )
		θ = 2π.*rand(N)
		V₀ = (rvs(T,N) |> ustrip)
		vel = map(V₀,θ) do velocity,angle
			velocity.* [cos(angle) , sin(angle)]
		end |> v -> hcat(v...)'
		vel *= ( rscale * 200 )
		pos = (lo == 1 ? rand(Uniform(0,0.5), N) : rand(Uniform(0.5,1),  N) ) |> x -> hcat(rand(N), x)
		new(N,pos,vel,rscale*radius,mass, T)
	end
end





mutable struct Simulation
	steps :: Int
	step :: Int
	gases :: Union{Gas,Vector{Gas}}
	vel :: Array{Float64}
	pos :: Array{Float64}
	mass :: Array{Float64}
	radius :: Array{Float64}
	function Simulation(steps ,gases...)
		vel = vcat([g.vel for g in gases]...) 
		pos = vcat([g.pos for g in gases]...) 
		mass = vcat([fill(g.mass, g.N) for g in gases]...)
		radius = vcat([fill(g.radius, g.N) for g in gases]...)
		gases = [gas for gas in gases]
		new(steps, 1, gases, vel, pos, mass, radius)
	end
end


# v(v) = sqrt.(v[:,1].^2 + v[:,2].^2)
# Eₖ(s::Simulation) = s.mass.*mv(s.vel).^2 |> sum




argon = Gas(;N = 100,
	 radius = 2e-10,
	 mass = 10^-7,
	 T = 298,
	 lo = 1)

argon2 = Gas(;N = 100,
	 radius = 2e-10,
	 mass = 10^-9,
	 T = 298,
	 lo = nothing)


s = Simulation(1200,argon,argon2)


function colisiones!(s:: Simulation)
	distancia = .0001
	#la distancia que elegi como prueba de colisiones es random
	dist = pairwise(Euclidean(), s.pos', dims=2) |> d -> findall(d .< distancia) 
	colisiones = [index for index in dist if index[1] < index[2] ]

	for cl ∈ colisiones
		i,j = (cl[1], cl[2])
		pos_i, vel_i = s.pos[i], s.vel[i] 
		pos_j, vel_j = s.vel[j], s.vel[j]
		rel_pos, rel_vel = pos_i - pos_j, vel_i - vel_j
		r_rel = dot(rel_pos, rel_pos)
		v_rel = dot(rel_vel, rel_pos)
		v_rel = (2 .* rel_pos .* v_rel)./(r_rel .- rel_vel)
		v_cm = (vel_i .+ vel_j)./ 2
		s.vel[i] = v_cm - v_rel/2
		s.vel[j] = v_cm + v_rel/2	
	end

	hit_left_wall = s.pos[:, X] .< s.radius
	hit_right_wall = s.pos[:, X] .> 1 .- s.radius
	hit_bottom_wall = s.pos[:, Y] .< s.radius
	hit_top_wall = s.pos[:, Y] .> 1 .- s.radius
	s.vel[hit_left_wall .|| hit_right_wall, X] *= -1
	s.vel[hit_bottom_wall .|| hit_top_wall, Y] *= -1

end

function advance_time!(s::Simulation)
	s.pos += s.vel * dt
	colisiones!(s)
	v(v) = sqrt.(v[:,1].^2 + v[:,2].^2)
	Eₖ(s::Simulation) = s.mass.*mv(s.vel).^2 |> sum
	# @info "Energía cinética en el paso $(s.step) = $(Eₖ(s))"
	s.step += 1
end

function main(s::Simulation)


	anim = @animate for step in 1:s.steps
		global p = Plots.plot(legend=:false)
		#global h = Plots.plot(legend=:false)
		xlims!(p,0,1)
		ylims!(p,0,1)
		Nᵢ = 1
		Nₛ = 0 
		cg = cgrad(:acton)
		for (col,gas) ∈ zip(cg, s.gases)
			Nₛ += size(gas.pos)[1]
			scatter!(p, s.pos[ Nᵢ:Nₛ ,1], s.pos[ Nᵢ:Nₛ ,2], color = col, label = :false)
			advance_time!(s) ###Esto no va acá ni en p2
			# labels!(" $(s.step*dt) s")
			Nᵢ += Nₛ
		end
		#histogram!(h, s.vel, label = :false)
		#plot(p,h, size = [1920, 1080])
	end



	gif(anim, "anim.gif", fps= 15)

	Plots.savefig(p, "grafico")
end


main(s)
