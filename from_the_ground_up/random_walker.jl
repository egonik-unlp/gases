using Plots


const init_speed = rand(2)


struct Map
    bounds :: Array{Float64}
    steps :: Int64
    step :: Int64
    function Map(bounds, steps)
        new( bounds, steps, 0 )
    end
end



struct RandomWalker
    pos :: Array{Float64}
end

function main(steps, agents)
    world = Map([0 0; 1920 1080], 100)
    walker = RandomWalker(rand(2)* Map.bounds[2,:])
    anim = @animate for step ∈ steps
        step!(world, walker)
    end
end







function step!(world, walker)
    walker.pos += rand(-1:1, 2)
end




