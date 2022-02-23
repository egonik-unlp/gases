#  Implementanción en Julia del generador de números aleatorios  
# desarrollado en el trabajo:
#  "Efficient Algorithm for Generating Maxwell Random Variables" de Nader M.A. Mohamed
# 
# J Stat Phys (2011) 145:1653–1660
# DOI 10.1007/s10955-011-0364-y
#
module MaxwellBoltzmann

export rvs
using Unitful



const k = 1.37 
const C = 0.5
const g = (2)/(k*C*sqrt(π))


function rvs(T :: Number) 
	r1,r2 = rand(), rand()
	y = - 2*log(r1)
	( (g^2*y) ≤ (r2/r1)^2 ? rvs(T) : (ustrip(T)*y)u"m/s") 
end


function rvs(T::Number, N::Int) 
	map( x -> rvs(T), 1:N) 
end


end

