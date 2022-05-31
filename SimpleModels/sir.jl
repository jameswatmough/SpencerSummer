
using DifferentialEquations, Plots

function sirs!(dx,x,p,t)
	β,ω,γ = p
	S,I = x
	dx[1] = -β*S*I + ω*(1-S-I)
	dx[2] =  β*S*I - γ*I
end
