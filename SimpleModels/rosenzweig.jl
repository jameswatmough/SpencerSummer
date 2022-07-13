using BifurcationKit, Plots, Setfield, Parameters, LinearAlgebra, ForwardDiff, Revise
const BK = BifurcationKit

# sup norm
norminf(x) = norm(x, Inf)

function rosenzweig!(dx,x,p,t)
	@unpack μ,γ,κ = p
	dx[1] = x[1]*(1 - x[1]/κ) - x[1]*x[2]/(1.0+x[1])
	dx[2] = μ*x[2]*(γ*x[1]/(1.0+x[1]) - 1.0)
	return(dx)
end

rosenzweig(z, p) = rosenzweig!(similar(z), z, p, 0)

# sample parameter values 
par_base = (μ=1., γ=2.0, κ=0.5)

x0 = [par_base.κ,0.]
tspan = (0.,100.)

jet = BK.getJet(rosenzweig; matrixfree=false)

# continuation options
opts_br = ContinuationPar(pMin = 0.0, pMax = 5.,
	# parameters to have a smooth result
	ds = 0.04, dsmax = 0.05,
	# this is to detect bifurcation points precisely with bisection
	detectBifurcation = 3,
	# Optional: bisection options for locating bifurcations
	nInversion = 8, maxBisectionSteps = 25, nev = 3)

# continuation of equilibria
br, = continuation(jet[1], jet[2], x0, par_base, (@lens _.κ), opts_br;
	recordFromSolution = (x, p) -> (prey = x[1], predator = x[2]),
	tangentAlgo = BorderedPred(),
	plot = true, normC = norminf)

diagram = bifurcationdiagram(jet...,
	# initial point and parameter
	x0, par_base,
	# specify the continuation parameter
	(@lens _.κ),
	# very important parameter. This specifies the maximum amount of recursion
	# when computing the bifurcation diagram. It means we allow computing branches of branches
	# at most in the present case.
	3,
	(args...) -> setproperties(opts_br; pMin = 0.0, pMax = 5., ds = 0.04, dsmax = 0.05, nInversion = 8, detectBifurcation = 3, dsminBisection =1e-18, maxBisectionSteps=20);
	recordFromSolution = (x, p) -> (predator = x[2], prey = x[1]))

plot(diagram; putspecialptlegend=false, markersize=2, plotfold=false, title = "#branches = $(size(diagram))")
