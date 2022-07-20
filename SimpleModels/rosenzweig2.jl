using BifurcationKit, Plots, Setfield, Parameters, LinearAlgebra, ForwardDiff, Revise
const BK = BifurcationKit

# sup norm
norminf(x) = norm(x, Inf)

recordFromSolution(x, p) = (x1 = x[1], x2 = x[2])

function rosenzweig!(dz,z,p,t)
    x1,x2 = z
	@unpack μ,γ,κ = p
	dz[1] = x1*(1 - x1/κ) - x1*x2/(1.0+x1)
	dz[2] = μ*x2*(γ*x1/(1.0+x1) - 1.0)
	dz
end

rosenzweig(z, p) = rosenzweig!(similar(z), z, p, 0)

# sample parameter values 
par_base = (μ=1., γ=2.0, κ=0.5)

x0 = [par_base.κ,1.]
tspan = (0.,100.)

jet = BK.getJet(rosenzweig; matrixfree=false)

# continuation options
opts_br = ContinuationPar(pMin = 0., pMax = 8.,
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
	normC = norminf)

diagram = bifurcationdiagram(jet...,
	# initial point and parameter
	x0, par_base,
	# specify the continuation parameter
	(@lens _.κ),
	# very important parameter. This specifies the maximum amount of recursion
	# when computing the bifurcation diagram. It means we allow computing branches of branches
	# at most in the present case.
	3,
	(args...) -> setproperties(opts_br; pMin = 0.0, pMax = 12., ds = 0.04, dsmax = 0.05, nInversion = 8, detectBifurcation = 3, dsminBisection =1e-18, maxBisectionSteps=20);
	recordFromSolution = (x, p) -> (predator = x[2], prey = x[1]),
    verbosity = 0, plot = false)

# branch of the diagram with Hopf point
brH = BK.getBranch(diagram, (1,)).γ

# newton parameters
optn_po = NewtonPar(tol = 1e-8,  maxIter = 25)

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.1, ds= -0.001, dsmin = 1e-4,
 newtonOptions = (@set optn_po.tol = 1e-8), precisionStability = 1e-2,
 detectBifurcation = 1, saveSolEveryStep=1)

Mt = 51 # number of time sections
	br_po, utrap = continuation(
	jet..., brH, 2, opts_po_cont,
	PeriodicOrbitTrapProblem(M = Mt);
	# help branching from Hopf
	usedeflation = true,
	# specific linear solver for ODEs
	linearPO = :Dense,
	recordFromSolution = (x, p) -> (xtt=reshape(x[1:end-1],2,Mt);
		return (max = maximum(xtt[1,:]),
			min = minimum(xtt[1,:]),
			period = x[end])),
	finaliseSolution = (z, tau, step, contResult; prob = nothing, kwargs...) -> begin
		# limit the period
		z.u[end] < 100
		end,
	normC = norminf)

plot(diagram, markersize=2, plotfold=false, title = "#branches = $(size(diagram))")



