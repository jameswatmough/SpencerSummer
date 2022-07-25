using Revise, Parameters, Setfield, Plots, LinearAlgebra
using BifurcationKit
const BK = BifurcationKit

# define the sup norm
norminf(x) = norm(x, Inf)

# function to record information from a solution
recordFromSolution(x, p) = (u1 = x[1], u2 = x[2])
####################################################################################################
function pp2!(dz, z, p, t)
	@unpack p1, p2, p3, p4 = p
	u1, u2 = z
	dz[1] = p2 * u1 * (1 - u1) - u1 * u2 - p1 * (1 - exp(-p3 * u1))
	dz[2] =	-u2 + p4 * u1 * u2
	dz
end

pp2(z, p) = pp2!(similar(z), z, p, 0)
jet  = BK.getJet(pp2; matrixfree=false)

# parameters of the model
par_pp2 = (p1 = 1., p2 = 3., p3 = 5., p4 = 3.)

# initial condition
z0 = zeros(2)

# continuation options
opts_br = ContinuationPar(pMin = 0.1, pMax = 1.0, dsmax = 0.01,
	# options to detect bifurcations
	detectBifurcation = 3, nInversion = 8, maxBisectionSteps = 25,
	# number of eigenvalues
	nev = 2,
	# maximal number of continuation steps
	maxSteps = 1000,
	# parameter theta, see `? continuation`. Setting this to a non
	# default value helps passing the transcritical bifurcation
	theta = 0.3)

diagram = bifurcationdiagram(jet...,
	# initial point and parameter
	z0, par_pp2,
	# specify the continuation parameter
	(@lens _.p1),
	# very important parameter. This specifies the maximum amount of recursion
	# when computing the bifurcation diagram. It means we allow computing branches of branches of branches
	# at most in the present case.
	3,
	(args...) -> setproperties(opts_br; ds = -0.001, dsmax = 0.01, nInversion = 8, detectBifurcation = 3);
	# δp = -0.01,
	recordFromSolution = recordFromSolution,
	verbosity = 0, plot = false)

# branch of the diagram with Hopf point
brH = BK.getBranch(diagram, (2,1)).γ

# newton parameters
optn_po = NewtonPar(tol = 1e-8,  maxIter = 25)

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.1, ds= -0.001, dsmin = 1e-4,
 newtonOptions = (@set optn_po.tol = 1e-8), precisionStability = 1e-2,
 detectBifurcation = 1, saveSolEveryStep=1)

Mt = 101 # number of time sections
	br_po, utrap = continuation(
	jet..., brH, 1, opts_po_cont,
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


plot(diagram); plot!(br_po, label = "Periodic orbits", legend = :bottomright)




