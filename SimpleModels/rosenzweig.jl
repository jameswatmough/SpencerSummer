using BifurcationKit, Plots, Setfield, Parameters, LinearAlgebra, ForwardDiff, Revise
const BK = BifurcationKit

# sup norm
norminf(x) = norm(x, Inf)

function rosenzweig!(dx,x,p,t)
	@unpack μ,γ,κ = p
	N,P = x
	dx[1] = N*(1 - N/κ) - N*P/(1.0+N)
	dx[2] = μ*P*(γ*N/(1.0+N) - 1.0)
	return(dx)
end

rosenzweig(z, p) = rosenzweig!(similar(z), z, p, 0)

# sample parameter values 
par_base = (μ=1., γ=2.0, κ=0.5)

x0 = [par_base.κ,0.]
tspan = (0.,100.)

prob = BifurcationProblem(
					rosenzweig, x0, par_base, (@lens _.κ);
					recordFromSolution = (x, p) -> (prey = x[1], predator = x[2])
			 )


# continuation options
opts_br = ContinuationPar(pMin = 0.0, pMax = 5.,
	# parameters to have a smooth result
	ds = 0.04, dsmax = 0.05,
	# this is to detect bifurcation points precisely with bisection
	detectBifurcation = 3,
	# Optional: bisection options for locating bifurcations
	nInversion = 8, maxBisectionSteps = 25, nev = 3)

br = continuation(prob,PALC(tangent=Bordered()), opts_br)

# found a transcritcal bifurcation, continue that
br2 = continuation(br,1,opts_br)

# found a Hopf bifurcation, plot the periodic orbits
opt_po = NewtonPar(tol = 1e-10, verbose = true, maxIter = 15)
opts_po_cont = ContinuationPar(
	dsmin = 0.001, dsmax = 0.04, ds = 0.03, pMax = 5., 
	maxSteps = 200, newtonOptions = opt_po, saveSolEveryStep = 2,
	plotEveryStep = 1, nev = 11, tolStability = 1e-6,
	detectBifurcation = 3, dsminBisection = 1e-6, maxBisectionSteps = 15, tolBisectionEigenvalue = 0.)

brH = continuation(
		br2,2, opts_po_cont,
		PeriodicOrbitTrapProblem(M=51);
		δp = 0.01, ampfactor=1,
		verobsity = 3, plot = false
		)

# now to extract the maximums and minimums of the periodic orbits
# presumably this can be cleaned up nicely
brH_len = length(brH.γ.sol)
brH_p = [brH.γ.sol[i].p for i in 1:brH_len]
# x component of solution is the flattened solution
# reshape(brH.γ.sol[i].x[1:end-1],2,:) gives the solution in matrix form
# not clear in documentation what the last entry of x is
# we'll use just the max/min of prey density along the solution
brH_max = [maximum(brH.γ.sol[i].x[1:2:end-1]) for i in 1:brH_len]
brH_min = [minimum(brH.γ.sol[i].x[1:2:end-1]) for i in 1:brH_len]
plot(br,br2)
plot!(brH_p,brH_max)
plot!(brH_p,brH_min)

