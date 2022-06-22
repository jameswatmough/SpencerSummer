using BifurcationKit, Plots, Setfield, Parameters, LinearAlgebra, ForwardDiff, Revise
const BK = BifurcationKit

function budworm!(dz, z, p, t)
    @unpack Rs, Re, Ks, Ke, P, B = p
    S, E = z
    dz[1] = Rs*S*(1-(S*Ke)/(E*Ks))
    dz[2] = Re*E*(1-E/Ke)-(P*B)/S
    dz
end

budworm(z, p) = budworm!(similar(z), z, p, 0)

dbudworm(z,p) = ForwardDiff.jacobian(x -> budworm(x,p), z)
jet = BK.getJet(budworm, dbudworm)

opts = ContinuationPar(dsmin = 1.e-4, dsmax = 1.e-2, ds = -1.e-4, pMin = -17., pMax = 20., detectBifurcation = 3)
#in order to reverse the direction of the continuation, the sign of ds must be flipped.

params = (Rs = 100., Re = 95., Ks = 1.1, Ke = 1.2, P = 1., B = 1.)

br, = continuation(budworm, dbudworm, [0.1,0.1], params, (@lens _.B), opts; printSolution = (x,p) -> (S = x[1], E = x[2]))

plot(br)