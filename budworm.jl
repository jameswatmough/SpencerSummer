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

dbudworm(z, p) = ForwardDiff.jacobian(x -> budworm(x,p), z)
jet = BK.getJet(budworm, dbudworm)

br = []

dssign = -1.
param = 18.
init = 1.2
bifurccount = 0

function cont(s, param, init)
    
    plot()
    
    for i in 1:15
        
        println(s)
        opts = ContinuationPar(pMin = param - 10., pMax = param + 10., theta = 1.,
            dsmin = 1.0e-3, ds = s*4.0e-3, dsmax = 5.0e-3,
            detectBifurcation = 3, 
            nInversion = 8, maxBisectionSteps = 25, nev = 3)
        #in order to reverse the direction of the continuation, the sign of ds must be flipped.

        params = (Rs = 100., Re = 95., Ks = 1.1, Ke = 1.2, P = 1., B = param)

        br, = continuation(budworm, dbudworm, [init,init], params, (@lens _.B), opts; printSolution = (x,p) -> (S = x[1], E = x[2]))

        for n in 1:length(br)-1
            if br[n+1][7] != br[n][7]
                println("bifurcation detected! (x = ", br[n+1][1], ", B = ", br[n+1][2], ")")
                s = -1*s
                global bpind
                bpind = n+1
            else
                global bpind
                bpind = n+1
            end
        end
        
        bfinfo = br[bpind]
        
        init = bfinfo[1]
        param = bfinfo[2]
        plot!(br, plotfold = true, linecolor = :blue)
        
        println("continuation trial ", i, " started at (x = ", br[1][1], ", B = ", br[1][2], ")")
        println("and terminated at (x = ", br[length(br)][1], ", B = ", br[length(br)][2], ").")

    end
    
    plot!()
    
end


cont(1.,0.1,0.1)



