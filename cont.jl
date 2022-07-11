using BifurcationKit, Plots, Setfield, Parameters, LinearAlgebra, ForwardDiff, Revise
const BK = BifurcationKit

function rosenzweig!(dz,z,p,t)
    @unpack b,d,K,a,w,c,m = p
    x1, y1 = z
    dz[1] = (b-d)*x1*(1 - x1/K) - a*w*x1*y1/(w+x1)
    dz[2] = c*a*w*x1*y1/(w+x1) - m*y1
    dz
end

rosenzweig(z, p) = rosenzweig!(similar(z), z, p, 0)
drosenzweig(z, p) = ForwardDiff.jacobian(x -> rosenzweig(x,p), z)

jet = BK.getJet(rosenzweig, drosenzweig)

z0 = [410.,410.]
z1 = [500.,500.]

printDetails = false
#setting this to 'true' prints every bifurcation detection as well as the start and end of every iteration
turningPointDetection = false
#what this should be modified to do is instead pick a threshold for sloping, and change direction accordingly, instead of hoping that a
#bifurcation always results in a change in the sign of the differential.


function cont!(s, param, init, iter)

    for i in 1:iter
        
        if printDetails == true
            println(s)
        end
        opts = ContinuationPar(pMin = param - 1000., pMax = param + 1000., theta = 1.,
            dsmin = 1.0e-3, ds = s*4.0e-3, dsmax = 5.0e-3,
            detectBifurcation = 3, 
            nInversion = 8, maxBisectionSteps = 25, nev = 3)

        params = (b=2., d=1., K=param, a=0.005, w=400., c=2., m=2.)

        br, = continuation(rosenzweig, drosenzweig, init, params, (@lens _.K), opts; printSolution = (x,p) -> (x1 = x[1], x2 = x[2]))

        for n in 1:length(br)-1
            if br[n+1][7] != br[n][7] || br[n+1][8] != br[n][8]
                if printDetails == true
                    println("bifurcation detected! (x = ", br[n+1][1], ", B = ", br[n+1][2], ")")
                end
                if turningPointDetection == true
                    s = -1*s
                end
                global bpind
                bpind = n+1
            else
                global bpind
                bpind = n+1
            end
        end
        
        bfinfo = br[bpind]
        
        init = [bfinfo[1],bfinfo[1]]
        param = bfinfo[2]
        plot!(br, plotfold = true, linecolor = :blue)
        
        if printDetails == true
            println("continuation trial ", i, " started at (x = ", br[1][1], ", K = ", br[1][2], ")")
            println("and terminated at (x = ", br[length(br)][1], ", K = ", br[length(br)][2], ").")
        end
    end
    
    plot!(legend = :topleft)
    
end



function cont(s, param, init, iter)
    
    plot()
    cont!(s,param,init,iter)
    
end



cont(-1.,500.,z0,200)
cont!(1.,1000.,z1,70)

#first argument: initial direction (-1 (left) or 1 (right))
#second arguent: initial parameter value
#third argument: initial condition estimate
#fourth argument: # of iterations