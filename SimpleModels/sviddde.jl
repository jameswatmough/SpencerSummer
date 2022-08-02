using DifferentialEquations
using Plots

# main dde model 
# max's put in to deal with possible negative solutions
# not clear that negatives are a problem, or that the max is the solution
function svid_model(dx,state,h,param,t)
	β, p, d_I, d_V, d_D, tau = param
	S, I, V, D = state
	infection = max(β*S*V)
	delayedinf =max(β*h(p,t-tau)[1]*h(p,t-tau)[3])
	dx[1] = (1-(S+I+D))*S - infection
	dx[2] = delayedinf - d_I*I
	dx[3] = min(p*I) - d_V*V
	dx[4] = d_I*I - d_D*D
end

# variation on above using a nonautonomous influx to build initial conditions
# after sometesting, using x0=[1,0,1,0] and h0 = [1,0,0,0] seems to work just as well
function svid_setics(dx,state,h,param,t)
	β, p, d_I, d_V, d_D, tau = param
	S, I, V, D = state
	inflow = .1*(t<0.1)
	infection = max(β*S*V,0.)
	delayedinf =max(β*h(p,t-tau)[1]*h(p,t-tau)[3] , 0.)
	dx[1] = (1-(S+I+D))*S - infection
	dx[2] = delayedinf - d_I*I
	dx[3] = min(p*I,0.) - d_V*V + inflow
	dx[4] = d_I*I - d_D*D
end


h0(p,t) = [1.,0.,0.,0.]

par_base = (β=5., p=30.0, d_I=2., d_V=20., d_D=1., tau=0.3)

x0 = [1., 0., 1., 0.]

lags = [par_base.tau]

tspan = (0.,10.)

prob = DDEProblem(svid_model,x0,h0,tspan,par_base;constant_lags=lags)

alg = MethodOfSteps(Tsit5())

sol = solve(prob,alg)
