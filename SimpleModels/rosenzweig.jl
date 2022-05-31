using DifferentialEquations, Plots

function rosenzweig!(dx,x,p,t)
	b,d,K,a,w,c,m = p
	dx[1] = (b-d)*x[1]*(1 - x[1]/K) - a*w*x[1]*y[1]/(w+x[1])
	dx[2] = c*a*w*x[1]*y[1]/(w+x[1]) - m*y[1]
end

# sample parameter values 
b=2
d=1
K=1000
a=0.005
w=400
c=2
m=2

p = [b,d,K,a,w,c,m]
x0 = [500.,500.]
tspan = (0.,100.)

# the equilibrium prey population is
Nstar = m*w/(c*a*w-m)  
# a hopf bifurcation occurs as K is increased through
Khopf = 2*Nstar+w
