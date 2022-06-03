
using DifferentialEquations, Plots

# Simple IGP ode model with saturating numerical responses

function igpode!(dx,x,p,t)
	c1,c2,c3,e1,e2,e3,k2,k3,m1,m2 = p
	dx[1] = x[1]*((1 - x[1]) - c1*x[2] - c2*x[3])
	dx[2] = x[2]*(e1*c1*x[1] - c3*x[3]/(1+k3*x[2]) - m1)
	dx[3] = x[3]*(e2*c2*x[1]/(1+k2*x[1]) + e3*c3*x[2]/(1+k3*x[2]) - m2)
end

function Figpode(x,p)
	c1,c2,c3,e1,e2,e3,k2,k3,m1,m2 = p
	return([
    x[1]*((1 - x[1]) - c1*x[2] - c2*x[3]),
    x[2]*(e1*c1*x[1] - c3*x[3]/(1+k3*x[2]) - m1),
    x[3]*(e2*c2*x[1]/(1+k2*x[1]) + e3*c3*x[2]/(1+k3*x[2]) - m2)
	 ])
end

function Jigpode(x,p)
	c1,c2,c3,e1,e2,e3,k2,k3,m1,m2 = p
	return([
    [ 1 - 2*x[1] - c1*x[2] - c2*x[3] ;;
		  -c1*x[1];;
			-c2*x[3]];
		[ e1*c1*x[2] ;;
		  e1*c1*x[1] - c3*x[3]*x[2]/(1+k3*x[2])^2 - m1 ;;
			-c3*x[2]/(1+k3*x[2]) ];
		[ x[3]*e2*c2*x[1]/(1+k2*x[1])^2 ;;  
		  x[3]*e3*c3*x[2]/(1+k3*x[2])^2 ;;
			e2*c2*x[1]/(1+k2*x[1]) + e3*c3*x[2]/(1+k3*x[2]) - m2 ];
		])
end

# parameter values for torus solution,
# but probably better to start with k2 = k3 =0
r=1;
k1=1;
c1=α2=1;
c2=α2=0.048;
e1=σ1=50;
c3=1;
m1=μ1=10;
e2=σ2=150;
e3=σ3=0.4;
m2=μ2=1;
k3=β3=0.2;
k2=β2=5;

p = [c1,c2,c3,e1,e2,e3,k2,k3,m1,m2]

# setting up the torus
x0 = [0.270, 0.105, 1.82]
tspan = (0.0,1000.0)

prob = ODEProblem(igpode!,x0,tspan,p)
sol = solve(prob)

anim = @animate for i in 10:2:tspan[2]
	plot(sol,vars=(1,2,3),xlabel="resource",ylabel="prey",zlabel="predator",legend=false)
	plot!(sol,vars=(1,2,3),tspan=(i-10,i),linecolor=:red)
       end

# saved this, using a variation of
gif(anim,"testanimcol.gif",fps=30)

#cool picture
solmid = solve(remake(prob;u0=[0.2, 1.05, .3875],tspan=(0.,3000)))
plot(solmid,vars=(1,2,3),tspan=(0,1080),lc=:red)
plot!(sol,vars=(1,2,3),lc=:blue,xlabel="resource",ylabel="prey",zlabel="predator",legend=false)
plot!(solmid,vars=(1,2,3),tspan=(1180,2000),lc=:red)
