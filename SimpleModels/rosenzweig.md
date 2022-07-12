# The Rosenzweig McArthur model

\begin{align}
	N' &= (b-d)N(1 - \frac{N}{K}) - \frac{awNP}{w+N}, \\
	P' &= \frac{cawNP}{w+N} - mP
\end{align}

We are interested in examining the effect of the carrying capacity, $K$, on the dynamics, hence suitable rescalings are $x = N/w$, $y = ((b-d)/aw)P$ and a timescale $1/(b-d)$.   This results in the equations

\begin{align}
	x' &= x(1 - \frac{x}{κ}) - \frac{xy}{1+x}, \\
	y' &= μ\left( \frac{γxy}{1+x} - y \right)
\end{align}

with $γ = caw/m$ and $μ = m/(b-d)$.

Note this is equivalent to setting $w=1$, $a=1$, $b=d+1$, $K=κ$, $m=μ$ and $c$ as $γμ$.

There are three possible equilibria:
(1) the trivial equilibrium, $(0,0)$, 
(2) the prey-only equilibrium, $(κ,0)$, and
(3) an interior equilibrium, $(\bar{x},\bar{y})$, with

$\bar{x} = 1/(γ-1)$ and $\bar{y} = (1+\bar{x})(1-\bar{x}/κ)$

It is convenient to introduce functions $f$ and $g$ so that we can write the system as 

\begin{align}
	x' &= xf(x,y), \\
	y' &= yg(x,y)
\end{align}

The Jacobian at $(x,y)$ is 
$$\begin{pmatrix}
 f + xf_x & xf_y \\
 yg_x & g + yg_y
\end{pmatrix}$$


\begin{align*}
	f(x,y) &= 1 - \frac{x}{κ} - \frac{y}{1+x}, \\
	f_x(x,y) &=  - \frac{1}{κ} + \frac{y}{(1+x)^2}, \\
	f_y(x,y) &= - \frac{1}{1+x}, \\
	g(x,y) &= μγ\frac{x}{1+x} - μ \\
	g_x(x,y) &= μγ\frac{1}{(1+x)^2} \\
	g_y(x,y) &= 0
\end{align*}


For the trivial equilibrium, 
$$J(0,0) = \begin{pmatrix}  1 & 0 \\ 0 & -μ \end{pmatrix}$$
Always a saddle whose stable and unstable manifolds are the two axes.

For the prey-only equilibrium,
$$J(κ,0) = \begin{pmatrix}  -1 & -κ/(1+κ) \\ 0 & μγκ/(1+κ)-μ \end{pmatrix}$$
This is stable if either $γ<μ$ or $κ < 1/(γ-1)$
and unstable if $γ>μ$ and $κ < 1/(γ-1)$.
Note that these conditions correspond to positivity of the interior equilibrium and a transcritical bifurcation of the two equilibria.

For the interior equilibria,

$$J(\bar{x},\bar{y}) = \begin{pmatrix}  
    \bar{x}f_x & \bar{x}f_y \\ 
    \bar{y}g_x & 0
	\end{pmatrix}$$

which has eigenvalues
$$\frac{xf_x}{2} 
     \pm \sqrt{ \frac{(xf_x)^2}{4} + \bar{x}\bar{y}f_yg_x }$$

since $f_yg_x < 0$ the eigenvalues are pure imaginary if $f_x=0$. 
Thus a Hopf bifurcation occurs as $f_x$ changes sign.

\begin{gather}
	f_x < 0  \\
- \frac{1}{κ} + \frac{y}{(1+x)^2} < 0 \\
 κy < (1+x)^2  \\
 (κ-\bar{x}) < 1+\bar{x} \\
 κ < 1+2\bar{x} \\
	κ < \frac{γ+1}{γ-1} \\
	Kb < \frac{caw+m}{caw-m}
\end{gather}

Fixing $γ>1$ and setting $\bar{x} = 1/(γ-1)$ we can follow the two interesting branches of equilibria as $κ$ increases.  
$(κ,0)$ which switch from stable to unstable as $κ$ increases through $\bar{x}$,
and $(\bar{x}, \bar{y}(κ))$ which bifurcates from the first branch at $κ=\bar{x}$ and switches from stable to unstable via a Hopf bifurcation as $κ$ increases through $1+2\bar{x} = (1+γ)/(1-γ)$.

The simplest plot is probably a plot of the predator component ($y$ or $P$), versus $κ$.
