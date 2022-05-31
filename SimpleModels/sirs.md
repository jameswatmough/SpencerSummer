# Quick and dirty look at the SIRS model 

\begin{align*}
	S' &= -β SI + ω R\\
	I' &=  β SI - γ I\\
	R' &=  γ I  - ω R
\end{align*}

Population remains constant if we neglect disease induced mortality.
Assume $S+I+R=1$.

\begin{align*}
	S' &= -β SI + ω (1-S-I)\\
	I' &=  β SI - γ I
\end{align*}

Endemic equilibrium found by solving 
\begin{align*}
	0 &= -β SI + ω (1-S-I)\\
	0 &=  β SI - γ I
\end{align*}

This leads to 
\begin{align*}
	S_e &= \frac{γ}{β} \\ 
	I_e &= (1-γ/β)\frac{ω}{ω+γ} = (1-γ/β)\left( 1 - \frac{γ}{ω+γ} \right)
\end{align*}

This is stable, with possible decaying oscillations.

Sensitivity of $I_e$ to changes in $ω$.

$$\frac{ω}{I_e}\frac{\partial I_e}{\partial ω} 
   = \left( (1-γ/β)\frac{1}{ω+γ} \right)^{-1} \left(-(1-γ/β)\frac{γ}{(ω+γ)^2}\right)
   = -\frac{γ}{(ω+γ)}$$

For SARS-CoV-2 and CoViD-19, a ball park estimate for the ratio of the duration of immunity to the duration of infection  is somewhere between 10:1 and 50:1.  With Ro = β/γ = 7, this translates to equilibrium prevalences of 0.08 and 0.02 respectively.  Thus the equilibrium prevalence and (presumably) peak prevalence of any subsequent waves of decaying oscillations are very sensitive to the duration of immunity.  A 10% change in duration of immunity induces an 9% to 10% change in peak and equilibrium prevalence.



