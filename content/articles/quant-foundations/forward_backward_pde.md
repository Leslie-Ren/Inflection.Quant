---
title: "Forward and Backward Kolmogorov PDEs: An Intuitive Look at Their Duality"
date: 2026-07-14
draft: false
math: true
tags: ["kolmogorov equations", "fokker-planck", "duality", "local volatility", "dupire"]
---

## Why This Matters

Most practitioners have seen the Black-Scholes PDE. It is a backward equation, closely related to the backward Kolmogorov PDE: fix a payoff at maturity, and the PDE propagates its value back to today.

There is also a forward equation, which may be less familiar. The backward equation has current spot and current time as its variables and takes the payoff as a terminal condition. The forward equation instead propagates a probability density forward from today, with the terminal value of the underlying and the maturity as its variables. This is the Fokker-Planck equation.

My interest in the forward equation comes from a plan to write about the local volatility model, where the Dupire calibration formula rests on Fokker-Planck. Rather than compress the forward equation into a few lines of that article, I want to treat it properly first. This article builds an intuitive picture of the forward and backward PDE and the duality that links them. Rather than open with the derivation, we start from a simple case where time and space are both discrete.

## Discrete Time, Discrete Space

Start with a finite Markov chain $X_n$. States are labeled $1, \dots, m$, and one step is described by a transition matrix $P$, where $P_{ij}$ is the probability of moving from state $i$ to state $j$ over one period. Each row sums to one. That a single matrix suffices is the Markov property: where the chain goes next depends on where it is now and nothing else, so its history can be discarded once the current state is known.

As a concrete example, take three weather states, sun, cloud, and rain, with $P$ the one-day transition probability matrix:

$$
P = \begin{array}{cc}
 & \begin{array}{ccc} \text{sun} & \text{cloud} & \text{rain} \end{array} \\
\begin{array}{c} \text{sun} \\ \text{cloud} \\ \text{rain} \end{array} &
\begin{pmatrix} 0.6 & 0.3 & 0.1 \\ 0.4 & 0.4 & 0.2 \\ 0.2 & 0.5 & 0.3 \end{pmatrix}
\end{array}.
$$

Row one says a sunny day is followed by sun, cloud, or rain with probabilities $0.6$, $0.3$, $0.1$.

Two natural questions can be asked of this chain. Given where it is now, what is the expected value of something that depends on where it ends? And given a distribution over where it is now, where is it likely to be later? The first is answered by a value function stepping backward from the terminal date, the second by a distribution stepping forward from today. The two equations run in opposite directions.

**The backward equation.** Assign a payoff to each state, so that $g(i)$ is the amount received if the chain ends in state $i$, and write $g$ for the column vector of these payoffs. The value at step $n$ is the expected payoff given where the chain is now,

$$
u_n(i) = \mathbb{E}\big[g(X_N) \mid X_n = i\big].
$$

Conditioning on the first move,

$$
u_n(i) = \sum_j P_{ij}\, u_{n+1}(j), \qquad \text{that is} \qquad u_n = P u_{n+1}.
$$

Suppose the payoff pays $10$ on a sunny terminal day and $2$ otherwise, so $u_N = g = (10, 2, 2)^\top$. One step back,

$$
u_{N-1} = P u_N = (6.8,\ 5.2,\ 3.6)^\top,
$$

the expected payoff at each state today.

**The forward equation.** Fix a distribution over states today, as a column vector $\pi_n$, so that $\pi_n(i) = \mathbb{P}(X_n = i)$. What is the distribution one step later? Summing over the intermediate state,

$$
\pi_{n+1}(j) = \sum_i P_{ij}\, \pi_n(i), \qquad \text{that is} \qquad \pi_{n+1} = P^\top \pi_n.
$$

Starting from a sunny day today, $\pi_0 = (1, 0, 0)^\top$, one step gives $\pi_1 = (0.6, 0.3, 0.1)^\top$, the first row of $P$, as it must.

**The duality.** The two updates sit side by side: the value steps through $P$,

$$
u_n = P\, u_{n+1},
$$

while the distribution steps through its transpose,

$$
\pi_{n+1} = P^\top \pi_n.
$$

What ties them together is the expected payoff, which can be computed two ways. Step the value back to today and average it under today's distribution:

$$
\pi_n^\top \big(P\, u_{n+1}\big).
$$

Or push today's distribution forward and average the next step's value under it:

$$
\big(P^\top \pi_n\big)^\top u_{n+1}.
$$

The two differ only in how the product is grouped and transposed:

$$
\pi_n^\top \big(P\, u_{n+1}\big) = \big(P^\top \pi_n\big)^\top u_{n+1}.
$$

This is the duality relationship. What follows from it is an invariance. Substituting the two update equations, each side becomes a product at a single step,

$$
\pi_n^\top u_n = \pi_{n+1}^\top u_{n+1},
$$

and this chains from one step to the next,

$$
\pi_0^\top u_0 = \pi_1^\top u_1 = \cdots = \pi_N^\top u_N.
$$

The expected payoff is the same number at every step, so it can be computed wherever it is convenient, at step 0, at step 5, or at maturity. That is what makes the forward and backward equations interchangeable as pricing methods rather than merely algebraically related.

## Continuous Time, Discrete Space

The Kolmogorov equations live in continuous time and continuous space. Rather than take both limits at once, we take them one at a time, starting with time. Space stays discrete, so the chain still moves among finitely many states, but transitions can now happen at any moment.

In the discrete chain, $P$ was the transition matrix over a fixed period, a day in the weather example. In continuous time there is no natural period to quote probabilities over. Over a short interval, $P$ must be close to the identity, since there is little time to change states. What remains is to see how the chance of leaving a state grows with the length of the interval.

Write $f(h)$ for the probability of having left state $i$ by time $h$. Over a window of length $2h$, the chain either leaves during the first half, or it does not and then leaves during the second. By the memoryless property of the Markov chain, the chain does not track how long it has been sitting in state $i$, so the second half is governed by the same $f(h)$ as the first, giving

$$
f(2h) = \underbrace{f(h)}_{\text{leaves in first half}} + \underbrace{\big(1 - f(h)\big)}_{\text{survives first half}}\underbrace{f(h)}_{\text{leaves in second}} = 2f(h) - \underbrace{f(h)^2}_{\text{two departures}}.
$$

For short $h$ the quadratic term is negligible, so $f(2h) \approx 2 f(h)$: doubling the window doubles the chance. The same argument applied to any subdivision gives $f(h) \approx q\,h$ for a constant rate $q$. So the chance of leaving is proportional to the length of the interval, and

$$
P = I + Q\, dt.
$$

The entries of $Q$ are exactly these rates: off the diagonal, $q_{ij} \geq 0$ is the rate of moving from $i$ to $j$, and the diagonal is minus the total rate of leaving each state, $Q_{ii} = -\sum_{j \neq i} q_{ij}$, so that $P_{ii} = 1 + Q_{ii}\, dt$ is the chance of staying. The rows of $Q$ sum to zero as the rows of $P$ sum to one.

**The backward equation.** Substituting into the backward recursion, $u(t) = P\, u(t + dt) = (I + Q\,dt)\, u(t+dt)$, rearranging and taking $dt \to 0$ gives

$$
\frac{du}{dt} = -Q u.
$$

**The forward equation.** With $\pi(t+dt) = P^\top \pi(t) = (I + Q^\top dt)\,\pi(t)$, the same limit gives

$$
\frac{d\pi}{dt} = Q^\top \pi.
$$

Notice that $Q$ is the sole input to both equations. Fix the rates and the whole process is determined, which is why $Q$ is called the generator.

**The duality.** In the discrete case the duality was:

$$
\pi^\top \big(P\, u\big) = \big(P^\top \pi\big)^\top u.
$$

Substituting $P = I + Q\,dt$ leaves the same identity for $Q$:

$$
\pi^\top \big(Q\, u\big) = \big(Q^\top \pi\big)^\top u.
$$

This is the duality in continuous time, and it has a physical interpretation.

**The left side: the value moving, distribution fixed.** Writing out the generator against a value vector,

$$
(Qu)_i = \sum_{j \neq i} q_{ij}\,(u_j - u_i),
$$

the expected rate of value gain from state $i$: each possible jump weighted by its rate and by the value it gains. Then $\pi^\top(Qu)$ averages that over where the chain is now.

**The right side: the distribution moving, value fixed.** $Q^\top\pi$ is the rate at which probability accumulates in each state. Pairing it with $u$ gives the rate at which expected payoff accrues as probability shifts between states.

The consequence of the duality is that the expected payoff does not move at all. By the product rule and the two evolution equations,

$$
\frac{d}{dt}\big(\pi^\top u\big) = \left(\frac{d\pi}{dt}\right)^\top u + \pi^\top \frac{du}{dt} = \big(Q^\top \pi\big)^\top u + \pi^\top \big(-Q u\big) = 0,
$$

This is the invariance in continuous time. In the discrete case the expected payoff was the same number at every step; here it is the same number at every instant.

## Continuous Time, Continuous Space

Finally, let the state space also become continuous. The chain $X_n$ becomes a process $X_t$ moving on a continuum, the distribution $\pi$ becomes a density $p(x,t)$ with $p(x,t)\,dx$ the probability of being near $x$ at time $t$, and the sums over states become integrals.

**The backward equation.** The value is defined as a conditional expectation,

$$
u(x,t) = \mathbb{E}\big[g(X_T) \mid X_t = x\big],
$$

and conditioning on the move over $[t, t+dt]$ leaves it unchanged,

$$
u(x,t) = \mathbb{E}\big[u(X_{t+dt}, t+dt) \mid X_t = x\big],
$$

the continuous version of $u_n = P u_{n+1}$, with the sum over destinations replaced by an expectation. The expected change in $u$ over $[t, t+dt]$ is therefore zero. Two things move over that interval, the clock and the state, and they can be separated,

$$
\mathbb{E}\big[u(X_{t+dt}, t+dt)\big] - u(x,t) = \underbrace{\mathbb{E}\big[u(X_{t+dt}, t+dt) - u(X_{t+dt}, t)\big]}_{\text{the clock advances}} + \underbrace{\mathbb{E}\big[u(X_{t+dt}, t)\big] - u(x,t)}_{\text{the state moves}},
$$

all expectations conditional on $X_t = x$. Divide by $dt$ and let $dt \to 0$. The first bracket holds the state fixed across its two terms, so it is $\partial_t u$; strictly it is $\mathbb{E}\big[\partial_t u(X_{t+dt}, t)\big]$, but $X_{t+dt} \to x$ as $dt \to 0$, so the evaluation point collapses onto $x$. The second holds time fixed at $t$ and gives the expected rate of change of the value from the state alone. The second bracket, divided by $dt$, can be denoted $Lu$: the instantaneous change of the value from the state moving. Here $L$ is the generator, the analogue of $Q$ in the discrete state case. Since the left side is zero, the backward equation takes the same form as in the discrete case,

$$
\partial_t u = -L u, \qquad L u(x,t) = \lim_{dt \to 0} \frac{\mathbb{E}\big[u(X_{t+dt}, t) \mid X_t = x\big] - u(x,t)}{dt}.
$$

This is the Kolmogorov backward equation, and it holds for any continuous-time Markov process. But it is not yet a PDE: $L$ is defined by a limit, not by derivatives, and nothing so far says what operator it actually is. To get a PDE, $L$ has to be computed, and that requires committing to how the process moves.

**The backward PDE.** Take the process to be a diffusion,

$$
dX_t = \mu(X_t)\, dt + \sigma(X_t)\, dW_t.
$$

Apply Itô's lemma to find $du$, with $t$ held fixed so that only $X$ moves — not the full lemma, so it does not carry the $\partial_t u\,dt$ term:

$$
du(X_t) = \Big(\mu\, \partial_x u + \tfrac{1}{2}\sigma^2 \partial_{xx} u\Big) dt + \sigma\, \partial_x u\, dW_t,
$$

and taking expectations kills the $dW$ term, leaving

$$
L = \mu(x)\, \partial_x + \tfrac{1}{2}\sigma^2(x)\, \partial_{xx}.
$$

Substituting into $\partial_t u = -Lu$ gives the backward equation as a PDE,

$$
\partial_t u + \mu(x)\, \partial_x u + \tfrac{1}{2}\sigma^2(x)\, \partial_{xx} u = 0.
$$

**The forward equation.** The forward equation describes how the density evolves, $\partial_t p$. If the underlying process has a known transition density, such as a diffusion with constant drift and volatility, then $p$ is available in closed form and taking its $t$-derivative is straightforward. Here we want the general case, where the underlying process may not have a closed-form density. There are two routes to tackle the problem. The first works from the definition, tackling $\partial_t p$ directly by working out how the density evolves over a short time interval; it reaches the forward equation without reference to the backward one. This approach is more involved and is laid out in [Appendix A](#appendix-a-deriving-the-forward-equation-directly-from-the-densitys-evolution). The second takes advantage of the duality: we already have the generator $L$ from the backward equation, and the duality is what converts it into the forward equation. We follow the second route here.

The duality route mirrors the discrete case. There the forward equation was the transpose of the backward one, $d\pi/dt = Q^\top\pi$ against $du/dt = -Qu$, and a matrix has a transpose immediately to hand. A differential operator does not, so what we need is the analogue of that transpose: the adjoint of $L$.

Start from what the transpose meant discretely. The duality established above,

$$
\pi^\top \big(Q\, u\big) = \big(Q^\top \pi\big)^\top u,
$$

is the statement that $Q$ can be moved from one side of the pairing $\pi^\top u = \sum_i \pi_i u_i$ to the other, at the cost of transposing it. Written in coordinates,

$$
\sum_i (Q u)_i\, \pi_i = \sum_i u_i\, (Q^\top \pi)_i.
$$

This is not a result to prove but the definition of the transpose: $Q^\top$ is whatever reproduces the left side with the operator acting on $\pi$ instead of $u$.

The continuum uses the same definition with the sum replaced by an integral and the density $p$ in the role of $\pi$. So $L^*$ is defined as the operator satisfying

$$
\int (L u)\, p \, dx = \int u\, (L^* p)\, dx
$$

for all suitable $u$ and $p$. This is the continuous duality, and $L^*$ is called the adjoint of $L$. With it, the forward equation is

$$
\partial_t p = L^* p,
$$

the exact analogue of $d\pi/dt = Q^\top\pi$ from the discrete case, with $L^*$ in place of $Q^\top$. This must hold because the expected payoff $\int u\,p\,dx$ is constant in time by definition, $\tfrac{d}{dt}\int u\,p\,dx = 0$, and pairing that invariance against the backward equation leaves $L^*$ as the only operator that can drive $p$. Here $L^*$ is still defined by a property rather than by derivatives; making it explicit is what turns the equation into a PDE.

To find $L^*$ explicitly, take $L = \mu(x)\,\partial_x + \tfrac12\sigma^2(x)\,\partial_{xx}$ from the backward equation and substitute it into the left side of the adjoint identity. Recovering $L^*$ is then a matter of moving the derivatives off $u$ and onto $p$ by integration by parts.

Take the drift term first. One integration by parts gives

$$
\int \mu(x)\, \partial_x u\, p(x)\, dx = -\int u(x)\, \partial_x\big(\mu(x) p(x)\big)\, dx.
$$

For the diffusion term, integrate by parts twice:

$$
\int \tfrac{1}{2}\sigma^2(x)\, \partial_{xx} u\, p(x)\, dx = \int u(x)\, \tfrac{1}{2}\partial_{xx}\big(\sigma^2(x) p(x)\big)\, dx.
$$

so the forward equation, written out, becomes

$$
\partial_t p = L^* p = -\partial_x\big(\mu(x) p\big) + \tfrac{1}{2}\partial_{xx}\big(\sigma^2(x) p\big).
$$


### Summary


| | Discrete time, discrete space | Continuous time, discrete space | Continuous time, continuous space |
|---|---|---|---|
| **Backward** | $u_n = P u_{n+1}$ | $\tfrac{du}{dt} = -Q u$ | $\partial_t u = -L u$ |
| **Forward** | $\pi_{n+1} = P^\top \pi_n$ | $\tfrac{d\pi}{dt} = Q^\top \pi$ | $\partial_t p = L^* p$ |
| **Duality** | $\pi^\top(P u) = (P^\top\pi)^\top u$ | $\pi^\top(Q u) = (Q^\top\pi)^\top u$ | $\int (Lu)\,p\,dx = \int u\,(L^*p)\,dx$ |
| **Invariance** | $\pi_n^\top u_n$ constant in $n$ | $\pi^\top u$ constant in $t$ | $\int u\,p\,dx$ constant in $t$ |

The transition object turns into a generator and then into a differential operator, its transpose turns into an adjoint, and the pairing that stays constant turns a sum into an integral. Backward and forward are adjoint throughout: $P$ against $P^\top$, $Q$ against $Q^\top$, and $L$ against $L^*$.

## Application in Finance

The most recognizable application is option pricing. Under the risk-neutral measure, a stock follows $dS_t = r S_t\, dt + \sigma S_t\, dW_t$, so the generator is

$$
L = r S\, \partial_S + \tfrac{1}{2}\sigma^2 S^2\, \partial_{SS}.
$$

The value of a claim with payoff $g(S_T)$ is $V(S, t) = e^{-r(T-t)}\,\mathbb{E}[g(S_T) \mid S_t = S]$, and it satisfies the backward equation with a discounting term,

$$
\partial_t V + L V - r V = 0.
$$

This is the Black-Scholes PDE: Kolmogorov backward under the pricing measure, with the extra $-rV$ accounting for discounting. Its variables are current spot and current time, and the terminal condition is the payoff. That the pricing expectation solves this backward PDE is the content of Feynman-Kac.

When is the forward equation the one you want? Whenever the object of interest is the density itself, as a function of where and when the underlying lands. Applying the adjoint recipe to $L$ gives

$$
L^* p = -\partial_S\big(r S\, p\big) + \tfrac{1}{2}\partial_{SS}\big(\sigma^2 S^2\, p\big),
$$

and the risk-neutral transition density $p(S, T \mid S_0, t_0)$ satisfies the Fokker-Planck equation $\partial_T p = L^* p$ in the forward variables $S$ and $T$. Whenever the density is the object you need, the forward equation is the natural direction: a single forward solve produces $p(S,T)$ across all terminal states at once, so vanilla prices at every strike and maturity follow by integration against payoffs, without a separate backward solve per claim. That efficiency is what makes it the workhorse of calibration and of density recovery for risk.

A prominent use is local volatility calibration. The call price is the density integrated against the call payoff,

$$
C(K, T) = e^{-rT}\int_K^\infty (S - K)\, p(S, T)\, dS,
$$

and because $p$ evolves forward in $T$ while the payoff depends on $K$, differentiating in $K$ and $T$ turns Fokker-Planck into a relation among the strike and maturity derivatives of $C$. That relation is the Dupire formula, and it is available because the density obeys a forward equation in exactly the variables the option surface is quoted in. The full construction is the subject of a follow up article on local vol.

## Appendix A: Deriving the Forward Equation Directly from the Density's Evolution

The starting point is the continuous analogue of $\pi_{n+1} = P^\top \pi_n$: the density at $t + dt$ is the density at $t$ pushed through one short transition,

$$
p(y, t+dt) = \int p(y \mid x, dt)\, p(x, t)\, dx,
$$

where $p(y \mid x, dt)$ is the probability density of landing near $y$ having started at $x$ an interval $dt$ earlier. As stated in the main text, if the transition density were known in closed form we would substitute it and be done. For a general diffusion it is not known. What we know about the transition is its first two moments:

$$
\mathbb{E}[\Delta X \mid X_t = x] = \int (y-x)\, p(y\mid x, dt)\, dy = \mu(x)\,dt + o(dt), \qquad \mathbb{E}\big[(\Delta X)^2 \mid X_t = x\big] = \int (y-x)^2\, p(y\mid x, dt)\, dy = \sigma^2(x)\,dt + o(dt),
$$

where $\Delta X = X_{t+dt} - x$, with all higher moments $o(dt)$.

So we have the unconditional density expressed as an integration over $x$, and two conditional moments integrated over $y$. **Taylor expansion** is the way to use them.

Since Taylor works for any function, let $\phi$ be an arbitrary smooth function. Expanding $\phi(x + \Delta X)$ in powers of the step and integrating against the transition gives

$$
\mathbb{E}\big[\phi(x + \Delta X) \mid X_t = x\big] = \phi(x) + \partial_x \phi(x)\,\mathbb{E}[\Delta X \mid X_t = x] + \tfrac{1}{2}\partial_{xx}\phi(x)\,\mathbb{E}\big[(\Delta X)^2 \mid X_t = x\big] + o(dt).
$$

For a continuous diffusion process, every term past the second carries a higher moment that is $o(dt)$; note that for a jump process this is not the case, and a different equation is reached. We focus on the continuous case here. Substituting the two moments,

$$
\mathbb{E}\big[\phi(x + \Delta X) \mid X_t = x\big] = \phi(x) + \mu(x)\,\partial_x \phi(x)\,dt + \tfrac{1}{2}\sigma^2(x)\,\partial_{xx}\phi(x)\,dt + o(dt).
$$

So even though the transition itself is unknown, its action on any smooth $\phi$ is computable from the two moments.

**Average over the starting point.** The expression above is conditional on the start $x$. To reach the unconditional density $p(y, t+dt)$, average it over where the process is now, against $p(x,t)$. This is the tower property: the unconditional expectation splits into the conditional inner expectation, just computed, averaged over the present density,

$$
\int \phi(y)\, p(y, t+dt)\, dy = \int \mathbb{E}\big[\phi(x + \Delta X) \mid X_t = x\big]\, p(x,t)\, dx = \int \Big(\phi(x) + \mu(x)\,\partial_x \phi\,dt + \tfrac{1}{2}\sigma^2(x)\,\partial_{xx}\phi\,dt\Big) p(x,t)\, dx + o(dt).
$$

Subtracting $\int \phi\, p(x,t)\,dx$ from both sides, dividing by $dt$ and letting $dt \to 0$ leaves

$$
\int \phi(y)\, \partial_t p(y,t)\, dy = \int \Big(\mu(x)\,\partial_x \phi(x) + \tfrac{1}{2}\sigma^2(x)\,\partial_{xx}\phi(x)\Big)\, p(x,t)\, dx.
$$

This is an integral identity, not yet the forward equation. The left integral is over the dummy variable $y$, which we rename $x$ so both sides share a variable and their integrands can be compared. To extract $\partial_t p$, we have to get rid of the $\phi$ derivatives on the right side, which is done by integration by parts: once for the drift term, twice for the diffusion term. Each pass requires $\phi$ to vanish at infinity so the boundary term drops, which is fine because $\phi$ is an arbitrary function. After $\phi$ is isolated on both sides, we get

$$
\int \phi\, \partial_t p \, dx = \int \phi \left(-\partial_x\big(\mu p\big) + \tfrac{1}{2}\partial_{xx}\big(\sigma^2 p\big)\right) dx.
$$

Two functions with the same integral against every such $\phi$ are equal, so the integrands agree, giving the forward equation

$$
\partial_t p = -\partial_x\big(\mu(x) p\big) + \tfrac{1}{2}\partial_{xx}\big(\sigma^2(x) p\big).
$$

This is the same Fokker-Planck equation the main text obtained in fewer steps.