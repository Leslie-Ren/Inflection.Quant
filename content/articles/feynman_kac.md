---
title: "How Randomness Solves a Deterministic Equation: An Intuitive Look at the Feynman-Kac Theorem"
date: 2026-04-28
draft: false
math: true
tags: ["Feynman-Kac", "PDEs", "Stochastic Processes", "Risk-Neutral Pricing", "Monte Carlo"]
---

## Why This Matters

The first time I encountered the Feynman-Kac theorem, I found it fascinating but
unintuitive. The theorem claims that a deterministic PDE and the expectation of
a stochastic process are two representations of the same object. A PDE is smooth and
deterministic. A stochastic expectation involves randomness, probability measures, and
averaging over infinitely many paths. How could these be the same thing? I understood the steps of the proof, but I still didn’t have a clear intuition for why this equivalence should exist.

I also found myself slightly confused about its role in practice. In derivative pricing, we often work directly with risk-neutral expectations. The PDE formulation and the stochastic formulation both appear natural, so it is not immediately obvious what additional insight Feynman-Kac is adding.

This article is my attempt to answer both. Starting from a simple random walk, I hope
the equivalence feels less like a coincidence by the end, and that it becomes clear
where Feynman-Kac sits in derivative valuation and why it matters even when it may
appear unnecessary.

---

## The Intuition: A Random Walk

Before stating any theorem, I want to show through the simplest possible example why
a deterministic equation and a stochastic expectation are naturally the same thing.
The key is to start from neither, and instead start from something more primitive.

### Setup

Consider a particle that sits on a grid of points with spacing $\Delta x$ on a line.
Starting from position $x$ at time $t$, at each discrete time step of size $\Delta t$,
the particle moves up by $\Delta x$ or down by $\Delta x$ with equal probability
$\frac{1}{2}$. At the final time $T$, we collect a payoff $g(X_T)$ depending on where
the particle ends up.

We want to find a function $u(x, t)$ that tells us the fair value of this payoff at
any position $x$ and time $t$ before expiry.

### The Averaging Property

We do not yet say what $u$ is: not an expectation, not the solution to a PDE. We
only impose one requirement: $u$ must be consistent with the random walk. That is,
the value at $(x, t)$ must equal the average of the values at the two positions the
particle could reach at the next step:

$$u(x, t) = \frac{1}{2}u(x + \Delta x, t + \Delta t) +
\frac{1}{2}u(x - \Delta x, t + \Delta t)$$

with the boundary condition $u(x, T) = g(x)$.

This is the only thing we are asking of $u$. If we know the value at every position
at time $t + \Delta t$, the value at time $t$ must be the average of the two possible
next positions.

---

## Two Consequences of the Same Property

This single averaging requirement has two very different looking consequences, and
this is the heart of the intuition.

### Consequence 1: $u$ satisfies a PDE.

Rearranging the averaging equation and subtracting $u(x, t + \Delta t)$ from both
sides:

$$u(x, t) - u(x, t + \Delta t) =
\frac{1}{2}\left[u(x + \Delta x, t + \Delta t) - 2u(x, t + \Delta t) +
u(x - \Delta x, t + \Delta t)\right]$$

The left side is a difference in time. The right side is a second difference in
space. Dividing through by $\Delta t$, using the diffusion scaling
$\frac{(\Delta x)^2}{\Delta t} = 1$ (equivalently $\Delta x = \sqrt{\Delta t}$, the
same scaling condition established in
[Brownian Motion: From Random Walks to Option Prices]({{< relref "understanding_brownian_motion.md" >}})),
and taking $\Delta t \to 0$, $\Delta x \to 0$, we obtain:

$$\frac{\partial u}{\partial t} + \frac{1}{2}\frac{\partial^2 u}{\partial x^2} = 0
\quad \text{with } u(x, T) = g(x)$$

The key point is not the algebra itself, but the structure: the local averaging rule
forces a second-order spatial structure in the limit. That structure is the PDE.

### Consequence 2: $u$ equals a stochastic expectation.

The averaging property also tells us how to compute $u$ by working forward in time.
Starting from position $x$ at time $t$, at each step the particle moves up or down
with equal probability. After two steps there are four possible paths, after three
steps there are eight, and so on. This generates a binary tree of possible paths,
where each branch represents one possible realization of the particle's journey from
$t$ to $T$.

Each path through the tree has a probability: since every step is equally likely,
a path consisting of $k$ up-moves and $n-k$ down-moves over $n$ total steps has
probability $\left(\frac{1}{2}\right)^n$. The value $u(x, t)$ is the average of
$g(X_T)$ weighted by these path probabilities, which is exactly the expectation of
$g(X_T)$ over all paths:

$$u(x, t) = \mathbb{E}\left[g(X_T) \mid X_t = x\right]$$

In the continuous limit, as $\Delta t \to 0$ and $\Delta x \to 0$, the binary tree
of discrete paths converges to Brownian motion, and the sum over tree paths becomes
an expectation over continuous paths. The tree has not disappeared; it has become the
probability measure over continuous paths that defines the expectation.

### Summary

We started from one primitive
requirement: consistency with a local averaging rule. That single condition leads to
a deterministic PDE when viewed infinitesimally, and a stochastic expectation when
viewed globally. The PDE and the expectation are not two different models that happen
to agree. They are two ways of reading the same underlying structure: one in the
language of calculus, one in the language of probability. This is the intuition
behind Feynman-Kac.

---

## Feynman-Kac: The Theorem

Let $X_t$ be a stochastic process under measure $\mathbb{Q}$ defined by:

$$dX = \mu(X, t) dt + \sigma(X, t) dW^{\mathbb{Q}}$$

Consider the function $u(x, t)$ defined as the stochastic expectation:

$$u(x, t) = \mathbb{E}^{\mathbb{Q}}\left[e^{-\int_t^T r(X_s, s)\,ds}
g(X_T) \middle| X_t = x\right]$$

and the PDE:

$$\frac{\partial u}{\partial t} + \mu(x, t)\frac{\partial u}{\partial x} +
\frac{1}{2}\sigma^2(x, t)\frac{\partial^2 u}{\partial x^2} - r(x, t)u = 0$$

with terminal condition $u(x, T) = g(x)$.

**Feynman-Kac states that these two representations are equivalent.[^1]**

[^1]: The theorem holds under standard regularity conditions on $\mu$, $\sigma$, $r$,
and $g$; we assume these are satisfied throughout.

- If $u(x,t)$ is defined by the expectation above, then it satisfies the PDE.
- If $u(x,t)$ is defined as the solution to the PDE above, then it has the
  stochastic representation given by the expectation.

| PDE component | Stochastic counterpart |
|---|---|
| Drift $\mu \frac{\partial u}{\partial x}$ | Drift of $X_t$ |
| Diffusion $\frac{1}{2}\sigma^2 \frac{\partial^2 u}{\partial x^2}$ | Diffusion of $X_t$ |
| Discounting $-r u$ | Discount factor $e^{-\int_t^T r\,ds}$ inside the expectation |
| Terminal condition $u(x, T) = g(x)$ | Payoff function $g(X_T)$ |

---

## Two Paths to Derivative Valuation

When I first learned derivatives pricing, the PDE approach and the martingale approach
were presented as two separate tools to reach for depending on the problem. It took me
a while to appreciate that they are not just compatible but provably equivalent, and
that Feynman-Kac is precisely what makes that equivalence rigorous. To see why, it
helps to understand what each approach delivers on its own.

### Path 1: The PDE Approach

The PDE approach starts from no-arbitrage. We construct a delta-hedged portfolio,
eliminate the stochastic term, and impose that any risk-free portfolio must earn the
risk-free rate. In the Black-Scholes setting for a futures option, this gives the
Black (1976) PDE:

$$\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 F^2
\frac{\partial^2 V}{\partial F^2} - rV = 0, \quad V(F, T) = g(F)$$

Compared with the Black-Scholes PDE for a spot asset, there is no first-order term
$rF\frac{\partial V}{\partial F}$, because a futures position costs nothing to enter
and the hedge therefore carries no financing cost.

The PDE is grounded in no-arbitrage from the start. Any function that solves it is,
by construction, consistent with the requirement that a delta-hedged portfolio must
earn exactly the risk-free rate. Solve it once on a grid in $(F, t)$ space and we
obtain prices across all underlying levels and all times before expiry in a single
pass.

### Path 2: The Martingale Approach

The martingale approach starts from a different principle. Under the risk-neutral
measure $\mathbb{Q}$, the no-arbitrage condition is equivalent to discounted prices of
traded assets being martingales. A futures contract costs nothing to enter, so for
futures it is the futures price $F_t$ itself, not its discounted value, that is a
martingale under $\mathbb{Q}$. From this, any derivative can be priced as the expected
discounted payoff:

$$V(F, t) = \mathbb{E}^{\mathbb{Q}}\left[e^{-r(T-t)}g(F_T) \mid \mathcal{F}_t\right]$$

This is a clean and flexible framework. Prices can be computed by Monte Carlo, by
numerical integration, or analytically in some cases. The price it defines is
arbitrage-free by construction. But there is something this formula does not
immediately provide: the hedge.

Defining $V$ as a conditional expectation tells us the price. It does not show that
this price is the cost of a replicating strategy, or that it is the same function the
delta-hedging argument produces. For that, we need to know that this $V$ satisfies
the PDE from Path 1. If it did not, the two approaches would give different prices
for the same derivative.

### Where Feynman-Kac Comes In

Applying Feynman-Kac to the martingale pricing formula, where $F_t$ under $\mathbb{Q}$
has zero drift and diffusion $\sigma F$, tells us that $V$ defined by the expectation
satisfies:

$$\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 F^2
\frac{\partial^2 V}{\partial F^2} - rV = 0, \quad V(F, T) = g(F)$$

This is exactly the PDE from Path 1. The two approaches are not just compatible in the cases we can solve by hand. They are guaranteed to produce the same function for any well-posed diffusion model, whether or not an analytical solution exists. In simple models like Black-Scholes this equivalence can feel almost unnecessary, but in more complex models such as stochastic volatility settings, where closed-form solutions are no longer available, the multidimensional version of Feynman-Kac provides the rigorous link that ensures the PDE formulation and the expectation formulation remain consistent representations of the same quantity.

One caveat in that setting: volatility is not a traded asset, so the pricing PDE contains a market price of volatility risk, and each choice of that quantity corresponds to a choice of $\mathbb{Q}$. Feynman-Kac links the PDE and the expectation for a given choice; it does not make the choice.

---

## When Does It Matter Which Representation We Use?

Both representations are mathematically equivalent but not equally convenient for
every problem. In practice, choosing between a PDE solver and Monte Carlo is one of
the more common decisions in quantitative work, and the right answer depends on the
structure of the problem.

| Situation | Preferred approach | Reason |
|---|---|---|
| Computing smooth Greeks | PDE | Finite differences on the grid are stable; Monte Carlo differentiation is noisy |
| Model calibration | Forward PDE or semi-analytic pricer | Each calibration iteration requires fast, deterministic prices across many strikes and expiries; a forward (Dupire or Fokker-Planck) PDE delivers all strikes in a single solve, while Monte Carlo is slower and introduces noise into the objective function |
| Pricing across a range of underlying scenarios | PDE | A single grid solve covers all $F$ at once; Monte Carlo requires a separate simulation per scenario |
| High-dimensional underlyings (basket options) | Monte Carlo | PDE grid grows exponentially in dimension; simulation cost scales with paths |
| Path-dependent payoffs (Asian, lookback) | Monte Carlo | Path history requires extra state variables, turning a 2D grid into 3D or higher; simulation handles it naturally by following the full path |
| Validating a PDE implementation | Monte Carlo | Feynman-Kac guarantees a simulation-based check that should agree with the grid |

## Beyond Monte Carlo and PDE

Both the PDE and Monte Carlo perspectives assume a fixed underlying random evolution and differ only in how that evolution is computed. The PDE approach propagates structure deterministically, while Monte Carlo propagates it through sampled paths. Feynman-Kac tells us these are not competing methods but two representations of the same object.

This naturally leads to another question: we are computing expectations over paths, so why should we be committed to a single way of assigning probabilities to those paths in the first place? In many problems, the same physical or financial system can be described with different probabilistic weightings of the same trajectories, and some of these representations make computation or analysis significantly simpler than others. Understanding how such reweightings can change the apparent dynamics without changing the value of expectations is the subject of Girsanov’s theorem, which I will discuss in the [next article]({{< relref "girsanov.md" >}}).
