---
title: "Finite Difference Methods: Marching Forward or Solving Together"
date: 2026-06-09
draft: false
math: true
tags: ["finite-difference", "pde", "black-scholes", "numerical-methods", "heat-equation"]
---

## Why This Matters

A derivative price can be computed two equivalent ways: as a risk-neutral expectation, or as the solution of a PDE. This is the [Feynman-Kac result]({{< relref "feynman_kac.md" >}}), which I explored in the earlier article. Monte Carlo is the natural way to handle the expectation, and the [previous article]({{< relref "mc_variance_reduction.md" >}}) worked through techniques for making it more efficient. Here I want to look at the other side, where the price is the solution of a PDE and we solve it on a grid.

What drew me to the PDE approach is how much falls out of that single grid solve. The price at any spot is just a point on the grid, so one solve covers a whole range of scenarios at once. The Greeks are numerical derivatives of the same grid, so they come almost for free. Early exercise, as in an American option, and knock-out barriers both slot in naturally, the first as a comparison at each grid point and the second as a boundary condition. Monte Carlo can do all of these too, but on the PDE side they feel less like add-ons and more like things the method was already set up to handle.

The workhorse for solving these PDEs numerically is the finite difference method, and that is what the rest of the article is about. Rather than start with the Black-Scholes PDE directly, I will build the method up on the heat equation from physics. The two share the same diffusion structure, but the heat equation has constant coefficients where Black-Scholes has variable ones, so it shows the numerical ideas without the algebraic clutter. Everything we develop carries over to Black-Scholes with those extra terms added back in.

---

## The Heat Equation

Picture a thin metal rod. Heat a spot on it, then leave it alone. The warm spot cools and spreads, and the rod drifts toward a uniform temperature. We want an equation for how the temperature $u(x, t)$ at position $x$ and time $t$ evolves.

Rather than reach for a known equation, let us ask what properties it must have. Heat conduction is local: a point exchanges heat only with its immediate neighbours, not with the rod as a whole. So the rate of change $\partial_t u$ at a point should be driven by how that point compares to its neighbours on each side.

Consider a small segment around a point. If the left neighbour is one degree cooler and the right one degree warmer, the point loses heat to the left at the same rate it gains from the right. The flows balance and the temperature holds steady, despite the clear slope through the point. So the slope $\partial_x u$ is not what drives the change; a constant slope produces balanced flow and no change at all. What drives the change is an imbalance between the two sides, the slope being steeper on one side than the other. That imbalance is the curvature $\partial_{xx} u$, and it gives the heat equation:

$$
\partial_t u = \kappa\, \partial_{xx} u, \qquad \kappa > 0,
$$

with $\kappa$ a positive constant, set by the material, controlling how fast heat spreads. The sign of the curvature sets the direction: a point that rises above its neighbours has negative curvature in $x$ and cools, a dip below them has positive curvature and warms, a straight run holds still.

The equation governs how each interior point evolves, but it says nothing about the two ends of the rod, which have a neighbour on only one side. On its own it does not determine a unique solution; we also have to say what happens at the ends. The simplest choice is to hold each end at a fixed temperature: imagine the ends sitting in an ice bath, pinned at zero while the interior evolves. This is a boundary condition, and the equation needs one at each end before the problem is complete. It is why every curve in the figure below is tied down to zero at the two edges.

{{< heat_rod_evolution >}}

The heat starts concentrated at a single point in the middle of the rod. The figure shows the profile shortly after, spreading outward into a bump that widens and flattens over time, conserving total heat while smoothing it out. This diffusion is the behaviour we will be approximating numerically.

---

## Discretising the Equation

To solve the heat equation numerically we replace the continuous derivatives with finite differences on a grid. Space is divided into points $x_0, x_1, \ldots, x_M$ spaced $\Delta x$ apart, and time into levels $t_0, t_1, \ldots, t_N$ spaced $\Delta t$ apart. We write $u_i^n$ for the approximation to $u(x_i, t_n)$.

For the time derivative, the natural finite difference uses the values one step apart in time:

$$
\partial_t u \approx \frac{u_i^{n+1} - u_i^n}{\Delta t}.
$$

This is the left side of the equation.

The right side is the spatial curvature $\kappa\, \partial_{xx} u$, with the standard central-difference approximation:

$$
\partial_{xx} u(x_i) \approx \frac{u_{i-1} - 2 u_i + u_{i+1}}{(\Delta x)^2}.
$$

This is the discrete form of the curvature: it measures how far $u_i$ sits above or below the average of its two neighbours, the same quantity that drove the heat equation.

The equation sets the time difference equal to this curvature, but it leaves one thing open: at which time level do we measure the curvature? At $t_n$, where we already know all the values, or at $t_{n+1}$, where we do not? Nothing so far decides it. We have to choose, and that single choice produces two schemes with very different character.

---

## The Explicit Scheme

The natural choice is to enforce the equation at $t_n$, the time level we already know. The spatial derivative is computed from known values, and the update reads:

$$
\frac{u_i^{n+1} - u_i^n}{\Delta t} = \kappa\, \frac{u_{i-1}^n - 2 u_i^n + u_{i+1}^n}{(\Delta x)^2}.
$$

Solving for the one unknown gives a direct formula:

$$
u_i^{n+1} = u_i^n + \lambda\,(u_{i-1}^n - 2 u_i^n + u_{i+1}^n), \qquad \lambda = \frac{\kappa\, \Delta t}{(\Delta x)^2}.
$$

Everything on the right is known, so we read off the new value at each grid point directly. Rewriting slightly shows what the update is doing:

$$
u_i^{n+1} = \lambda\, u_{i-1}^n + (1 - 2\lambda)\, u_i^n + \lambda\, u_{i+1}^n.
$$

The value at the next time step is a weighted combination of the current value at that point and its two neighbours. Each point is updated independently, which makes the scheme cheap and easy to implement.

{{< fdm_stencil_explicit >}}

The stencil shows what the update reaches for: the new value at $t_{n+1}$ is built from the three known values in the previous time column at $t_n$.

This update applies only to the interior points. The two end points have a neighbour on just one side, so the three-point stencil would run off the grid. We do not compute them from the update at all: the boundary condition fixes them directly, setting the edge values to zero at every step. The first interior point still reaches for its outer neighbour, but that neighbour is now the known boundary value, so its update is well defined.

---

## When the Explicit Scheme Breaks

The explicit scheme works well when the time step is small, but the simplicity hides a trap. Recall that the update writes the next value as a combination of three current values, with weights $\lambda$, $1 - 2\lambda$, $\lambda$ that always sum to one. As long as all three are non-negative, which holds when $\lambda \leq 1/2$, the next value is a weighted average of the three, so it lands somewhere between the coolest and warmest of them. A hot point cools toward its neighbours, never dropping below them in a single step, which is how heat should behave. Starting from a single point at temperature 1 between two neighbours at 0, with $\lambda = 1/3$, the heat simply spreads and fades, every value staying between 0 and 1.

{{< stability_heatmap_stable >}}

Once $\lambda > 1/2$, the central weight turns negative and the update is no longer an average. Now a point can be pushed past its neighbours instead of settling between them. Let us look at the same starting point with $\lambda = 1$, and watch what the update does step after step.

{{< stability_heatmap >}}

Follow the centre row. The point starts at 1, and a single step sends it to -1, already colder than the neighbours it was meant to settle toward. The next step overshoots it back to 3, then to -7, then to 19. Each step flips the sign and grows the magnitude, in space and in time, producing the intensifying checkerboard above. Left to run, it swamps the real solution within a few steps.

The condition for the weights to stay non-negative is

$$
\lambda = \frac{\kappa\, \Delta t}{(\Delta x)^2} \leq \frac{1}{2},
$$

which caps the time step at $\Delta t \leq (\Delta x)^2 / (2\kappa)$. The cost of this is steep. The maximum step shrinks with the square of the spatial spacing, so using ten times as many spatial points forces the time step to shrink by a factor of a hundred. Ten times the spatial points and a hundred times the time steps is a thousandfold increase in total work, all to refine the spatial grid tenfold. To stay stable on a fine grid, the explicit scheme is forced into a very large number of small steps.

---

## The Implicit Scheme

Given the stability constraint of the explicit scheme, can we find an alternative that does not suffer from it? The explicit scheme derives the future value directly from the present: it reads the current curvature and steps forward, and when the step is too large, that forward step overshoots. Instead of deriving the future value, what if we constrain it, requiring the future value to satisfy the equation? We enforce the equation at $t_{n+1}$, evaluating the spatial derivative at the level we are solving for:

$$
\frac{u_i^{n+1} - u_i^n}{\Delta t} = \kappa\, \frac{u_{i-1}^{n+1} - 2 u_i^{n+1} + u_{i+1}^{n+1}}{(\Delta x)^2}.
$$

Now the right-hand side contains three unknowns, so we cannot solve for $u_i^{n+1}$ on its own. Collecting the unknowns on one side:

$$
-\lambda\, u_{i-1}^{n+1} + (1 + 2\lambda)\, u_i^{n+1} - \lambda\, u_{i+1}^{n+1} = u_i^n.
$$

{{< fdm_stencil_implicit >}}

The stencil shows why: the new value connects to its spatial neighbours in the same time column at $t_{n+1}$, which are themselves unknown, rather than only to known values at $t_n$.

There is one such equation at every interior grid point, giving a linear system in all the unknowns at $t_{n+1}$. In matrix form, with $\mathbf{u}^n$ the vector of grid values,

$$
A\, \mathbf{u}^{n+1} = \mathbf{u}^n,
$$

where $A$ is tridiagonal: $1 + 2\lambda$ on the diagonal, $-\lambda$ on the two off-diagonals. Each step requires solving this system rather than reading off an answer. We solve it directly with a standard tridiagonal algorithm rather than forming the inverse $A^{-1}$: the inverse is dense, and building it would discard the structure that makes the problem cheap. A general dense system of $N$ unknowns costs on the order of $N^3$ operations to solve, but the tridiagonal system, being almost all zeros, can be solved in work proportional to $N$, nearly as cheap as a single sweep over the grid.

The payoff is that the implicit scheme is stable for any time step. There is no $\lambda \leq 1/2$ constraint, so the step is no longer limited by stability. Each step costs more than an explicit step, but the freedom to take large steps wins comfortably on the fine grids where the explicit scheme would be forced into a crawl.

---

## Why the Implicit Scheme Cannot Blow Up

The explicit update cools a hot point by an amount that grows with the step size. When the step is large it removes more than the whole gap to the neighbours, so the point overshoots past them. Subtraction can remove too much.

The implicit update does not work by subtraction. Isolating the centre term in its equation,

$$
(1 + 2\lambda)\, u_i^{n+1} = u_i^n + \lambda\left(u_{i-1}^{n+1} + u_{i+1}^{n+1}\right),
$$

the new value carries a factor of $1 + 2\lambda$, so finding it means *dividing* by that factor. The denominator is always greater than one, no matter how large $\lambda$ grows, so the division can only shrink a value toward zero, never flip its sign or amplify it. A sharp feature is divided down harder as the step grows, which is what diffusion should do to it. No step size turns the division into an overshoot.

The same fact has a physical reading, and it comes down to which state the equation is enforced on. The implicit scheme enforces the equation on the new state, so the new values are required to satisfy it. This acts as a check: if the solver tried to overshoot a hot point into a dip below its neighbours, that dip would have positive curvature, and the equation says a point with positive curvature warms rather than sits cold. The overshoot would violate the equation the new state must obey, so it cannot be the solution. The explicit scheme enforces the equation only on the old, known state and computes the new value from it directly, with nothing constraining the new value itself. That is the difference: implicit holds the new state to the equation, explicit does not.

---

## Reflection: Marching Forward or Solving Together

Looking past the formulas, the two schemes hold different views of how the system evolves.

- **Explicit:** the known state generates the next one. For the heat equation, where we solve forward from the present, this is literally the present generating the future.
- **Implicit:** the step is a constraint. The new state is the one that satisfies the equation at every point at once, reached by solving rather than generating.

Both converge to the same continuous solution as the grid is refined, since the gap between $t_n$ and $t_{n+1}$ vanishes and the choice of where to enforce the equation stops mattering. What differs is the stance: marching forward, or solving a constraint. That stance is what made one scheme able to overshoot and the other not.

---

## Crank-Nicolson

So far we have only considered stability. The other concern is accuracy: how close the answer is to the true solution. A Taylor expansion sizes the error. For the time difference, expanding $u^{n+1}$ about $t_n$,

$$
\frac{u^{n+1} - u^n}{\Delta t} = \partial_t u + \frac{\Delta t}{2}\, \partial_{tt} u + \cdots,
$$

so the error is proportional to $\Delta t$; the same expansion on the spatial neighbours gives a central difference with error proportional to $(\Delta x)^2$. Both vanish as the grid is refined. The time error is only first-order, and once we take the large steps the implicit scheme allows, it dominates: driving it down forces us back into many small steps, the very thing going implicit was meant to avoid.

So how can we make the scheme converge faster in time? Let us borrow the idea behind computing a delta numerically. To estimate the slope of an option value in the underlying, a one-sided difference $\frac{f(x+h) - f(x)}{h}$ is first-order, but a centred difference $\frac{f(x+h) - f(x-h)}{2h}$ is second-order, which is why the centred form is the usual choice: its errors on the two sides cancel. The lesson is that a central estimate beats a one-sided one. Explicit and implicit are both one-sided in time, each evaluating the spatial term at one end of the step, which is what makes their time error first-order.

Crank-Nicolson applies the same fix in time: it evaluates the spatial term at both ends of the step and averages them, centring it at the midpoint.

$$
\frac{u_i^{n+1} - u_i^n}{\Delta t} = \frac{1}{2}\Big( \underbrace{\kappa\, \tfrac{\delta^2 u_i^n}{(\Delta x)^2}}_{\text{explicit}} + \underbrace{\kappa\, \tfrac{\delta^2 u_i^{n+1}}{(\Delta x)^2}}_{\text{implicit}} \Big), \qquad \delta^2 u_i = u_{i-1} - 2u_i + u_{i+1}.
$$

{{< fdm_stencil_cranknicolson >}}

The stencil reaches into both time levels: the three points at $t_n$ that the explicit scheme uses, and the coupled points at $t_{n+1}$ that the implicit scheme uses. Centring this way makes the time error second-order, so a given accuracy needs far fewer steps. The unknowns at $t_{n+1}$ are still coupled, so each step is a tridiagonal solve, and the scheme keeps the implicit scheme's unconditional stability.

---

## From the Heat Equation to Black-Scholes

The same three schemes carry over to the Black-Scholes PDE,

$$
\partial_t V + \tfrac{1}{2}\sigma^2 S^2\, \partial_{SS} V + r S\, \partial_S V - r V = 0.
$$

It has the same diffusion structure as the heat equation, but not in pure form. The drift term $rS\,\partial_S V$ and the state-dependent diffusion $\tfrac{1}{2}\sigma^2 S^2$ break the symmetry of the heat equation's clean diffusion: the coefficient $\tfrac{1}{2}\sigma^2 S_i^2$ now varies from point to point rather than being constant, and the first-derivative drift adds a directional push that the heat equation does not have. Neither changes how the schemes are built; both just enter the grid as additional terms.

As before, the choice of scheme is a choice of which level to evaluate the PDE terms at. Explicit reads the second derivative, the first derivative, and the discount term $rV$ all at the known level, so each unknown value is computed directly. Implicit takes all three at the unknown level, which couples the unknowns into a tridiagonal system. The coefficient $S_i$ is the same either way; only the level of the $V$ values changes. Explicit, implicit, and Crank-Nicolson all carry over, with only the matrix entries differing.

One thing does change when we move from the heat equation to Black-Scholes. Pricing is a terminal-value problem: we know the payoff at expiry and solve backward toward today, so the computation marches backward in time. This sharpens what the explicit scheme is really doing. It does not march forward in time; it marches from the known level to the unknown one. For the heat equation the known level is the present, so the two look the same. For Black-Scholes the known level is at expiry, so the march runs backward. Time direction is incidental; known to unknown is the point.

Two examples show what the PDE approach buys us, each handled naturally on the grid: the early-exercise decision of an American option, and the knock-out level of a barrier option.

---

## Early Exercise: American Options

An American option can be exercised at any time before expiry, which adds a decision at every point: hold the option, or exercise it now. The PDE handles this naturally because it already solves backward in time.

At each time level, the backward solve produces the continuation value, what the option is worth if held, computed from the same finite difference machinery as the European case. Write it $V_i^{\text{cont}}$. Against this we compare the immediate exercise value, the payoff from exercising now: for an American put with strike $K$, that is $\max(K - S_i, 0)$. The option value at the node is whichever is larger:

$$
V_i^n = \max\left(\,V_i^{\text{cont}},\ \max(K - S_i, 0)\,\right).
$$

So $V_i^{\text{cont}}$, the continuation value, is what the PDE step alone gives, the same value a European option would have. The American value $V_i^n$ is that with the early-exercise choice layered on top. The distinction matters for the next step: it is $V_i^n$, the post-maximum value, that is carried backward into the following solve, so the right to exercise early at a later node feeds into the value at earlier ones.

This comparison is applied pointwise at every grid point and time level. The grid already holds the value at every spot and time, so the question "is it worth exercising here and now" can be asked everywhere at once. The boundary between the hold region and the exercise region falls out of the solve, with no special treatment beyond the pointwise maximum.

---

## Path Dependence: Barrier Options

A down-and-out call is worthless if the underlying ever falls to the barrier $B$ before expiry. Take the barrier to be continuously monitored, live at all times. On the grid this is a boundary condition: the value is held at zero at every point at or below the barrier, at every time level.

$$
V(S, t) = 0 \quad \text{for } S \leq B.
$$

This is the ice bath again. Just as the rod's ends were held at zero and the chill spread inward, the value is held at zero along the barrier and the backward solve carries that zero into the living region above it. A spot just above the barrier has a neighbour pinned at zero, which drags its value down, and that influence spreads further up the grid at each step. By the time the solve reaches today, the value at a spot well above the barrier has already been reduced by the chance of drifting down and knocking out, without any path ever being traced.

This is where the grid picture and the path picture part ways most sharply. In Monte Carlo the barrier is an event in time: a trajectory crosses $B$ and the option dies. The PDE has no trajectories, and it cannot, since the equation relates each point only to its neighbours, leaving nowhere for a path's past to live. Instead the barrier becomes geometry: the edge of the region the equation is solved on. The question "did the path cross the barrier" is replaced by "where is the edge of the domain," and the backward solve prices the knockout by letting that zero-edge bleed inward.

This also avoids the monitoring bias that troubles Monte Carlo, where the barrier is checked only at discrete sample times and a path can dip below and back between checks, the same bias the [Monte Carlo article]({{< relref "mc_variance_reduction.md" >}}) had to correct for on this very option. The grid has no unobserved trajectory between steps; the barrier condition applies across the whole domain at once. The PDE still carries its own time-discretisation error, but not this bias. One practical point: the barrier should land on a grid line. If $B$ falls between two spatial points, the effective barrier is misplaced by up to a grid spacing, so the $S$ grid is usually aligned to put a node exactly on $B$.

This clean treatment is for a continuously-monitored barrier. A discretely-monitored barrier, checked only at set dates, is handled by imposing the zero condition only on those dates and letting the region below $B$ stay alive in between, which is closer to what Monte Carlo does and loses some of the tidiness.

---

## Conclusion

Explicit and implicit schemes solve the same equation in two ways: marching forward off the values we know, or solving for a level jointly consistent with them. The first is simple but carries a stability limit; the second costs a linear solve per step and removes the stability constraint. Crank-Nicolson sits between them and gains an order of accuracy. It all traces back to one choice, at which time level the equation is enforced, the same choice that reappears, flipped in direction, when we move from the forward-solved heat equation to the backward-solved Black-Scholes PDE. And what makes the approach worth the trouble is everything that comes after the price: Greeks fall out of the grid as finite differences, a single solve values the option across every spot, early exercise reduces to a pointwise maximum, and level-based barriers to a boundary condition. The grid holds the whole solution surface, and most of the questions worth asking can be read directly off it.

**One last thread, left for you to pull.** The explicit scheme is stable only when $\Delta t$ scales with $(\Delta x)^2$. If that relationship looks familiar, it is the same scaling that turns a random walk into Brownian motion, the one explored in the [Brownian motion article]({{< relref "brownian_bridge.md" >}}): step variance grows like $\Delta t$, step size like $\Delta x$, and the two must balance as $(\Delta x)^2 \sim \Delta t$ for the walk to converge to a diffusion. That is probably not a coincidence. By [Feynman-Kac]({{< relref "feynman_kac.md" >}}), the heat equation is the deterministic face of a Brownian motion, so a scheme for it is, in some sense, simulating that process. Write the stable update with its weights summing to one, and read them as the chances of stepping up, staying, and stepping down. Then ask what the stability condition is really demanding of those weights, and whether a scheme that violates it is still describing a random walk at all. I think the answer is worth finding yourself.