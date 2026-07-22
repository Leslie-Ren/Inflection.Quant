---
title: "Local Volatility: From the Implied Vol Surface to Risk-Neutral Dynamics"
date: 2026-07-23
draft: true
math: true
tags: []
---

# Local Volatility: From the Implied Vol Surface to Risk-Neutral Dynamics

In the previous article we constructed the implied vol surface, which was primarily used to price vanilla options. But a great deal of what trades is not vanilla. Products like barriers and autocallables depend on the path the underlying takes, not only where it lands.

Suppose I price a barrier by Monte Carlo. At each step the spot sits at some level, and I need a volatility to advance it. What vol do I use? The surface gives me a vol for every strike, but simulation does not ask about strikes. It asks what volatility the spot experiences at this level, at this moment, which the surface cannot answer.

Local volatility is one solution. It gives the instantaneous vol of the underlying at each level and time. This article is my attempt to explain the model and its calibration, and to see how price and risk change when we price a barrier under local volatility instead of a constant one.

## The Local Volatility Model

Start with the model we already have. Black-Scholes assumes the underlying follows

$$dS_t = r S_t \, dt + \sigma \, S_t \, dW_t$$

with a single constant $\sigma$. That constant is why the model cannot fit the smile. Local volatility keeps the diffusion and the single Brownian driver, and lets the volatility vary with level and time:

$$dS_t = r S_t \, dt + \sigma(S_t, t) \, S_t \, dW_t$$

The function $\sigma(S_t, t)$ is the local volatility, the volatility the underlying diffuses with when the spot is at level $S$ at time $t$. The volatility itself carries no randomness.

## Calibrating Local Volatility

Once we know $\sigma(S, t)$, the model is fully specified: we can simulate paths by Monte Carlo or solve the pricing PDE. The hard part is getting $\sigma(S, t)$ in the first place.

My natural instinct is to reuse the approach from the implied vol article. There, we guessed a volatility, priced the option, compared to the market, and iterated until the price matched. Could we do the same here? Guess $\sigma(S, t)$, price the whole vanilla surface, compare, and adjust.

The trouble is that $S$ and $t$ are continuous, so every point $(S, t)$ the spot might visit needs its own volatility, and each vanilla price depends on the local vol at every point its paths pass through. The unknowns cannot be solved one at a time. Even discretized onto a grid, matching the surface by trial and error means solving for all the values together and repricing everything at each step. We need something better than guess-and-match.

Here is the shift in view that provides it. Instead of asking what prices a given $\sigma(S, t)$ produces, ask the reverse: given the prices, what $\sigma(S, t)$ must have produced them? If the quoted surface already contains enough information to determine the local volatility, we can read it off directly, with no iteration at all.

To see why that is even possible, borrow a simpler problem. Consider the heat equation,

$$u_t = a \, u_{xx},$$

which describes how a temperature profile $u(x, t)$ diffuses over time with diffusivity $a$. Normally we are given $a$ and solve for $u$. But suppose instead we could observe $u$ everywhere: the full profile at every position and time. Then $u_t$ and $u_{xx}$ are both things we can measure, and the equation rearranges to

$$a = \frac{u_t}{u_{xx}}.$$

The diffusivity is no longer something to guess. It is pinned by the solution itself. Observe enough of $u$ and the coefficient that generated it is determined.

Dupire's insight is that option prices offer exactly this situation. The quoted surface is not a single price; it is a price for every strike and every maturity, which is the analogue of observing $u$ everywhere. And the call price, viewed as a function of strike $K$ and maturity $T$, satisfies a forward equation in which the local variance plays the role of the unknown coefficient. Rearranging that equation for the coefficient gives Dupire's formula:

$$\sigma^2(K, T) = \frac{\dfrac{\partial C}{\partial T}}{\tfrac{1}{2} K^2 \dfrac{\partial^2 C}{\partial K^2}}.$$

The structure is the same as the heat equation read-off. The numerator is a time derivative of the observed surface. The denominator is a second derivative in the space variable, here the strike. Their ratio is the coefficient that generated the surface, evaluated at the level $S = K$ and time $t = T$. No guessing, no iteration: the local volatility is a direct calculation on the quotes.

Two pieces of that claim are asserted rather than shown: that the call surface satisfies a forward equation at all, and that its second strike derivative is the risk-neutral density. Both are derived in the appendix; here it is enough that the surface determines the coefficient, exactly as in the heat equation.