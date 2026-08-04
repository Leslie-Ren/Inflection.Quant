---
title: "Local Volatility: From the Implied Vol Surface to Risk-Neutral Dynamics"
date: 2026-07-23
draft: true
math: true
tags: []
---

# Local Volatility: From the Implied Vol Surface to Risk-Neutral Dynamics

## Why This Matters

In the [earlier article]({{< ref "vol_surface_calibration.md" >}}) we constructed the implied volatility surface and used it primarily to price vanilla options. But a great deal of what trades is not vanilla. Products like barriers and autocallables depend on the path the underlying takes, not only where it lands.

Suppose I price a barrier by Monte Carlo. At each step the spot sits at some level, and I need a volatility to advance it. What vol do I use? The surface gives me a vol for every strike, but simulation does not ask about strikes. It asks what volatility the spot experiences at this level, at this moment, which the surface cannot answer.

Local volatility is one solution. It gives the instantaneous vol of the underlying at each level and time. This article explains how the model works, and is my exploration of four questions:

- How should the local vols be calibrated in theory?
- What are the challenges of calibrating it in practice?
- If we price a barrier option under local volatility, how do the results differ from a constant vol, and why?
- Why can local volatility not be applied directly to commodities, and what adaptations does it need?

## Set-up

Start from what we already have. Black-Scholes gives the underlying one volatility for all levels and all times,

$$dS_t = r S_t \, dt + \sigma \, S_t \, dW_t$$

and that single $\sigma$ is exactly what leaves it unable to fit the smile. The smallest change we can make is to stop treating it as a constant:

$$dS_t = r S_t \, dt + \sigma(S_t, t) \, S_t \, dW_t$$

It remains a diffusion with a single Brownian driver, and $\sigma(S_t, t)$ carries no randomness of its own. It is the vol the spot diffuses with at level $S$ and time $t$.

## Local Vol Calibration and the Dupire Formula

Once we know $\sigma(S, t)$, the model is fully specified: we can simulate paths by Monte Carlo or solve the pricing PDE. The hard part is getting $\sigma(S, t)$ in the first place.

My default instinct is to reuse the approach from the implied vol article. There, we guessed a volatility, priced the option, compared to the market, and iterated until the price matched. Could we do the same here? Guess $\sigma(S, t)$, price the whole vanilla surface, compare, and adjust.

The trouble is that $S$ and $t$ are continuous, so every point $(S, t)$ the spot might visit needs its own volatility, and each vanilla price depends on the local vol at every point its paths pass through. The unknowns cannot be solved one at a time. Even discretized onto a grid, matching the surface by trial and error means solving for all the values together and repricing everything at each step. We need something better than guess-and-match.

### What the Market Already Gives Us

The market does give us the implied density directly. Based on Breeden-Litzenberger, derived [in the earlier article]({{< ref "vol_surface_calibration.md#breeden-litzenberger" >}}),

$$\frac{\partial^2 C}{\partial K^2}(K, T) = e^{-rT} p(K, T)$$

so for every expiry we quote, we know the market implied distribution of the underlying.

That is not yet enough. The density at time $t$ is the result of the spot diffusing from today through every level it passed on the way, so a single density reflects the accumulated effect of local vol everywhere before $t$, not the local vol precisely at $t$. Similarly, if you know an object's position at a single instant, you cannot infer its velocity, which is the instantaneous rate of change, and we need to take the derivative of position to get the velocity. Here if we know how the density changes with time, $\partial p / \partial t$, we can back out the volatility driving it.

Luckily, we already know how $\partial p / \partial t$ behaves. In the [forward and backward Kolmogorov article]({{< ref "forward_backward_pde.md" >}}) we covered how the density evolves with respect to time and derived the Fokker-Planck equation. For our diffusion it reads

$$\frac{\partial p}{\partial t} = -\frac{\partial}{\partial S}\left[ r S \, p \right] + \frac{1}{2} \frac{\partial^2}{\partial S^2}\left[ \sigma^2(S, t) S^2 p \right]$$

### Linking Densities to Option Quotes

Note the Fokker-Planck equation is in terms of $S$, and we have to link it to the information available from the market, which is option prices. Since the call price can be written as an integral over $S$,

$$C(K, T) = e^{-rT} \int_K^\infty (S - K) \, p(S, T) \, dS$$

multiplying Fokker-Planck by the payoff $(S - K)$ and integrating over the same range moves the whole relation into quote space. Evaluating at $t = T$, the maturity of the option we are pricing,

$$\int_K^\infty (S - K) \frac{\partial p}{\partial T} \, dS = \int_K^\infty (S - K) \left( -\frac{\partial}{\partial S}\left[ r S \, p \right] + \frac{1}{2} \frac{\partial^2}{\partial S^2}\left[ \sigma^2(S, T) S^2 p \right] \right) dS$$

The left side is straightforward. Differentiate the call price with respect to $T$, then solve for the integral:

$$\int_K^\infty (S - K) \frac{\partial p}{\partial T} \, dS = e^{rT} \left( \frac{\partial C}{\partial T} + r C \right)$$

The right side involves partial derivatives in $S$, so we evaluate both terms by integration by parts. Boundary terms vanish at $S = K$ because the payoff is zero there, and as $S \to \infty$ provided $p$ decays fast enough to suppress the polynomial growth of $S^2 p$, which holds for any density with a finite second moment.

The drift term takes one pass, giving

$$\int_K^\infty (S - K) \left( -\frac{\partial}{\partial S}\left[ r S p \right] \right) dS = r e^{rT} C - r K e^{rT} \frac{\partial C}{\partial K}$$

The diffusion term takes two passes, and the second collapses the integral onto the single point $S = K$, which is what isolates the local vol at one strike out of a relation that held across all levels:

$$\int_K^\infty (S - K) \cdot \frac{1}{2} \frac{\partial^2}{\partial S^2}\left[ \sigma^2 S^2 p \right] dS = \frac{1}{2} \sigma^2(K, T) K^2 e^{rT} \frac{\partial^2 C}{\partial K^2}$$

Equating the two sides, the $rC$ terms cancel:

$$\frac{\partial C}{\partial T} = -r K \frac{\partial C}{\partial K} + \frac{1}{2} \sigma^2(K, T) K^2 \frac{\partial^2 C}{\partial K^2}$$

Rearranged, this is the **Dupire formula**:

$$\sigma^2(K, T) = \frac{\dfrac{\partial C}{\partial T} + r K \dfrac{\partial C}{\partial K}}{\dfrac{1}{2} K^2 \dfrac{\partial^2 C}{\partial K^2}}$$

Beyond the derivation we just did above, why does the local variance take on such a form? I find it more intuitive to see it in discrete form:

$$e^{rT} \Delta C \approx r K \, \Delta T \cdot \mathbb{Q}(S_T > K) + \frac{1}{2} \sigma^2(K, T) K^2 \, \Delta T \cdot p(K, T)$$

Dupire states that the additional value of the call option from extending its maturity slightly, from $T$ to $T + \Delta T$, is driven by two factors: how the process drifts, assuming the underlying is in the money, and the diffusion, how violently the spot moves, weighted by the density at the strike.

Dupire gives us a powerful tool to compute the local vol at every $(K, T)$ from call prices. No iteration, no joint solve needed.

## How Local Vol is Calibrated in Practice

To apply Dupire directly, we need to know how the call price changes for a small change in $K$ and $T$. That is not available from the market, since quotes are not given at continuous $K$ or $T$. But we can calculate what the call price should be if an implied vol surface is in place, which translates any $K$ and $T$ into a call price. That is exactly what the earlier implied vol article was for.

A surface is usually displayed against strike or delta and maturity, so it might be tempting to differentiate against $K$ and $T$. But they are not the natural coordinates to view the vols in. What drives the dispersion of the price is the variance, and that variance is measured on the return rather than on the absolute price level. So the quantities the problem actually depends on are total variance and log-moneyness, and that is what we will use to derive the local vol form from the vol surface.

With $x = \log(F_T/K)$ the log-moneyness and $w(x, T)$ the total implied variance, substitute the Black-Scholes formula for $C$ and apply the chain rule to convert the strike and maturity derivatives into derivatives of $x$ and $w$. We get

$$\sigma^2(x, T) = \frac{\dfrac{\partial w}{\partial T}}{1 - \dfrac{x}{w}\dfrac{\partial w}{\partial x} + \dfrac{1}{4}\left( -\dfrac{1}{4} - \dfrac{1}{w} + \dfrac{x^2}{w^2} \right)\left( \dfrac{\partial w}{\partial x} \right)^2 + \dfrac{1}{2}\dfrac{\partial^2 w}{\partial x^2}}$$

This is **Dupire in Gatheral's form**, and it provides a bridge to compute the local vol from the implied vol surface. The numerator is how the total variance changes with respect to time. The denominator appears convoluted but has a simple meaning: it is the ratio of the modeled density from the implied vol surface against the baseline model, the constant vol lognormal, which is the function $g(x)$ from [Appendix B]({{< ref "vol_surface_calibration.md#appendix-b-deriving-the-total-variance-butterfly-condition" >}}) of the implied vol article. When the vol is constant, $\partial w / \partial x$ and $\partial^2 w / \partial x^2$ are both zero, the denominator is $1$, and we are left with $\sigma^2 = \partial w / \partial T$, which is the constant vol itself.

### Taking the Derivatives

The Gatheral formula needs three derivatives of $w$, and the naive way to get them is to bump $x$ and $T$ and take finite differences. It works, but it is not the default way to do it. Differencing amplifies the noise in the implied vol fit, especially for the second derivative in the denominator. The result is a local vol grid that looks jagged even when the input surface is smooth, and the pricing model, Monte Carlo or PDE, inherits that noise.

For $\partial w / \partial x$ and $\partial^2 w / \partial x^2$, analytical derivatives are usually preferred, since they give a closed form on a smooth calibration. The cubic spline from the implied vol article is piecewise polynomial, so its first and second derivatives are available in closed form on each segment and are continuous across the knots by construction. Parametric models such as SSVI are also differentiable in closed form, and parametric surfaces are often the input used in practice.

The maturity derivative, $\partial w / \partial T$, works similarly. The analytical form follows from the total variance interpolation set by the implied vol model. One nuance is worth noting. If total variance is linear in $T$ at fixed $x$, the derivative is piecewise constant and undefined at the quoted expiries themselves, since the intervals on either side give different slopes. On each interval it is well defined, and since variance accumulates as an integral over time, the value assigned at an isolated point does not affect any price. The usual convention is to take the slope of the interval ahead, which leaves the local vol surface right-continuous in maturity and jumping as the simulation crosses each strip.

The analytical approach tells us that in order to use the implied vol surface as an input, instead of just saving the vol at all the grid points, we should store the model itself, spline coefficients or parameters, and take the derivatives from it directly. Local vol calibration still runs downstream of the implied vol fit, but what crosses between them is the model rather than a set of grid points.

In some settings, passing the model along is not possible. Traders may mark the surface directly, sometimes with manual adjustments, or the marking system and the calibration library may be separate enough that only a grid of vols crosses between them. In that case a surface has to be fitted to the grid before anything can be differentiated. The risk here is that the refitted model may differ from the one underlying the original marked vol surface.

### What the Input Surface Must Satisfy

Local variance is a ratio, so it exists and is positive only when the numerator and denominator are. Both conditions are the arbitrage checks from the earlier article, arriving here as requirements rather than as diagnostics.

The numerator is positive when total variance is non-decreasing in maturity, which is precisely the absence of calendar arbitrage. A calendar violation produces a negative local variance and the calibration fails outright at that point.

The denominator is positive when the smile satisfies the butterfly condition, which is the Gatheral $g$-function test. A butterfly violation is a negative density, and dividing by it produces a negative local variance for the same reason. This is what the planted violation in the Dec 2026 strip does if it is left in.

Two further conditions are less visible. The surface must be twice differentiable in $x$, so a fit that is merely continuous will not do. And the wings must not grow so fast that the density loses the moments the derivation assumed, which is the condition the linear total variance extrapolation was chosen to respect.

There is also a softer failure. Far into the wings the density is nearly zero, so the denominator is nearly zero, and local variance becomes the ratio of two small numbers. It may stay positive while being numerically meaningless. The market has almost no information about local vol where the spot rarely goes, and no amount of care in the fit creates it.

### Precomputing the Grid

Evaluating Gatheral's formula inside a Monte Carlo loop would mean recomputing derivatives at every step of every path. In practice the local vol is computed once on a grid in $(x, T)$, stored, and interpolated during pricing. The grid is chosen to cover the region the paths will actually visit, with enough resolution near the money and in maturity to keep the interpolation error below the discretization error of the scheme.

