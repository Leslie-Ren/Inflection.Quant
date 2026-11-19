---
title: "Local Volatility: From the Implied Vol Surface to Risk-Neutral Dynamics"
date: 2026-07-23
draft: False
math: true
tags: ["local volatility", "Dupire", "Gatheral", "implied volatility", "Gyöngy", "barrier options", "marginal density", "commodities", "Monte Carlo"]
ShowToc: true
TocOpen: false
---

## Why This Matters

In an [earlier article]({{< ref "vol_surface_calibration.md" >}}) we constructed the implied volatility surface and used it primarily to price vanilla options. But a great deal of what trades is not vanilla. Products like barriers and autocallables depend on the path the underlying takes, not only where it lands.

Suppose I price a barrier by Monte Carlo. At each step the spot sits at some level, and I need a volatility to advance it. What vol do I use? The implied vol surface gives me a vol for every strike. A simulation needs something else: a vol for the level the simulated path has reached at that step.

The local vol model provides a solution. It assigns a vol to every spot level and time, chosen so that the model reprices every vanilla option on the implied vol surface. This article starts with how local vol is calibrated, both in theory and in practice. I then explore a few questions:

- Local vol and Black-Scholes give the same price for a vanilla option, but do they agree on its risk and P&L attribution?
- Given that the two models agree on vanilla prices, where does their difference in a barrier option price come from?
- Local vol is widely used for equities, but why can't it be applied directly to some commodities such as natural gas?

## From Constant Vol to Local Vol

Start from what we already have. Black-Scholes gives the underlying one volatility for all levels and all times,

$$dS_t = r S_t \, dt + \sigma \, S_t \, dW_t$$

so every step of the MC simulation would use the same vol, no matter where the path is. Black-Scholes can still match every vanilla by using a different $\sigma$ for each strike, which is what the implied vol surface does, but then each option is priced under its own model. A path dependent product needs one model for the whole path. The smallest change we can make is to let the vol depend on the level and time:

$$dS_t = r S_t \, dt + \sigma(S_t, t) \, S_t \, dW_t$$

It remains a diffusion with a single Brownian driver, and $\sigma(S_t, t)$ carries no randomness of its own. It is the vol the spot diffuses with at level $S$ and time $t$.

Throughout the article we assume the stock pays no dividends, to keep the focus on the local vol model itself. Real dividends are paid in cash on discrete dates, and their future amounts are not known in advance, so modelling them properly is a topic of its own. For anyone interested, Buehler (2007) shows how local vol can be set up with discrete dividends, each a mix of a fixed cash amount and a proportion of the stock price.

## Local Vol Calibration and the Dupire Formula

Once we know $\sigma(S, t)$, the model is fully specified: we can simulate paths by Monte Carlo or solve the pricing PDE. The hard part is getting $\sigma(S, t)$ in the first place.

My default instinct is to reuse the approach from the implied vol article. There, we guessed a volatility, priced the option, compared to the market, and iterated until the price matched. Could we do the same here? Guess $\sigma(S, t)$, price the whole vanilla surface, compare, and adjust.

The trouble is that $S$ and $t$ are continuous, so every point $(S, t)$ the spot might visit needs its own volatility, and each vanilla price depends on the local vol at every point its paths pass through. The unknowns cannot be solved one at a time. Even discretized onto a grid, matching the surface by trial and error means solving for all the values together and repricing the whole surface at each step. We need something better than guess-and-match.

### What the Market Already Gives Us

The market does give us the implied density directly. Based on Breeden-Litzenberger, derived [in the earlier article]({{< ref "vol_surface_calibration.md#breeden-litzenberger" >}}),

$$\frac{\partial^2 C}{\partial K^2}(K, T) = e^{-rT} p(K, T)$$

so for every expiry we quote, we know the market implied distribution of the underlying.

That is not yet enough. The density at time $t$ is the result of the spot diffusing from today through every level it passed on the way, so a single density reflects the accumulated effect of local vol everywhere before $t$, not the local vol precisely at $t$. Knowing where an object is at one instant does not tell you its velocity; for that you need how the position changes. Likewise, knowing how the density changes with time, $\partial p / \partial t$, lets us back out the volatility driving it.

Luckily, we already know how $\partial p / \partial t$ behaves. In the [forward and backward Kolmogorov article]({{< ref "forward_backward_pde.md" >}}) we covered how the density evolves with respect to time and derived the Fokker-Planck equation. For our diffusion it reads

$$\frac{\partial p}{\partial t} = -\frac{\partial}{\partial S}\left[ r S \, p \right] + \frac{1}{2} \frac{\partial^2}{\partial S^2}\left[ \sigma^2(S, t) S^2 p \right]$$

### Linking Densities to Option Quotes

Note the Fokker-Planck equation is in terms of $S$, and we have to link it to the information available from the market, which is option prices. Since the call price can be written as an integral over $S$,

$$C(K, T) = e^{-rT} \int_K^\infty (S - K) \, p(S, T) \, dS$$

multiplying Fokker-Planck by the payoff $(S - K)$ and integrating over the same range moves the whole relation into quote space. Evaluating at $t = T$, the maturity of the option we are pricing,

$$\int_K^\infty (S - K) \frac{\partial p}{\partial T} \, dS = \int_K^\infty (S - K) \left( -\frac{\partial}{\partial S}\left[ r S \, p \right] + \frac{1}{2} \frac{\partial^2}{\partial S^2}\left[ \sigma^2(S, T) S^2 p \right] \right) dS$$

The left side is straightforward. Differentiate the call price with respect to $T$, then solve for the integral:

$$\int_K^\infty (S - K) \frac{\partial p}{\partial T} \, dS = e^{rT} \left( \frac{\partial C}{\partial T} + r C \right)$$

The right side involves partial derivatives in $S$, so we evaluate both terms by integration by parts.

The drift term takes one pass, giving

$$\int_K^\infty (S - K) \left( -\frac{\partial}{\partial S}\left[ r S p \right] \right) dS = r e^{rT} C - r K e^{rT} \frac{\partial C}{\partial K}$$

Boundary terms vanish at $S = K$ because the payoff is zero there, and as $S \to \infty$ provided the density decays sufficiently fast at infinity, as is satisfied under the usual regularity and moment conditions.

The diffusion term takes two passes, and the second collapses the integral onto the single point $S = K$, which is what isolates the local vol at one strike out of a relation that held across all levels:

$$\int_K^\infty (S - K) \cdot \frac{1}{2} \frac{\partial^2}{\partial S^2}\left[ \sigma^2 S^2 p \right] dS = \frac{1}{2} \sigma^2(K, T) K^2 e^{rT} \frac{\partial^2 C}{\partial K^2}$$

Equating the two sides, the $rC$ terms cancel:

$$\frac{\partial C}{\partial T} = -r K \frac{\partial C}{\partial K} + \frac{1}{2} \sigma^2(K, T) K^2 \frac{\partial^2 C}{\partial K^2}$$

Rearranged, this is the **Dupire formula**:

$$\sigma^2(K, T) = \frac{\dfrac{\partial C}{\partial T} + r K \dfrac{\partial C}{\partial K}}{\dfrac{1}{2} K^2 \dfrac{\partial^2 C}{\partial K^2}}$$

Beyond the derivation we just did above, why does the local variance take on such a form? I find it more intuitive to see it in discrete form:

$$e^{rT} \Delta C \approx r K \, \Delta T \cdot \mathbb{Q}(S_T > K) + \frac{1}{2} \sigma^2(K, T) K^2 \, \Delta T \cdot p(K, T)$$

Here $\Delta C$ is the extra value of the call from extending its maturity from $T$ to $T + \Delta T$, and it comes from two sources. The first is the strike: paying $K$ later saves interest of about $rK\,\Delta T$, and that saving only matters if the option is exercised, with probability $\mathbb{Q}(S_T > K)$. The second is the extra variance. A call's payoff is flat on either side of the strike and bends only at $K$, so extra variance adds value only when the spot is near the strike at maturity. That is why the term is weighted by the density at $K$, and why only the local vol at $(K, T)$ appears: the way call prices change around one strike and one maturity reveals the local vol at that single point.

The derivation does not need the rate to be constant. With a deterministic rate $r(t)$, the same steps go through with $e^{-rT}$ replaced by the discount factor $e^{-\int_0^T r(u)\,du}$, and Dupire's formula keeps its form with $r$ replaced by $r(T)$, the rate at the option's maturity.

## Calibration in Practice

Dupire gives us a powerful tool to compute the local vol at every $(K, T)$ from call prices, with no iteration and no joint solve needed. To apply it directly, we need to know how the call price changes for a small change in $K$ and $T$. That is not available from the market, since quotes are not given at continuous $K$ or $T$. The most common way to fill this gap is to fit an implied vol surface, which translates any $K$ and $T$ into a call price. This is why we introduced the implied vol surface in the earlier article first, and it is the route we follow here.

A surface is usually displayed against strike or delta and maturity, so it might be tempting to differentiate against $K$ and $T$. But they are not the natural coordinates to view the vols in. What drives the dispersion of the price is the variance, and that variance is measured on the return rather than on the absolute price level. So the quantities the problem actually depends on are total variance and log-moneyness, and that is what we will use to derive the local vol form from the vol surface.

With $x = \ln(F_T/K)$ the log-moneyness and $w(x, T)$ the total implied variance, substitute the Black-Scholes formula for $C$ and apply the chain rule to convert the strike and maturity derivatives into derivatives of $x$ and $w$. We get

$$\sigma^2(x, T) = \frac{\dfrac{\partial w}{\partial T}}{1 - \dfrac{x}{w}\dfrac{\partial w}{\partial x} + \dfrac{1}{4}\left( -\dfrac{1}{4} - \dfrac{1}{w} + \dfrac{x^2}{w^2} \right)\left( \dfrac{\partial w}{\partial x} \right)^2 + \dfrac{1}{2}\dfrac{\partial^2 w}{\partial x^2}}$$

This is **Dupire in Gatheral's form**, and it provides a bridge to compute the local vol from the implied vol surface. The numerator is how the total variance changes with respect to time. The denominator appears convoluted but has a simple meaning: it is the ratio of the modeled density from the implied vol surface against the baseline model, the constant vol lognormal, which is the function $g(x)$ from [Appendix A]({{< ref "vol_surface_calibration.md#appendix-a-sketch-of-the-total-variance-butterfly-condition" >}}) of the implied vol article. When the vol is constant, $\partial w / \partial x$ and $\partial^2 w / \partial x^2$ are both zero, the denominator is $1$, and we are left with $\sigma^2 = \partial w / \partial T$, which is the constant variance itself.

### Differentiating the Implied Vol Surface

The Gatheral formula needs three derivatives of $w$, and the naive way to get them is to bump $x$ and $T$ and take finite differences. It works, but it is not the default way to do it. Differencing amplifies the noise in the implied vol fit, especially for the second derivative in the denominator. The result is a local vol grid that looks jagged even when the input surface is smooth, and the pricing model, Monte Carlo or PDE, inherits that noise.

For $\partial w / \partial x$ and $\partial^2 w / \partial x^2$, analytical derivatives are usually preferred, since they give a closed form on a smooth calibration. The cubic spline from the implied vol article is piecewise polynomial, so its first and second derivatives are available in closed form on each segment and are continuous across the knots by construction. Parametric models such as SSVI (Gatheral and Jacquier, 2014) are also differentiable in closed form, and parametric surfaces are often the input used in practice.

The maturity derivative, $\partial w / \partial T$, works similarly. The analytical form follows from the total variance interpolation set by the implied vol model. One nuance is worth noting. If total variance is linear in $T$ at fixed $x$, the derivative is piecewise constant and undefined at the quoted expiries themselves, since the intervals on either side give different slopes. On each interval it is well defined, and since variance accumulates as an integral over time, the value assigned at an isolated point does not affect any price. The usual convention is to take the slope of the interval ahead, which leaves the local vol surface right-continuous in maturity and jumping as the simulation crosses each quoted expiry.

The analytical approach tells us that in order to use the implied vol surface as an input, instead of just saving the vol at all the grid points, we should store the model itself, spline coefficients or parameters, and take the derivatives from it directly. Local vol calibration still runs downstream of the implied vol fit, but what crosses between them is the model rather than a set of grid points.

In some settings, passing the model along is not possible. Traders may mark the surface directly, sometimes with manual adjustments, or the marking system and the calibration library may be separate enough that only a grid of vols crosses between them. In that case a surface has to be fitted to the grid before anything can be differentiated. The risk here is that the refitted model may differ from the one underlying the original marked vol surface.

### Is Your Input Surface Good Enough?

Not every smooth looking implied vol surface is suitable for local vol calibration. Local variance is a ratio, so it exists and is positive only when the numerator and denominator are, and what makes a surface a good input turns out to be exactly the arbitrage conditions covered in the [implied vol article]({{< ref "vol_surface_calibration.md" >}}).

**The numerator requires no calendar arbitrage.** It is positive when total variance is strictly increasing in maturity at fixed log-moneyness, which is precisely the calendar condition. A violation returns a negative local variance and the calibration fails outright at that point. This is a reasonable requirement for a single underlying such as an equity vol surface. For commodities, where each contract month is a separate underlying, it may not hold strictly, which is the subject of a [later section](#adapting-local-vol-to-commodities).

**The denominator requires no butterfly arbitrage.** To have a density at all, we need the first and second derivatives of $w$, so the surface must be twice differentiable in $x$. In the implied vol article we used a cubic spline in $\sigma$ against log-moneyness, which satisfies the condition by the chain rule, since squaring a twice differentiable function leaves it twice differentiable. Had we used a simple linear interpolation in $\sigma$ instead, the slope of $w$ would jump at every knot, so its second derivative, and therefore the denominator of Gatheral's formula, would spike there. Given that the density exists, requiring it to be positive is exactly the [no butterfly arbitrage condition]({{< ref "vol_surface_calibration.md#negative-density-and-butterfly-arbitrage" >}}).

This is also where the wing treatment matters. Beyond the last anchor there are no quotes to constrain the fit, so the extrapolation we choose is what decides $g$ in the far wings. In the implied vol article we settled on total variance growing linearly in $x$, with a slope consistent with Lee's asymptotic bound, and checked that the resulting density remained positive. The local variance in the wings therefore does not turn negative or blow up.

An implied vol surface that is twice continuously differentiable in $x$ makes $g$ exist, and with no arbitrage the local vol exists and is continuous, but that does not make it smooth. On each segment of the cubic spline the curvature is linear in $x$, so while it stays continuous across a knot, its slope changes there. The local vol is built from that curvature and inherits the same behavior: it is continuous, but has a small kink at every knot. The Greeks are derivatives of the price, and through the local vol they pick up its slope, which jumps at each knot, so they can show small bumps near the knots. A parametric surface such as SSVI is smooth to every order, so its local vol has no kinks and the Greeks are smooth, which is one reason a smooth parametric surface is often preferred as the input for local vol calibration.

### Precomputing the Grid

Evaluating Gatheral's formula inside the pricing loop would mean recomputing the local vol at every step. It is more efficient to compute it once on a grid in $(x, T)$ and store it. A PDE pricer needs it at every space node on every time step, and those points are fixed in advance, so the grid can be aligned with the pricer's mesh and read off directly. In space it sits on the pricer's nodes; in time, if the Crank-Nicolson method is used, it sits at the midpoints of the time steps. Monte Carlo evaluates at the points the paths reach, so the grid is sized to cover the region they will visit and values are interpolated between nodes. One of the interpolation choices is bilinear interpolation on local variance: it is cheap enough to call millions of times, and it does not overshoot. Each interpolated value stays within the range of the surrounding nodes.

How fine the grid needs to be depends on the accuracy of the pricing scheme it feeds. The price carries two errors, one from interpolating the local vol between nodes and one from the pricing scheme's own time and space stepping. If the first is the larger, refining the pricing scheme's time step will not improve the accuracy, since the local vols driving it are already off by more than the stepping error being removed. Whether that is the case can be settled empirically: halve the spacing of the local vol grid, reprice, and see whether the price and Greeks move by more than the error tolerance set for the model. If they do not, the local vol grid is fine enough and the remaining error belongs to the pricing scheme.

### An Alternative: Calibrating Directly From Option Prices

The route above fits a smooth implied vol surface first and differentiates it, so the quality of the local vol depends on the quality of that surface. Andreasen and Huge (2011) propose an alternative that calibrates the local vol directly from option prices, without going through an implied vol surface.

Both methods rest on the same Dupire equation, which links call prices across strikes and expiries to the local vol. The difference lies in how the partial derivatives in the formula are evaluated. Dupire's formula, and Gatheral's version of it, take them directly from the smooth fitted surface. Andreasen and Huge replace them with finite differences on a grid, so each step of their method is Dupire's formula in discrete form. They guess a local vol, step the prices forward with the discretized equation, and adjust the guess until the model prices match the market. They assume a local vol that is piecewise constant in strike and work from the shortest expiry, one expiry at a time. The original paper uses a single fully implicit step between two expiries, which is fast, but the fitted vol then reprices vanilla options only within that one-step scheme. Used in a finer PDE or a Monte Carlo simulation, it would no longer match the market quotes. Taking more steps brings it closer to the true local vol and remains cheap.

With Dupire's formula, the local variance is computed from the slopes and curvature of the fitted surface, so any arbitrage in that surface can turn the result negative. Andreasen and Huge reverse the roles: the local variance is a guessed input, kept positive, and the prices are the output. Prices generated by a positive variance describe a genuine diffusion, so they rise with maturity and stay convex in strike, and the fully implicit step preserves this in discrete form. The surface is therefore free of arbitrage by construction, which makes the method more robust than differentiating a fitted surface when quotes are sparse or noisy. I plan to cover this method in a future article.

## Vanilla Option Risk Under Local Vol and Black-Scholes

Now that we know how to calibrate local vols, we can put the model to work. We start from an arbitrage-free implied vol surface, use it to calibrate the local vols, and compare the two models on the simplest product first, a vanilla call. Both models must give it the same price, but we want to compare their risk.

The input vol surface is hypothetical and carries a realistic skew. It is constructed from the SSVI model. The details of the model are not important for this exercise, only that it has closed forms for $\partial w/\partial x$, $\partial^2 w/\partial x^2$ and $\partial w/\partial T$, and that it is arbitrage-free by construction, which lets us calibrate the local vols without any problems.

The calibrated local vols and the implied vols they came from are shown below, as surfaces and as the one year skew. The surfaces can be rotated for a better view.

{{< iv_lv_comparison >}}

One immediate difference is that the local vol skew is steeper: from the forward at 105 down to a strike of 70, the implied vol rises by $6.8$ points while the local vol rises by $13.4$. The slope of the local vol is roughly two times the slope of the implied vol, and that is not a coincidence. The implied vol is the single constant we solve for so that Black-Scholes returns the market price. Since the model holds that vol fixed for the option's whole life, the number we recover stands in for all the volatility the spot experiences along the way, which makes it an average of the local vols on the path it takes. The relationship is only roughly one half for two reasons. The average in the argument is not the arithmetic average of the vols; it is the harmonic average, taken over the levels between the forward and the strike. And the relation itself is exact only in the short maturity limit, whereas the skews compared here are at one year. The derivation of the factor of two, and where the approximations enter, is in [Appendix A](#appendix-a-why-the-local-vol-skew-is-roughly-twice-the-implied-vol-skew).


### Same Option Price, Different Distributions

The contract is a one year call struck at 100, with spot at 100 and a risk-free rate of $5\%$. Under constant vol it is priced at the implied vol of the 100 strike, $30.70\%$, and Black-Scholes gives a price of $14.497$. Under local vol there is no closed form. For a single asset the pricing equation is one dimensional, and a PDE is usually preferred over simulation: it is fast and precise, and the Greeks we want to compare come straight from the grid. The local vol pricer solves the backward equation on a grid in $\ln S$ with Crank-Nicolson time stepping, started with a few implicit steps to damp the kink in the payoff, as explained in the [finite difference article]({{< ref "fdm.md" >}}). The local vol at each node is read from the precomputed grid.

<details style="margin-bottom: 1.5em;">
<summary>Local vol PDE pricer</summary>

```python
import numpy as np
from scipy.linalg import solve_banded

def lv_call_pde(S0, K, T, r, local_var, n_space=1201, n_time=500, n_rannacher=4):
    """European call under local vol: Crank-Nicolson in y = ln S, marching back from T."""
    y = np.linspace(np.log(S0 / 5), np.log(S0 * 5), n_space)
    dy, S = y[1] - y[0], np.exp(y)
    V = np.maximum(S - K, 0.0)                      # payoff at maturity
    times = np.linspace(T, 0.0, n_time + 1)
    for n in range(n_time):
        t_hi, t_lo = times[n], times[n + 1]
        dt = t_hi - t_lo
        v = local_var(S, 0.5 * (t_hi + t_lo))       # local variance at each node
        a = 0.5 * v / dy**2                         # diffusion weight from the second difference in y
        b = (r - 0.5 * v) / (2 * dy)                # drift weight from the first difference in y
        lo, di, up = a - b, -2 * a - r, a + b      # weights on V[i-1], V[i] and V[i+1]
        th = 1.0 if n < n_rannacher else 0.5        # implicit start-up, then Crank-Nicolson
        rhs = V.copy()
        rhs[1:-1] += (1 - th) * dt * (lo[1:-1] * V[:-2] + di[1:-1] * V[1:-1] + up[1:-1] * V[2:])
        ab = np.zeros((3, n_space)); ab[1] = 1.0
        ab[0, 2:] = -th * dt * up[1:-1]
        ab[1, 1:-1] = 1 - th * dt * di[1:-1]
        ab[2, :-2] = -th * dt * lo[1:-1]
        rhs[0], rhs[-1] = 0.0, S[-1] - K * np.exp(-r * (T - t_lo))   # boundary values
        V = solve_banded((1, 1), ab, rhs)
    return np.interp(np.log(S0), y, V)
```

</details>


The two models price it consistently, as they must: the local vol PDE returns $14.497$, matching the Black-Scholes price to within $0.0001$. Although the price is the same, the two models imply different distributions for the spot at one year. Both carry the same mean, the forward at $105.13$, since both models are risk neutral, but they distribute the mass differently. The market implied density, which local vol reproduces, is heavier in the far left tail ($< 55$), which is what the downside skew prices. This is typical of equity markets, where downside protection is in higher demand and is priced richer, which lifts the implied vols on that side. With the same mean, that extra mass from the left tail has to be offset, so the market density is lighter on the left shoulder and in the far right tail, and heavier just above the forward, between about 100 and 150.

{{< terminal_density >}}

For the call, what matters is the density weighted by the payoff. The call price is the probability of finishing in the money times the discounted average payoff on those paths. With more mass just above the strike, local vol finishes in the money more often, but with less in the far right tail it pays less on those paths. The right panel shows the density difference weighted by the payoff: two lobes of equal area, so what local vol gains just above the strike it gives back in the far tail, where each path pays far more, and the call price is unchanged.

|  | constant vol | local vol | difference |
|---|---|---|---|
| P($S_T > 100$) | 50.4% | 56.3% | $+5.9$ pts |
| $e^{-rT}$ E[payoff \| in the money] | 28.77 | 25.77 | $-3.01$ |
| call price | 14.497 | 14.497 | $0.000$ |


### Deriving the Greeks Under Each Model

Under Black-Scholes, the call is quoted through the implied vol at its strike, $C = C_{BS}(S, \sigma_K)$. The Black-Scholes Greeks are the derivatives of $C_{BS}$ with $\sigma_K$ held fixed:

$$\Delta_{BS} = N(d_1), \qquad \Gamma_{BS} = \frac{\varphi(d_1)}{S \sigma_K \sqrt{T}}, \qquad \text{Vega} = S \varphi(d_1) \sqrt{T}, \qquad \Theta_{BS} = -\frac{S \varphi(d_1) \sigma_K}{2\sqrt{T}} - r K e^{-rT} N(d_2)$$

One immediate difference is how delta is defined. When the spot moves, the implied vol at strike $K$ usually moves with it. The partial derivative of the Black-Scholes formula, $\partial C_{BS}/\partial S$, holds the vol fixed, which assumes the vol at $K$ does not move. Local vol makes no such assumption and allows the vol at strike $K$ to depend on the spot, so we write it $\sigma_K(S)$. The delta that actually hedges the option is the total derivative,

$$\Delta = \frac{\partial C_{BS}}{\partial S} + \text{Vega}\, \frac{\partial \sigma_K}{\partial S}$$

To see how $\sigma_K$ responds to the spot, write the vol at a strike near the money as the at-the-money vol plus the skew times the distance from the money:

$$\sigma_K \approx \sigma_{ATM} + \text{skew} \cdot \ln(K/F), \qquad \text{skew} = \frac{\partial \sigma}{\partial \ln K}$$

Here $\sigma_{ATM}$ is the vol at the forward $F$, and the skew is the slope of the smile in log strike.[^skew] When the spot moves, the forward moves with it, $d\ln F = d\ln S$, so the distance between the fixed strike and the forward $F$ changes. The total delta needs $\partial \sigma_K / \partial S$. Since the skew is already defined in log strike, it is simpler to work in log spot as well, and the two derivatives are related by

$$\frac{\partial \sigma_K}{\partial S} = \frac{1}{S}\,\frac{\partial \sigma_K}{\partial \ln S}$$

With the strike fixed and the skew treated as constant near the money,

$$\frac{\partial \sigma_K}{\partial \ln S} = \frac{d\sigma_{ATM}}{d\ln S} + \text{skew} \cdot \frac{\partial \ln(K/F)}{\partial \ln S} = \frac{d\sigma_{ATM}}{d\ln S} - \text{skew}$$

This tells us the move in the vol at a fixed strike $K$ is driven by how the at-the-money vol moves, less the skew effect. The skew term is known from today's smile, so the only unknown is the at-the-money move. The skew is also the natural yardstick for it: if the at-the-money vol moves by exactly one skew, the two terms cancel and the vol at $K$ stays put. Practitioners therefore measure the at-the-money move in units of skew, a ratio called the skew stickiness ratio (SSR):

$$\text{SSR} = \frac{d\sigma_{ATM}/d\ln S}{\text{skew}}$$

so the vol at a fixed strike moves by

$$\frac{\partial \sigma_K}{\partial \ln S} = (\text{SSR} - 1)\,\text{skew}$$

and the delta becomes

$$\Delta_{\text{SSR}} = \Delta_{BS} + \text{Vega}\, \frac{(\text{SSR} - 1)\,\text{skew}}{S}$$

If the SSR feels abstract, it helps to relate it to three assumptions about how the smile moves. You have probably heard of the first two, and the third is the one local vol implies.

- **Sticky moneyness** assumes the vol depends only on log-moneyness $\ln(K/F)$. At the money, $\ln(K/F)$ is always $0$, so the at-the-money vol never changes, its sensitivity to $\ln S$ is zero, and $\text{SSR} = 0$.
- **Sticky strike** assumes the vol depends only on the strike $K$. A change in spot then leaves the vol at a fixed strike unchanged, so $\partial\sigma_K/\partial\ln S = 0$. From the derivation above, $\partial\sigma_K/\partial\ln S = (\text{SSR} - 1)\,\text{skew}$, and with a nonzero skew, $\text{SSR}$ must be 1.
- **Local vol** assumes a function $\sigma(S, t)$: for a given $t$, the vol at a fixed $S$ stays the same. This mirrors sticky strike, where each strike keeps its implied vol. In the short maturity limit the at-the-money implied vol and the local vol at the spot are the same ([Appendix A](#appendix-a-why-the-local-vol-skew-is-roughly-twice-the-implied-vol-skew)), so the at-the-money vol moves along the local vol curve, just as it moves along the implied smile under sticky strike. The only difference is the slope. Under sticky strike the slope is the implied skew, which gives $\text{SSR} = 1$, as explained in the previous bullet. The local vol skew is twice the implied skew (also in Appendix A), so the at-the-money vol moves twice as far, and $\text{SSR} \to 2$ as $T \to 0$.

In equity markets, sticky strike is a common assumption when the Black-Scholes model is applied, so let's compare the Greeks of the one year at-the-money call under sticky strike and under local vol. The local vol Greeks come from the [PDE pricer]({{< ref "fdm.md" >}}): we bump the spot with $\sigma(S, t)$ held fixed, and for vega we bump the implied vol surface in parallel and recalibrate.

|  | Black-Scholes<br>(sticky strike) | local vol |
|---|---|---|
| delta | 0.624 | 0.534 |
| gamma | 0.0124 | 0.0252 |
| vega (per vol point) | 0.379 | 0.379 |
| theta (per year) | $-8.22$ | $-12.21$ |

Let's go through the results one by one, starting with delta. From the total delta formula, $\Delta = \Delta_{BS} + \text{Vega}\,(\text{SSR} - 1)\,\text{skew}/S$, the local vol delta differs from the sticky strike delta by the second term. Local vol has the higher SSR, and the skew at the money is typically negative for equity ($-0.148$ in our example), so the local vol delta must be lower. Hedging a long position in this call option under local vol therefore requires selling less of the underlying.

Under local vol the option has almost twice the convexity, and we can see why by differentiating the total delta once more. Since the price is $C_{BS}(S, \sigma_K(S))$,

$$
\Gamma = \Gamma_{BS} + 2\,\text{Vanna}\,\frac{\partial \sigma_K}{\partial S} + \text{Volga}\left(\frac{\partial \sigma_K}{\partial S}\right)^2 + \text{Vega}\,\frac{\partial^2 \sigma_K}{\partial S^2}.
$$

Under sticky strike $\sigma_K$ does not depend on $S$, so only $\Gamma_{BS}$ survives. Under local vol, vanna and volga are close to zero for an at-the-money option: in this example the vanna term is $0.00006$ and the volga term $0.000002$, against a Black-Scholes gamma of $0.0124$. That leaves the last term, the product of the vega, $37.95$ per unit of vol, and $\partial^2\sigma_K/\partial S^2$. The second derivative is not negligible because $\sigma_K$ is a curved function of the spot. Under local vol, $\sigma_K$ is built from the local vols that paths see, so as $t \to 0$, $\sigma_K$ inherits the shape of $\sigma(\cdot, t)$ around the spot, and there the local vol is not a straight line but bends.[^shortend] On this surface $\partial^2\sigma_K/\partial S^2 \approx 0.000328$, so the vega term is about $0.0125$, roughly the size of $\Gamma_{BS}$ itself, which is why the gamma doubles. The doubling is specific to this example, but the direction holds when the local vol is convex around the spot.

Given the differences in delta and gamma, the change in theta is easy to confirm, because by the no-arbitrage argument both models satisfy the same pricing PDE:

$$
\Theta + rS\,\Delta + \tfrac12 \sigma^2 S^2\,\Gamma = rC,
$$

where $\sigma$ is a constant under Black-Scholes and depends on $S$ under local vol, $\sigma(S, t)$. Since both models are meant to reprice the market quotes exactly, $C$ is the same, and the larger gamma term in local vol has to be balanced by a more negative theta. The impact of delta is small compared with that of gamma.

Finally, the vega is the same in both columns. A parallel shift of the implied vol surface moves the vanilla price by the same amount under both models: Black-Scholes simply takes the bumped vol and plugs it into the formula, and local vol, recalibrated to the bumped surface, reprices the vanilla at that same value.

### P&L Attribution After a Spot and Vol Move

After making sense of the Greeks, let's put them to use in explaining the P&L after a market move. Suppose that over one trading day the spot falls $2\%$, from 100 to 98, and implied vols rise, as they usually do in a sell-off. Say the smile keeps its shape while its at-the-money vol moves with the spot at $\text{SSR} = 1.5$, and on top of that the whole surface rises by one vol point. This translates to a rise of $1.15$ points in the implied vol at the 100 strike, and the call loses $0.821$. We choose $\text{SSR} = 1.5$ because a value of 1 or 2 would favour one of the two models, and because the S&P 500 and Euro Stoxx 50 option markets have historically shown an SSR of around 1.5 for maturities from a month to a few years (Bergomi, 2016).

We show the P&L attribution under the two models in the table below.

<table>
<thead>
<tr><th style="border-bottom:1px solid #888"></th><th colspan="2" style="text-align:center;border-bottom:1px solid #888">sticky strike</th><th colspan="2" style="text-align:center;border-bottom:1px solid #888">local vol</th></tr>
<tr><th></th><th style="border-left:2px solid #888;">formula</th><th>P&amp;L</th><th style="border-left:2px solid #888;">formula</th><th>P&amp;L</th></tr>
</thead>
<tbody>
<tr><td>delta P&L</td><td style="border-left:2px solid #888;">$\Delta\, dS$</td><td>$-1.248$</td><td style="border-left:2px solid #888;">$\Delta_{LV}\, dS$</td><td>$-1.067$</td></tr>
<tr><td>gamma P&L</td><td style="border-left:2px solid #888;">$\tfrac12 \Gamma\, dS^2$</td><td>$+0.025$</td><td style="border-left:2px solid #888;">$\tfrac12 \Gamma_{LV}\, dS^2$</td><td>$+0.050$</td></tr>
<tr><td>theta P&L</td><td style="border-left:2px solid #888;">$\Theta\, dt$</td><td>$-0.033$</td><td style="border-left:2px solid #888;">$\Theta_{LV}\, dt$</td><td>$-0.048$</td></tr>
<tr><td>vega P&L</td><td style="border-left:2px solid #888;">$\text{Vega}\, \Delta\sigma_K^{\text{market}}$</td><td>$+0.436$</td><td style="border-left:2px solid #888;">$\text{Vega}\,(\Delta\sigma_K^{\text{market}} - \Delta\sigma_K^{LV})$</td><td>$+0.262$</td></tr>
<tr><td>unexplained</td><td style="border-left:2px solid #888;"></td><td>$-0.001$</td><td style="border-left:2px solid #888;"></td><td>$-0.018$<sup>*</sup></td></tr>
<tr><td>total</td><td style="border-left:2px solid #888;"></td><td>$-0.821$</td><td style="border-left:2px solid #888;"></td><td>$-0.821$</td></tr>
</tbody>
</table>

<small><sup>*</sup> The small unexplained term under local vol comes from truncating the expansion: $\sigma_K$ is a curved function of the spot, so delta and gamma miss slightly more than under sticky strike. It disappears under full revaluation.</small>

The delta, gamma and theta P&L are fairly self-explanatory, but the vega P&L needs some clarification. We saw earlier that vega is the same in both models, so why is the vega P&L different? Recall how we calculated the vega of $0.379$: we applied a parallel shift to the entire vol surface. The market move here is not purely parallel. Of the $1.15$ points at the 100 strike, $1$ point is the parallel shift, which neither model expects. Sticky strike expects no vol change at all. Local vol, with $\sigma(S, t)$ held fixed, can move implied vols only through the spot, so a rise in vol that has nothing to do with the spot move can only come from recalibrating the surface. For this part of the move, the two models record an identical vega P&L of $0.379$. The other $0.15$ points comes from the smile moving with the spot, and here the models disagree. The vega P&L of each model reflects only the part of the vol move that the model did not already expect from the spot move:

$$
\Delta\sigma = \Delta\sigma_K^{\text{market}} - \Delta\sigma_K^{\text{model}}.
$$

Sticky strike expects no change in the vol at a fixed strike, so the whole $1.15$ points is a surprise. Local vol is different. With $\sigma(S, t)$ held fixed, repricing the call from a spot of 98 gives an implied vol at the 100 strike that is $0.46$ points higher, so this rise is what local vol expects from the spot move alone. Recall that the local vol delta is $\Delta_{LV} = \Delta_{BS} + \text{Vega}\,\partial\sigma_K/\partial S$. Its delta P&L therefore already contains $\text{Vega} \times 0.46 \approx 0.17$, the value of the expected vol rise, which is why the local vol delta loses about $0.17$ less than the sticky strike delta. Counting the same $0.46$ points again in vega would double count it, so only the remaining $0.69$ points of vol movement is a surprise to the local vol model.

The two models agree on the total P&L, since it is set by the market, but they tell different stories. Sticky strike says the option lost $1.25$ on delta and the vol rally recovered $0.44$. Local vol says the option lost less on delta, because a vol rise was expected in a sell-off, and only a smaller vol surprise of $0.26$ was left for vega. So which model is more consistent with market behaviour? As it turns out, neither. As noted above, markets have historically shown an SSR of around 1.5, so the real smile moves somewhere between what sticky strike and local vol assume. This means a desk could build its own SSR assumption into risk management, hedging with a delta based on an SSR estimated from market history instead of the value either model implies:

$$
\Delta_{\text{desk}} = \Delta_{BS} + \text{Vega}\,(\text{SSR}_{\text{est}} - 1)\,\frac{\text{skew}}{S}.
$$


This example uses Greek-based attribution. In practice, desks can also shock one market factor at a time and fully revalue the option after each shock. This removes the small unexplained term, but it does not change the conclusion. The takeaway from this exercise is that sticky strike and local vol read the same market movements differently. They agree on the total P&L, but they disagree on how much of it is a response to the spot move and how much is a genuine move in vol. Under local vol, part of the vol rise is the expected reaction to the spot falling, so it sits in delta; under sticky strike, all of it sits in vega. This is why the two models should not be mixed for risk management. Exotic options are often hedged with vanilla options. If the exotic is valued under local vol while the vanillas are valued under sticky strike, the hedge ratios are off and the book shows more P&L noise. Both products should be valued within the same model, with adjustments to the vanilla hedge ratios reflecting the desk's assumptions about smile dynamics.

## How Does Local Vol Change the Barrier Option Price?

The vanilla example in the previous section served to build intuition for the local vol model, but vanillas are not where the model earns its place. Its natural application is to exotic options, whose payoffs depend on the entire path of the underlying rather than on its terminal value alone. The vanillas pin down the distribution of the spot at each maturity, but not how the spot travels between them, and that is where the two models can disagree. Let's take a barrier option as an example and analyse where the price difference between the two models comes from.

The contract is the at-the-money call from the previous section with an up-and-out barrier at 130 added:

| spot | strike | barrier (up-and-out) | maturity | risk-free rate |
|---|---|---|---|---|
| 100 | 100 | 130, monitored daily | 1 year | $5\%$ |

Its value comes only from paths that end between the strike and the barrier without touching 130, and the payoff jumps from nearly 30 to nothing at the barrier.

For a single asset barrier like this one, a PDE remains the standard tool: the knock-out is just a boundary condition at 130. For barriers on several assets the grid grows with each dimension and [Monte Carlo]({{< ref "mc_variance_reduction.md" >}}) becomes the practical choice. Here we go with Monte Carlo since tracking individual paths shows directly how often they survive and where they finish, which is what we want to explore. We price the option under constant vol and local vol with the same MC settings, so that no part of the gap comes from the numerical scheme: 100,000 paths, 252 time steps, and the same seed. Under constant vol, every path uses the implied vol of the 100 strike, $30.70\%$, since the barrier option is a call struck at 100 made conditional on survival. Under local vol, the grid is precomputed in $(x, t)$ at a spacing of $0.005$ in log-moneyness and $0.001$ years in time, and read by bilinear interpolation on local variance. The spacing follows the refinement test set out in [Precomputing the Grid](#precomputing-the-grid): halving it moves the price by less than the MC standard error.

The barrier option prices are significantly different: constant vol prices it at $1.619$ and local vol at $2.348$, $45\%$ higher!

The probability of survival differs by only $1.3\%$ between the two models, so the difference in barrier survival alone cannot account for a price gap of this size. The gap has to be driven largely by the terminal payoff. To see this, we plot the terminal density of the surviving paths under each model.[^density]

{{< barrier_density >}}

Above the barrier, $130$, the option is knocked out, so the surviving density is zero there. Well below the barrier the comparison looks much like the vanilla case: paths that end far below 130 are unlikely to have been knocked out, so local vol again puts more mass in the far left tail, less on the left shoulder, and more between the strike and the barrier. Close to 130 both densities are pulled down towards zero, since paths ending there have most likely touched the barrier along the way.

The right panel weights the density difference by the payoff. For the vanilla, it showed two lobes of equal area: what local vol gained just above the strike it gave back in the far right tail. The barrier removes that tail, since any path finishing above 130 is knocked out and pays nothing. What remains is the band between the strike and the barrier, exactly where local vol has the extra mass, so the whole area is positive. That area is the discounted price gap, $2.348 - 1.619 \approx 0.73$, and it combines two effects shown in the table below: the paying state is $5.8\%$ more likely under local vol, and the discounted average payoff in it is higher as well, by $0.70$.

|  | constant vol | local vol | difference |
|---|---|---|---|
| barrier price | 1.619 | 2.348 | $+45.0\%$ |
| standard error | 0.015 | 0.017 | |
| P(survive) | 62.5% | 63.8% | $+1.3$ pts |
| P(survive and $S_T > 100$) | 16.5% | 22.3% | $+5.8$ pts |
| $e^{-rT}$ E[payoff \| survive, in the money] | 9.82 | 10.52 | $+0.70$ |

The size and even the sign of this gap belong to this particular example. A different contract such as a down-and-out put would give a different comparison. Still, the example provides a general lesson: an exotic should be priced with a model that reproduces the smile. Black-Scholes with the constant vol of a single strike implies a lognormal distribution at maturity, which differs from the one the market implies, and most of the gap above comes from that alone.

## Local Vol Matches the Spot's Marginal Distributions, Not Its Path

The barrier example leaves an open question. A model that matches the market vol smile is the minimum requirement, and local vol meets it. But is that enough to capture how the spot moves between the quoted maturities? We can confirm that local vol matches the smile by checking that it reprices the market's vanilla options at each maturity, but those prices say nothing about the vol dynamics in between. To judge whether local vol models those dynamics realistically, we need a different perspective: compare what the model assumes about vol with how vol actually behaves in the market.

Local vol assumes the vol is a deterministic function of spot and time, so the vol carries no randomness of its own. The market does not treat vol that way. The VIX is the market's implied vol for the SPX, and options on the VIX are themselves quoted with an implied vol, a vol of vol, so the market prices vol itself as a random quantity. The return and the vol also tend to be negatively correlated in many markets: when the price drops, vol often rises, yet the size of the rise is not fixed by the size of the drop. Both observations point towards giving vol randomness of its own, correlated with the returns, which is the idea behind [stochastic vol models]({{< ref "heston.md" >}}).

Before discrediting the local vol model immediately due to its lack of a second driver, let's ask what simplifying assumption it is making. Suppose the true underlying dynamics is that the spot follows a diffusion whose instantaneous volatility is some process $\beta_t$,

$$dS_t = r S_t \, dt + \beta_t \, S_t \, dW_t$$

$\beta_t$ is unknown to us. It may be driven by its own Brownian motion, or by several other factors. It is not determined by $S_t$ the way the local vol function is. Take $f$ to be any smooth function of the spot. Ito's lemma gives

$$df(S_t) = \left( r S_t f'(S_t) + \frac{1}{2} \beta_t^2 S_t^2 f''(S_t) \right) dt + \beta_t S_t f'(S_t) \, dW_t$$

and taking expectations kills the stochastic integral, leaving

$$\frac{d}{dt}\mathbb{E}\left[ f(S_t) \right] = \mathbb{E}\left[ r S_t f'(S_t) \right] + \frac{1}{2} \mathbb{E}\left[ \beta_t^2 S_t^2 f''(S_t) \right]$$

The first term involves only $S_t$. The second couples $\beta_t$ to $S_t$, and this is where the unknown process has to be dealt with. Applying the tower property, only the average of $\beta_t^2$ at each level survives once it is multiplied by a function of the spot alone, so

$$\frac{d}{dt}\mathbb{E}\left[ f(S_t) \right] = \mathbb{E}\left[ r S_t f'(S_t) + \frac{1}{2} \mathbb{E}\left[ \beta_t^2 \mid S_t \right] S_t^2 f''(S_t) \right] \tag{1}$$

Now run the same expectation in the local vol model, where the volatility is the deterministic function $\sigma(S_t, t)$. Ito gives the same expression with $\sigma^2(S_t, t)$ in place of $\beta_t^2$, and no conditioning is needed since it is already a function of the level:

$$\frac{d}{dt}\mathbb{E}\left[ f(S_t) \right] = \mathbb{E}\left[ r S_t f'(S_t) + \frac{1}{2} \sigma^2(S_t, t) S_t^2 f''(S_t) \right] \tag{2}$$

The expectation in (1) is taken over the marginal density of $S_t$ under the true risk-neutral dynamics, and the one in (2) is taken over its marginal density under local vol. Local vol is calibrated to reprice the market's vanilla options, and by [Breeden-Litzenberger]({{< ref "vol_surface_calibration.md#breeden-litzenberger" >}}) the density is recovered directly from those prices by differentiating twice in $K$, so the density under local vol must match the one under the true dynamics at every maturity. Both expectations are then taken over the same density, and subtracting one equation from the other leaves

$$
\mathbb{E}\left[ \tfrac{1}{2} \left( \mathbb{E}\left[ \beta_t^2 \mid S_t \right] - \sigma^2(S_t, t) \right) S_t^2 f''(S_t) \right] = 0
$$

for every smooth $f$. Since $f''$ can be any smooth function, the difference $\mathbb{E}[\beta_t^2 \mid S_t] - \sigma^2(S_t, t)$ must be zero at every level the spot can reach, which gives

$$
\sigma^2(S, t) = \mathbb{E}\left[ \beta_t^2 \mid S_t = S \right]
$$

This is the **Gyöngy** result. No matter what the true volatility process is, the local vol model reproduces its marginal density at every maturity, and the local variance it needs, $\sigma^2(S, t)$, is the conditional average of the true instantaneous variance at that level and time, $\mathbb{E}[\beta_t^2 \mid S_t = S]$. In the barrier example, at $S = 130$ one year out the calibrated local vol is $26.1\%$, or $0.068$ in local variance. This value is sometimes misinterpreted as a forecast: if the spot reaches 130, this is the variance it will have there. Gyöngy gives a clear interpretation. It is the conditional average of the true instantaneous variance across all volatility states that can occur when the spot is at 130 at that time. Paths that reach 130 one year out arrive there in different conditions, some after a calm stretch with a low vol and some after a turbulent one with a high vol, and the local variance averages over all possible states, weighted by how likely each is.

Although the local vol model matches the marginal distribution by keeping the mean of the instantaneous variance, it throws away the distribution of that variance at a given level and time. A single asset barrier is path dependent in a relatively simple way, since its value depends on whether a particular level is breached, and local vol is still used often in practice for such products. A forward starting option is different, because its value depends on the smile seen from a future date rather than from today. In local vol, the vol used at every future level and time is already fixed by today's surface. The figure below plots slices of the local vol surface from the barrier example at 1M, 6M and 1Y.

{{< forward_skew >}}

The slices flatten quickly. The local vol at 90% of the forward exceeds the one at 110% by 12.0 vol points at 1M, 6.3 vol points at 6M and 4.7 vol points at 1Y. Today's 3M smile is built from the steep early slices. The same 3M option starting in a year only experiences the slices from 1Y onward, and since implied vol is an average of the local vols along the path, its smile can only come out much flatter than today's. One way to fix this is to give the vol randomness of its own, so that the skew can regenerate itself over time. This leads to the stochastic vol model, which we cover in the [Heston article]({{< ref "heston.md" >}}).

## Adapting Local Vol to Commodities

So far the article has assumed an equity setting, where the vol surface holds the vols of a single underlying seen at a range of maturities. Commodity options work differently. They are written on futures, and each contract month is its own underlying with its own futures price, so moving along the maturity axis changes the underlying as well as the horizon. The December contract and the following June contract respond to different supply and demand, and an option on one is not a direct statement about the other.

That can be a problem for Dupire, since the numerator of the local variance is $\partial w / \partial T$. Taking that derivative on an equity surface compares the same underlying at two horizons. Taking it on a commodity surface compares two different underlyings, and the difference between them reflects not only the passage of time but also the distinction between the contracts themselves.

Whether this matters in practice depends on how close the contract months are to being the same underlying. When they are highly correlated, treating the strip as a single underlying is a reasonable approximation. Precious metals are a good example: a gold futures curve is close to deterministic carry, and the contract months move almost together. For natural gas the assumption may not be appropriate. A January contract and a July contract are vastly different and not interchangeable: they are exposed to different weather, different storage positions, and different demand. Their prices are only loosely correlated, and their implied vol levels also differ. Differencing across them does not isolate the effect of time. Across a gas strip total variance is also not guaranteed to increase along the surface, and where the calendar condition breaks the formula returns a negative local variance.

Some markets do list several option expiries against a single contract, such as the Henry Hub weeklies, which are listed for four weeks. But that reaches only the nearest contract, and the rest of the strip is left with one expiry each. For most of the curve the maturity dimension therefore has to be modeled rather than read off the quotes. There is more than one way to do this, and the topic deserves its own article, so here I only summarize the main ideas of two approaches.

### Two-Factor Model

The first approach assumes a shape for how each contract's vol evolves in time. A common choice is a two factor model (Schwartz and Smith, 2000), in which the futures delivering at $T$ follows

$$\frac{dF_t(T)}{F_t(T)} = \sigma_1 e^{-\kappa (T-t)}\, dW^1_t + \sigma_2\, dW^2_t, \qquad dW^1_t\, dW^2_t = \rho\, dt.$$

The first factor captures short term shocks, which hit the nearby contracts hardest and fade along the curve at speed $\kappa$. The second factor captures long term shocks, which move all contracts in parallel. Together they give the contract delivering at $T$ a deterministic instantaneous vol

$$\sigma_T(t)^2 = \sigma_1^2 e^{-2\kappa (T-t)} + 2\rho\,\sigma_1 \sigma_2\, e^{-\kappa (T-t)} + \sigma_2^2,$$

which typically rises as the contract approaches delivery. This is often referred to as the Samuelson effect (Samuelson, 1965). Since the instantaneous vol is deterministic, the total variance accumulated up to time $t$, its integral over time, also has a closed form:

$$V_T(t) = \int_0^t \sigma_T(u)^2\, du = \frac{\sigma_1^2}{2\kappa}\left(e^{-2\kappa (T-t)} - e^{-2\kappa T}\right) + \frac{2\rho\,\sigma_1 \sigma_2}{\kappa}\left(e^{-\kappa (T-t)} - e^{-\kappa T}\right) + \sigma_2^2\, t.$$

The model produces a flat smile for each contract, since $\sigma_T(t)$ does not depend on the strike. To match the market skew, each log-moneyness $x = \ln(F_0(T)/K)$ gets a multiplier $M(x)$, chosen so that the scaled model variance up to the option expiry $T_E$ equals the implied total variance quoted there:

$$M(x)\, V_T(T_E) = \sigma_{\text{imp}}(x)^2\, T_E.$$

The same multiplier then fills the missing maturity dimension. The total implied variance at an intermediate time $t$ is taken to be $w(x, t) = M(x)\, V_T(t)$, and the local vol comes out of Dupire's formula as in the equity case: the numerator $\partial w / \partial t = M(x)\, \sigma_T(t)^2$ is driven by the two factor assumption, while the log-moneyness derivatives in the denominator are driven by the skew through $M(x)$. Since $M(x)$ and $\sigma_T(t)^2$ are both positive, the numerator can never turn negative, which removes the calendar problem described above. The drawback of the two factor model is that the result depends on the assumed dynamics, and its parameters, $\sigma_1$, $\sigma_2$, $\kappa$ and $\rho$, require careful consideration. For anyone interested in further reading, Ögetbil and Hientzsch (2022) discuss a similar idea.

### A Fictitious Spot Driving the Curve

The second approach drives all contract months with a single process, so that the one expiry each contract carries becomes one point on a maturity axis shared by all of them. Nastasi, Pallavicini and Sartorelli (2020) build this process from what they call a fictitious spot price, a process $s_t$ that is not traded and need not match any quoted spot. Futures prices $F_t(T)$ are derived from it as

$$\frac{F_t(T)}{F_0(T)} - 1 = (s_t - 1)\, e^{-\int_t^T a(u) \, du}$$

with $s_0 = 1$. The relative move of the contract delivering at $T$ equals the move in $s_t$, damped by a factor that depends on the time to delivery. At $t = 0$ the right side is zero, so the formula recovers today's curve exactly. A contract about to expire has a damping factor close to one and moves with $s_t$ one for one, which is why $s_t$ behaves like a spot, while a distant contract has a small damping factor and barely reacts. This again produces the Samuelson effect. For the futures to be martingales under the risk-neutral measure, $s_t$ must mean revert towards one at speed $a$,

$$ds_t = -a(t)\,(s_t - 1)\, dt + \eta(t, s_t)\, dW_t.$$

Instead of modelling the options on each futures contract separately, we can now model them all as options on the fictitious spot, by translating each contract's strike into a level of $s_t$. Inverting the futures formula, a strike $K$ on the contract delivering at $T$ corresponds to the level

$$s^* = 1 - \left(1 - \frac{K}{F_0(T)}\right) e^{\int_{T_E}^T a(u) \, du}$$

of the fictitious spot, where $T_E$ is the option expiry. An option on a futures contract, expiring in six months, therefore becomes an option on $s_t$ six months out, struck at $s^*$. Nastasi, Pallavicini and Sartorelli then extend Dupire's formula to the mean-reverting drift and calibrate it against the pooled set of options to give the local vol $\eta(t, s)$.

The drawback of this construction is that the whole curve moves through a single factor, so the contract months are perfectly correlated in their shocks and differ only in how strongly they are damped. The mean reversion speed $a$ also cannot be pinned down from vanilla quotes and needs to be calibrated from calendar spread options. In a future article, I plan to implement both models and compare how they behave in practice.

## Conclusions

Black-Scholes with a constant vol and local vol make different assumptions about how the vol moves with the spot. Under sticky strike, the vol at a fixed strike stays put when the spot moves; under local vol, it moves with the spot, at an SSR of about two for short maturities. The two models can agree on the price of a vanilla and still disagree on its Greeks, and so they attribute the same market move differently in P&L. This is why they should not be mixed in risk management: an exotic valued under local vol and hedged with vanillas valued under sticky strike leaves a book whose hedge ratios do not line up.

The barrier example showed that for path dependent products the minimum requirement is a model that reproduces the market's terminal density, which local vol does by construction. Gyöngy's result gives a clear interpretation of the model: the local variance is the conditional average of the true instantaneous variance at each level and time. This means local vol does not model the distribution of the instantaneous variance around that average. That is often acceptable for a single asset barrier, but not for products that depend on how the smile itself evolves, such as forward starting options. We will address this issue in the [next article]({{< ref "heston.md" >}}).

## Appendix A: Why the Local Vol Skew is Roughly Twice the Implied Vol Skew

To see where the factor of two comes from, let's start by examining how the option price behaves in the short maturity limit, where $T \to 0$.

The price of an option is set by how far the strike sits from the forward measured in standard deviations,

$$d = \frac{x}{\sigma\sqrt{T}}$$

The forward carries no drift, so $F_T = F e^{-\sigma^2 T/2 + \sigma\sqrt{T} Z}$ with $Z$ standard normal, and the call can be approximated as

$$C \approx e^{-rT} F \sigma \sqrt{T} \; \mathbb{E}\left[ (Z + d)^+ \right]$$

The details specific to the contract have been pushed into the two factors. The scale $e^{-rT} F \sigma \sqrt{T}$ is the discounted size of a one standard deviation move on the forward, which is what an option is worth per unit of standardized payoff, and the expectation is a pure function of $d$, the number of standard deviations the strike sits away. Writing it out,

$$\mathbb{E}\left[ (Z + d)^+ \right] = \int_{-d}^{\infty} (z + d)\,\varphi(z)\,dz$$

with $\varphi$ the normal density. At the money this is trivial, giving $\mathbb{E}[Z^+] = 1/\sqrt{2\pi}$, so

$$C \approx \frac{1}{\sqrt{2\pi}} \, e^{-rT} F \sigma \sqrt{T}$$

The price is proportional to the vol, and matching prices means matching vols. This shows $\sigma_{IV}(0) = \sigma_{LV}(0)$: the two must agree at the money as $T \to 0$. That tells us the vol levels, but nothing about how the vol skews of the two models relate. So we need to look at a non-trivial case, one involving strikes sitting away from the forward. Take an out-of-the-money call, where $d < 0$. The in-the-money case follows the same argument due to put-call parity. Changing variables to $z = -d + u$,

$$\int_{-d}^{\infty} (z + d)\,\varphi(z)\,dz \approx \varphi(d) \int_0^\infty u\, e^{-|d|u}\,du = \frac{\varphi(d)}{d^2}$$

so

$$C \approx e^{-rT} F \sigma \sqrt{T} \, \frac{\varphi(d)}{d^2} \qquad \Longrightarrow \qquad \ln C \approx -\frac{d^2}{2} + O(\ln|d|)$$

As $T \to 0$ the distance $|d|$ grows to infinity, so the price is dominated by the Gaussian term while the prefactor contributes only at logarithmic order. So two models agreeing on the price must agree on $d$. Under Black-Scholes the volatility is constant for a given strike, so the distance is the log-moneyness divided by the standard deviation over the option's life:

$$d = \frac{x}{\sigma_{IV}(x)\sqrt{T}}$$

Under local vol the volatility depends on the level, so the distance has to be accumulated level by level. Each stretch $du$ of log-moneyness between the forward and the strike is divided by the volatility available at that level, so the local vol analogue of $d$ is

$$d = \frac{1}{\sqrt{T}} \int_0^x \frac{du}{\sigma_{LV}(u)}$$

Equating the two, the $\sqrt{T}$ cancels:

$$\frac{x}{\sigma_{IV}(x)} = \int_0^x \frac{du}{\sigma_{LV}(u)}$$

This is the **Berestycki-Busca-Florent** result, exact in the limit $T \to 0$. It reveals that the implied vol is the harmonic average of the local vol over log-moneyness, from the forward ($u = 0$) to the strike ($u = x$). The factor of two follows immediately. Take the local vol to be linear near the money, $\sigma_{LV}(u) = \sigma_0 + bu$, and expand the reciprocal for small $bu / \sigma_0$,

$$\frac{1}{\sigma_{LV}(u)} \approx \frac{1}{\sigma_0}\left( 1 - \frac{bu}{\sigma_0} \right)$$

Averaging over $[0, x]$ gives

$$\frac{1}{\sigma_{IV}(x)} \approx \frac{1}{\sigma_0}\left( 1 - \frac{bx}{2\sigma_0} \right) \qquad \Longrightarrow \qquad \sigma_{IV}(x) \approx \sigma_0 + \frac{bx}{2}$$

The implied vol slope is $b/2$ against the local vol's $b$. Two approximations are in play. The relation holds in the short-maturity limit, and the local vol is treated as linear over the range traveled. Neither holds cleanly at the one year maturity discussed in the main text, where the local vol is convex rather than linear. The measured ratio comes out slightly below two, and the instantaneous slope at the money gives $1.8$. The factor of two is a rule of thumb, not an exact identity.

## References

Andreasen, J., and Huge, B. (2011). Volatility interpolation. *Risk*, 24(3), 76-79.

Berestycki, H., Busca, J., and Florent, I. (2002). Asymptotics and calibration of local volatility models. *Quantitative Finance*, 2(1), 61-69.

Bergomi, L. (2016). *Stochastic Volatility Modeling*. Chapman and Hall/CRC.

Breeden, D. T., and Litzenberger, R. H. (1978). Prices of state-contingent claims implicit in option prices. *Journal of Business*, 51(4), 621-651.

Buehler, H. (2007). Volatility and dividends: volatility modelling with cash dividends and simple credit risk. Working paper, TU Berlin.

Dupire, B. (1994). Pricing with a smile. *Risk*, 7(1), 18-20.

Gatheral, J. (2006). *The Volatility Surface: A Practitioner's Guide*. Wiley.

Gatheral, J., and Jacquier, A. (2014). Arbitrage-free SVI volatility surfaces. *Quantitative Finance*, 14(1), 59-71.

Gyöngy, I. (1986). Mimicking the one-dimensional marginal distributions of processes having an Itô differential. *Probability Theory and Related Fields*, 71(4), 501-516.

Lee, R. W. (2004). The moment formula for implied volatility at extreme strikes. *Mathematical Finance*, 14(3), 469-480.

Nastasi, E., Pallavicini, A., and Sartorelli, G. (2020). Smile modeling in commodity markets. *International Journal of Theoretical and Applied Finance*, 23(3), 2050019.

Ögetbil, O., and Hientzsch, B. (2022). A flexible commodity skew model with maturity effects. arXiv:2212.07972.

Samuelson, P. A. (1965). Proof that properly anticipated prices fluctuate randomly. *Industrial Management Review*, 6(2), 41-49.

Schwartz, E., and Smith, J. E. (2000). Short-term variations and long-term dynamics in commodity prices. *Management Science*, 46(7), 893-911.

[^shortend]: The local vol near $t = 0$ depends on how the surface is extended below the first quoted expiry. Here we hold the one day implied vols constant there. A different choice changes the local vol gamma and theta.

[^density]: A plain histogram of the simulated paths shows the same comparison but looks choppy, so we smooth it by using the last time step: given where a surviving path stands one day before maturity, its final step is lognormal with a known density, and averaging that density over the 100,000 paths gives a smooth estimate of the terminal density.

[^skew]: The skew is the slope of the implied vol curve against strike. Here it is measured in log strike, $\partial \sigma / \partial \ln K$, the vol change per percentage move in strike, as is common in the quant literature. Some desks use $\partial \sigma / \partial K$ instead, the vol change per dollar of strike. The two are related by $\partial \sigma / \partial \ln K = K\, \partial \sigma / \partial K$.
