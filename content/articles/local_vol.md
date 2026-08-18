---
title: "Local Volatility: From the Implied Vol Surface to Risk-Neutral Dynamics"
date: 2026-07-23
draft: False
math: true
tags: ["local volatility", "Dupire", "Gatheral", "implied volatility", "Gyongy", "barrier options", "marginal density", "commodities", "Monte Carlo"]
---

## Why This Matters

In an [earlier article]({{< ref "vol_surface_calibration.md" >}}) we constructed the implied volatility surface and used it primarily to price vanilla options. But a great deal of what trades is not vanilla. Products like barriers and autocallables depend on the path the underlying takes, not only where it lands.

Suppose I price a barrier by Monte Carlo. At each step the spot sits at some level, and I need a volatility to advance it. What vol do I use? The surface gives me a vol for every strike, but simulation does not ask about strikes. It asks what volatility the spot experiences at this level, at this moment, which the surface cannot answer.

Local volatility is one solution. It gives the instantaneous vol of the underlying at each level and time. This article explains how the model works, and is my exploration of five questions:

- How should local vol be calibrated in theory?
- What are the challenges of calibrating it in practice?
- If we price a barrier option under local vol, how do the results differ from constant vol, and why?
- What does local vol tell us about the path the underlying takes, and what does it leave out?
- Why can't local vol be applied directly to some commodities, and what adaptations does it need?

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

The right side involves partial derivatives in $S$, so we evaluate both terms by integration by parts. Boundary terms vanish at $S = K$ because the payoff is zero there, and as $S \to \infty$ provided the density decays sufficiently fast at infinity, as is satisfied under the usual regularity and moment conditions.

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

## How Local Vol is Calibrated in Practice

Dupire gives us a powerful tool to compute the local vol at every $(K, T)$ from call prices, with no iteration and no joint solve needed. To apply it directly, we need to know how the call price changes for a small change in $K$ and $T$. That is not available from the market, since quotes are not given at continuous $K$ or $T$. But we can calculate what the call price should be if an implied vol surface is in place, which translates any $K$ and $T$ into a call price. This is why we introduced the implied vol surface in the earlier article first.

A surface is usually displayed against strike or delta and maturity, so it might be tempting to differentiate against $K$ and $T$. But they are not the natural coordinates to view the vols in. What drives the dispersion of the price is the variance, and that variance is measured on the return rather than on the absolute price level. So the quantities the problem actually depends on are total variance and log-moneyness, and that is what we will use to derive the local vol form from the vol surface.

With $x = \log(F_T/K)$ the log-moneyness and $w(x, T)$ the total implied variance, substitute the Black-Scholes formula for $C$ and apply the chain rule to convert the strike and maturity derivatives into derivatives of $x$ and $w$. We get

$$\sigma^2(x, T) = \frac{\dfrac{\partial w}{\partial T}}{1 - \dfrac{x}{w}\dfrac{\partial w}{\partial x} + \dfrac{1}{4}\left( -\dfrac{1}{4} - \dfrac{1}{w} + \dfrac{x^2}{w^2} \right)\left( \dfrac{\partial w}{\partial x} \right)^2 + \dfrac{1}{2}\dfrac{\partial^2 w}{\partial x^2}}$$

This is **Dupire in Gatheral's form**, and it provides a bridge to compute the local vol from the implied vol surface. The numerator is how the total variance changes with respect to time. The denominator appears convoluted but has a simple meaning: it is the ratio of the modeled density from the implied vol surface against the baseline model, the constant vol lognormal, which is the function $g(x)$ from [Appendix B]({{< ref "vol_surface_calibration.md#appendix-b-deriving-the-total-variance-butterfly-condition" >}}) of the implied vol article. When the vol is constant, $\partial w / \partial x$ and $\partial^2 w / \partial x^2$ are both zero, the denominator is $1$, and we are left with $\sigma^2 = \partial w / \partial T$, which is the constant variance itself.

### Taking the Derivatives

The Gatheral formula needs three derivatives of $w$, and the naive way to get them is to bump $x$ and $T$ and take finite differences. It works, but it is not the default way to do it. Differencing amplifies the noise in the implied vol fit, especially for the second derivative in the denominator. The result is a local vol grid that looks jagged even when the input surface is smooth, and the pricing model, Monte Carlo or PDE, inherits that noise.

For $\partial w / \partial x$ and $\partial^2 w / \partial x^2$, analytical derivatives are usually preferred, since they give a closed form on a smooth calibration. The cubic spline from the implied vol article is piecewise polynomial, so its first and second derivatives are available in closed form on each segment and are continuous across the knots by construction. Parametric models such as SSVI are also differentiable in closed form, and parametric surfaces are often the input used in practice.

The maturity derivative, $\partial w / \partial T$, works similarly. The analytical form follows from the total variance interpolation set by the implied vol model. One nuance is worth noting. If total variance is linear in $T$ at fixed $x$, the derivative is piecewise constant and undefined at the quoted expiries themselves, since the intervals on either side give different slopes. On each interval it is well defined, and since variance accumulates as an integral over time, the value assigned at an isolated point does not affect any price. The usual convention is to take the slope of the interval ahead, which leaves the local vol surface right-continuous in maturity and jumping as the simulation crosses each strip.

The analytical approach tells us that in order to use the implied vol surface as an input, instead of just saving the vol at all the grid points, we should store the model itself, spline coefficients or parameters, and take the derivatives from it directly. Local vol calibration still runs downstream of the implied vol fit, but what crosses between them is the model rather than a set of grid points.

In some settings, passing the model along is not possible. Traders may mark the surface directly, sometimes with manual adjustments, or the marking system and the calibration library may be separate enough that only a grid of vols crosses between them. In that case a surface has to be fitted to the grid before anything can be differentiated. The risk here is that the refitted model may differ from the one underlying the original marked vol surface.

### Is Your Input Surface Good Enough?

Not every smooth looking implied vol surface is suitable for local vol calibration. Local variance is a ratio, so it exists and is positive only when the numerator and denominator are, and what makes a surface a good input turns out to be exactly the arbitrage conditions covered in the [implied vol article]({{< ref "vol_surface_calibration.md" >}}).

**The numerator requires no calendar arbitrage.** It is positive when total variance is strictly increasing in maturity at fixed log-moneyness, which is precisely the calendar condition. A violation returns a negative local variance and the calibration fails outright at that point. This is a reasonable requirement for a single underlying such as an equity vol surface. For commodities, where each contract month is a separate underlying, it may not hold strictly, which is the subject of a later section.

**The denominator requires no butterfly arbitrage.** To have a density at all, we need the first and second derivatives of $w$, so the surface must be twice differentiable in $x$. In the implied vol article we used a cubic spline in $\sigma$ against log-moneyness, which satisfies the condition by the chain rule, since squaring a twice differentiable function leaves it twice differentiable. Had we used a simple linear interpolation in $\sigma$ instead, the second derivative would be undefined at the knots and $g$ would not exist. Given that the density exists, requiring it to be positive is exactly the [no butterfly arbitrage condition]({{< ref "vol_surface_calibration.md#negative-density-and-butterfly-arbitrage" >}}).

This is also where the wing treatment matters. Beyond the last anchor there are no quotes to constrain the fit, so the extrapolation we choose is what decides $g$ in the far wings. In the implied vol article we settled on total variance growing linearly in $x$, with a slope consistent with Lee's asymptotic bound, and checked that the resulting density remained positive. The local variance in the wings therefore does not turn negative or blow up.

### Precomputing the Grid

Evaluating Gatheral's formula inside the pricing loop would mean recomputing the local vol at every step. It is more efficient to compute it once on a grid in $(x, T)$ and store it. A PDE solver needs it at every space node on every time step, and those points are fixed in advance, so the grid can be aligned with the solver's mesh and read off directly. Monte Carlo evaluates at whatever points the paths reach, so the grid is sized to cover the region they will visit and values are interpolated between nodes. One of the interpolation choices is bilinear interpolation on local variance: it is cheap enough to call millions of times, and it cannot overshoot, so a grid of positive values stays positive.

How fine the grid needs to be depends on the accuracy of the pricing scheme it feeds. The price carries two errors, one from interpolating the local vol between nodes and one from the pricing scheme's own time and space stepping. If the first is the larger, refining the pricing scheme's time step will not improve the accuracy, since the local vols driving it are already off by more than the stepping error being removed. Whether that is the case can be settled empirically: halve the spacing of the local vol grid, reprice, and see whether the price moves by more than the tolerance set for the model. If it does not, the local vol grid is fine enough and the remaining error belongs to the pricing scheme.

## How Does Local Vol Change the Barrier Option Price?

Now that we know how to calibrate local vols, we can put the model to work and see how local vol and the implied vol surface differ in pricing an exotic, which is the motivation for the entire topic. We start from an arbitrage-free implied vol surface, use it to calibrate the local vols, then price a simple barrier option under both models and compare.

The input surface is hypothetical but carries a realistic skew. It is constructed from the SSVI model. The details of the model are not important for this exercise, only that it has closed forms for $\partial w/\partial x$, $\partial^2 w/\partial x^2$ and $\partial w/\partial T$, and that it is arbitrage-free by construction, which lets us calibrate the local vols without any problems.

The calibrated local vols and the implied vols they came from are shown below, as surfaces and as the one year skew. The surfaces can be rotated for a better view.

{{< iv_lv_comparison >}}

One immediate difference is that the local vol skew is steeper: from $x = 0$ out to $x = +0.4$ the implied vol rises by $6.7$ points while the local vol rises by $13.2$. The slope of the local vol is roughly two times the slope of the implied vol, and that is not a coincidence. The implied vol is the single constant we solve for so that Black-Scholes returns the market price. Since the model holds that vol fixed for the option's whole life, the number we recover stands in for all the volatility the spot experiences along the way, which makes it an average of the local vols on the path it takes. The relationship is only roughly one half for two reasons. The average in the argument is not the arithmetic average of the vols; it is the harmonic average, taken over the levels between the forward and the strike. And the relation itself is exact only in the short maturity limit, whereas the skews compared here are at one year. The derivation of the factor of two, and where the approximations enter, is in [Appendix A](#appendix-a-why-the-local-vol-skew-is-roughly-twice-the-implied-vol-skew).

On the one year skew the local vol sits below the implied vol near the money, $29.8\%$ against $30.7\%$ at the 100 strike. This is driven by the vol smile, under which the market implied distribution has fatter tails than the lognormal. Recall from Gatheral's form that the local variance is the time slope of total variance divided by $g$, the ratio of the market implied density to a lognormal carrying the same total variance at that log-moneyness. That baseline fixes the variance, so the extra tail mass has to be paid for out of elsewhere. The result is fatter tails, thinner shoulders and a higher peak. At the money the market implied density exceeds the lognormal, $g$ is above one, and dividing by it puts the local vol below the implied vol.

**Valuing a Barrier Option**

The contract is an up-and-out call, spot and strike at 100, barrier at 130, one year to maturity, at a risk-free rate of $5\%$. It is an ordinary call that is cancelled if the spot ever trades at 130, so it pays only on paths that finish in the money without having reached the barrier. The payoff is concentrated in the band between the strike and the barrier, and it is discontinuous at the barrier itself, worth nearly 30 just below and nothing at it. Both models are priced with the same basic Monte Carlo so that no part of the gap comes from the numerical scheme: 100,000 paths, 252 time steps, and the same seed.[^mc] Under local vol the grid is precomputed in $(x, t)$ at a spacing of $0.005$ in log-moneyness and $0.001$ years in maturity, and read by bilinear interpolation on local variance. The spacing follows the refinement test set out in the previous section: halving it moves the vanilla price by less than the MC standard error.

The constant vol model has to be given a fair chance. The barrier option is a call struck at 100 made conditional on survival, so the constant vol is set to the implied vol of the 100 strike, $30.70\%$.

Before we compare the barrier prices, it is worth looking at the vanilla first. The two models price it consistently, as they must, and the difference is within the MC standard error. What they do not share is how they get there. The call price is the probability of finishing in the money times the average payoff when it does: local vol finishes in the money more often but pays less on those paths, and the constant vol model does the reverse. The skew is what separates them. Both densities carry the same mean, the forward at $105.13$, since both models are risk neutral, but they distribute the mass differently. The market implied density is heavier in the far left tail, which is what the downside skew prices, and thinner in the far right tail. This is typical of equity markets, where downside protection is in higher demand and is priced richer, which lifts the implied vols on that side. The right panel shows why the call price is unaffected: the density difference weighted by the payoff has two lobes of equal area, so what local vol gains just above the strike it gives back in the far tail, where each path pays far more.

{{< terminal_density >}}

|  | local vol | constant vol | difference |
|---|---|---|---|
| call price | 14.352 | 14.329 | $+0.2\%$ |
| standard error | 0.066 | 0.072 | |
| P($S_T > 100$) | 56.3% | 50.3% | $+6.0$ pts |
| E[payoff \| in the money] | 26.82 | 29.95 | $-3.13$ |

The barrier option prices are significantly different: local vol prices it at $2.348$ against $1.619$ under constant vol, $45\%$ higher.

|  | local vol | constant vol | difference |
|---|---|---|---|
| barrier price | 2.348 | 1.619 | $+45.0\%$ |
| standard error | 0.017 | 0.015 | |
| P(survive) | 63.8% | 62.5% | $+1.3$ pts |
| P(survive and $S_T > 100$) | 22.3% | 16.5% | $+5.8$ pts |
| E[payoff \| survive, in the money] | 11.06 | 10.32 | $+0.74$ |

The probability of survival differs by only $1.3$ points between the two models, so the difference in barrier survival alone cannot account for a price gap of this size. The main driver in this example is the truncation of the terminal payoff. The two lobes that cancelled for the vanilla no longer do: the up-and-out collects only the part between the strike and the barrier, which is the region where local vol has the extra mass, and the far tail that used to offset it is gone. The paying state is $5.8$ points more likely under local vol, and the average payoff in it is higher as well.

## Local Vol Matches the Marginals, Not the Path

From the last section, we see both constant implied vol and local vol repriced the vanilla option, but give different barrier prices. They both integrated their own dynamics honestly to reach the barrier price, so each is self-consistent. The vanilla option quotes cannot tell us which one to believe, because they provide only the distribution of the spot at each maturity, while the barrier depends on the path taken to get there.

What the two models do share is that both make the volatility deterministic, so the vol carries no randomness of its own, but the market does not behave that way. The VIX, a measure of SPX implied vol, may still move when the SPX does not, and there is a market in options on the VIX. Those option contracts would be unnecessary if vol were deterministic given the index level. The return and the vol also tend to be negatively correlated in many markets: when the price drops, vol often rises. These market behaviors indicate vol is not deterministic, and should be treated as a second source of randomness correlated to the returns.

Before discrediting the local vol model immediately due to its lack of a second driver, let's ask what simplifying assumption it is making. Suppose the true underlying dynamics is that the spot follows a diffusion whose instantaneous volatility is some process $\beta_t$,

$$dS_t = r S_t \, dt + \beta_t \, S_t \, dW_t$$

$\beta_t$ is unknown to us. It may be driven by its own Brownian motion, or by several other factors. It is not determined by $S_t$ the way the local vol function is. Take $f$ to be any smooth function of the spot. Ito's lemma gives

$$df(S_t) = \left( r S_t f'(S_t) + \frac{1}{2} \beta_t^2 S_t^2 f''(S_t) \right) dt + \beta_t S_t f'(S_t) \, dW_t$$

and taking expectations kills the stochastic integral, leaving

$$\frac{d}{dt}\mathbb{E}\left[ f(S_t) \right] = \mathbb{E}\left[ r S_t f'(S_t) \right] + \frac{1}{2} \mathbb{E}\left[ \beta_t^2 S_t^2 f''(S_t) \right]$$

The first term involves only $S_t$. The second couples $\beta_t$ to $S_t$, and this is where the unknown process has to be dealt with. Applying the tower property, whatever $\beta_t$ was doing, only its average at each level survives once it is multiplied by a function of the spot alone, so

$$\frac{d}{dt}\mathbb{E}\left[ f(S_t) \right] = \mathbb{E}\left[ r S_t f'(S_t) + \frac{1}{2} \mathbb{E}\left[ \beta_t^2 \mid S_t \right] S_t^2 f''(S_t) \right]$$

Now run the same expectation in the local vol model, where the volatility is the deterministic function $\sigma(S_t, t)$. Ito gives the same expression with $\sigma^2(S_t, t)$ in place of $\beta_t^2$, and no conditioning is needed since it is already a function of the level:

$$\frac{d}{dt}\mathbb{E}\left[ f(S_t) \right] = \mathbb{E}\left[ r S_t f'(S_t) + \frac{1}{2} \sigma^2(S_t, t) S_t^2 f''(S_t) \right]$$

We want the local vol model to return the same vanilla prices as the true dynamics, so we ask what $\sigma$ makes the two evolutions agree. Both start from the same spot today, so matching the evolution keeps the marginals together as $t$ advances, and agreeing on $\mathbb{E}[f(S_t)]$ for every smooth $f$ is agreeing on the marginal density itself. A call price is its discounted expected payoff,

$$C(K, T) = e^{-rT} \int \max(S - K, 0) \, p(S, T) \, dS$$

and by Breeden-Litzenberger the density is recovered from those prices by differentiating twice in $K$. Agreeing on the marginal density is therefore agreeing on the vanilla prices at every strike, and the two densities coincide if

$$\sigma^2(S, t) = \mathbb{E}\left[ \beta_t^2 \mid S_t = S \right]$$

This is the **Gyongy** result. No matter what the true volatility process is, there is a local vol model that reproduces its marginal density at every maturity, and the local variance it needs is the conditional average of the true instantaneous variance at that level and time. In the barrier example, at $S = 130$ one year out the calibrated local vol is $26.1\%$, or $0.068$ in local variance. This value is sometimes misinterpreted as a forecast: if the spot reaches 130, this is the variance it will have there. Gyongy gives a clear interpretation. It is the conditional average of the true instantaneous variance across all volatility states that can occur when the spot is at 130 at that time.

Although the local vol model matches the marginal distribution by keeping the mean of the instantaneous variance, it throws away the distribution of it at a given level and time. A single-asset barrier is path dependent in a relatively simple way: its value depends on whether a particular level is breached, and local vol still remains widely used for such products in practice. For a forward starting option, whose value depends on the smile seen from a future date rather than from today, local vol is inadequate, since it cannot model the forward vol skew. Fixing that requires giving the vol randomness of its own. This leads to the stochastic vol model, which we will cover in the next article.

## Adapting Local Vol to Commodities

Everything so far has assumed an equity setting, where the vol surface holds the vols of a single underlying seen at a range of maturities. Commodity options work differently. They are written on futures, and each contract month is its own underlying with its own futures price, so moving along the maturity axis changes the underlying as well as the horizon. The December contract and the following June contract respond to different supply and demand, and an option on one is not a direct statement about the other.

That can be a problem for Dupire, since the numerator of the local variance is $\partial w / \partial T$. Taking that derivative on an equity surface compares the same underlying at two horizons. Taking it on a commodity surface compares two different underlyings, and the difference between them reflects not only the passage of time but also the distinction between the contracts themselves.

Whether this matters in practice depends on how close the contract months are to being the same underlying. When they are highly correlated, treating the strip as a single underlying is a reasonable approximation. Crude is workable on this basis, and precious metals more so, since a gold futures curve is close to deterministic carry and the contract months move almost together. For natural gas the assumption may not be appropriate. A January contract and a July contract are vastly different and not interchangeable: they are exposed to different weather, different storage positions, and different demand. Their prices are only loosely correlated, and their implied vol levels also differ. Differencing across them does not isolate the effect of time. Across a gas strip total variance is also not guaranteed to increase along the surface, and where the calendar condition breaks the formula returns a negative local variance. Local vols fail to calibrate there.

Some markets do list several option expiries against a single contract, such as the Henry Hub weeklies, which are listed for four weeks. But that reaches only the nearest contract, and the rest of the strip is left with one expiry each. For most of the curve the maturity dimension therefore has to be modeled rather than read off the quotes. One idea is to stop treating each contract month as its own underlying and to pool them into a single background process, so that the one expiry each contract carries becomes one point on a maturity axis shared by all of them.

Nastasi, Pallavicini and Sartorelli take this route with what they call a fictitious spot price, a process $s_t$ that is not traded and need not match any quoted spot. Futures prices are derived from it as

$$F_t(T) = F_0(T)\left( 1 - (1 - s_t)\, e^{-\int_t^T a(u) \, du} \right)$$

which recovers today's curve exactly and can be inverted for the level of $s_t$ at which the contract trades at any given strike. An option on a contract month expiring in six months therefore becomes an option on $s_t$ six months out, struck at the corresponding level. Each contract contributes its single expiry at its own maturity, and the strip becomes a set of maturities on one process. Dupire is then extended to the drift this dynamics carries and calibrated against the pooled set.

In this model the spot is assumed to mean revert at speed $a$, and it cannot be pinned down from vanilla quotes and needs to be calibrated from calendar spread options. I have not implemented this model and have no desk experience with it, so I am pointing at the reference rather than vouching for how it behaves in production. It would however be an interesting topic to dive into, so I may write a follow up article on it.

## Appendix A: Why the Local Vol Skew is Roughly Twice the Implied Vol Skew

To see where the factor of two comes from, let's start by examining how the option price behaves in the short maturity limit, where $T \to 0$.

The price of an option is set by how far the strike sits from the forward measured in standard deviations,

$$d = \frac{x}{\sigma\sqrt{T}}$$

The forward carries no drift, so $F_T = F e^{-\sigma^2 T/2 + \sigma\sqrt{T} Z}$ with $Z$ standard normal, and the call can be approximated as

$$C \approx e^{-rT} F \sigma \sqrt{T} \; \mathbb{E}\left[ (Z + d)^+ \right]$$

Everything specific to the contract has been pushed into the two factors. The scale $e^{-rT} F \sigma \sqrt{T}$ is the discounted size of a one standard deviation move on the forward, which is what an option is worth per unit of standardized payoff, and the expectation is a pure function of $d$, the number of standard deviations the strike sits away. Writing it out,

$$\mathbb{E}\left[ (Z + d)^+ \right] = \int_{-d}^{\infty} (z + d)\,\varphi(z)\,dz$$

with $\varphi$ the normal density. At the money this is trivial, giving $\mathbb{E}[Z^+] = 1/\sqrt{2\pi}$, so

$$C \approx \frac{1}{\sqrt{2\pi}} \, e^{-rT} F \sigma \sqrt{T}$$

The price is proportional to the vol, and matching prices means matching vols. This shows $\sigma_{IV}(0) = \sigma_{LV}(0)$: the two must agree at the money as $T \to 0$. That tells us the vol levels, but nothing about how the vol skews of the two models relate. So we need to look at a non-trivial case, one involving strikes sitting away from the forward. Take an out-of-the-money call, where $d < 0$. The in-the-money case follows the same argument due to put-call parity. Changing variables to $z = -d + u$,

$$\int_{-d}^{\infty} (z + d)\,\varphi(z)\,dz \approx \varphi(d) \int_0^\infty u\, e^{-|d|u}\,du = \frac{\varphi(d)}{d^2}$$

so

$$C \approx e^{-rT} F \sigma \sqrt{T} \, \frac{\varphi(d)}{d^2} \qquad \Longrightarrow \qquad \log C \approx -\frac{d^2}{2} + O(\log|d|)$$

As $T \to 0$ the distance $|d|$ grows to infinity, so the price is dominated by the Gaussian term while the prefactor contributes only at logarithmic order. So two models agreeing on the price must agree on $d$. Under Black-Scholes the volatility is constant for a given strike, so the distance is the log-moneyness divided by the standard deviation over the option's life:

$$d = \frac{x}{\sigma_{IV}(x)\sqrt{T}}$$

Under local vol the volatility depends on the level, so the distance has to be accumulated level by level. Each stretch $du$ of log-moneyness between the forward and the strike is divided by the volatility available at that level, so the local vol analogue of $d$ is

$$d = \frac{1}{\sqrt{T}} \int_0^x \frac{du}{\sigma_{LV}(u)}$$

Equating the two, the $\sqrt{T}$ cancels:

$$\frac{x}{\sigma_{IV}(x)} = \int_0^x \frac{du}{\sigma_{LV}(u)}$$

This is the **Berestycki-Busca-Florent** result, exact in the limit $T \to 0$. It reveals that the implied vol is the harmonic average of the local vol over the levels between the forward and the strike. The factor of two follows immediately. Take the local vol to be linear near the money, $\sigma_{LV}(u) = \sigma_0 + bu$, and expand the reciprocal for small $bu / \sigma_0$,

$$\frac{1}{\sigma_{LV}(u)} \approx \frac{1}{\sigma_0}\left( 1 - \frac{bu}{\sigma_0} \right)$$

Averaging over $[0, x]$ gives

$$\frac{1}{\sigma_{IV}(x)} \approx \frac{1}{\sigma_0}\left( 1 - \frac{bx}{2\sigma_0} \right) \qquad \Longrightarrow \qquad \sigma_{IV}(x) \approx \sigma_0 + \frac{bx}{2}$$

The implied vol slope is $b/2$ against the local vol's $b$. Two approximations are in play. The relation holds in the short-maturity limit, and the local vol is treated as linear over the range traveled. Neither holds cleanly at the one year maturity discussed in the main text, where the local vol is convex rather than linear. The measured ratio comes out slightly below two, and the instantaneous slope at the money gives $1.8$. The factor of two is a rule of thumb, not an exact identity.

[^mc]: No conditional Monte Carlo or Brownian bridge correction is applied for barrier crossings between monitoring dates, so both prices carry the same discretization bias and the comparison stays like for like.