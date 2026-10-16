---
title: "Constructing the Implied Volatility Surface: From Market Quotes to an Arbitrage-Free Fit"
date: 2026-06-26
draft: false
math: true
tags: ["volatility-surface", "implied-volatility", "black-model", "no-arbitrage", "butterfly-arbitrage", "calendar-arbitrage", "commodities", "crude-oil"]
ShowToc: true
TocOpen: false
---

## Why This Matters

For vanilla options, the simple models are usually sufficient. A plain call or put can be priced off Black-Scholes directly; we often do not need to reach for local volatility or a stochastic volatility model. Those heavier models earn their place with exotics, where the payoff depends on how the smile behaves rather than just its level today. For a vanilla option, we take the market's implied volatility at the relevant strike and maturity and feed it into Black-Scholes. That assumes a volatility surface already exists. Before Black-Scholes can price anything, the surface it reads from has to be built, and building it can be less straightforward than it appears.

Fitting a smile to the quotes is the routine part, but a surface that looks smooth can still imply an arbitrage, so we need to check the fit rather than trust it unconditionally. The deep wings also need proper treatment, since that is where the market quotes thin out or vanish. This article works through the vol surface construction end to end, from the raw option quotes to an arbitrage-free surface.

## Crude Oil as a Working Example

It helps to fix a concrete product and carry it through every step. WTI, traded on CME, is the most liquid crude benchmark, and its options surface is richly quoted. Brent trades on ICE with good liquidity in the near-dated contracts, thinning in the back months and deep wings. There is more than one way to build its surface. One approach calibrates directly from Brent's own quotes, which reflects Brent's own skew and supply-demand. Another builds Brent relative to WTI through a Brent-WTI volatility spread. The two benchmarks are highly correlated, so a spread can be more stable than an independent fit, and it leans on WTI's deeper quotes where Brent runs sparse. The cost is that the spread itself has to be modeled, and its assumptions carry their own risk. Which approach is more suitable depends on the book, and on whether the risk reaches into the regions where Brent's quotes run out.

For this article I calibrate the Brent surface directly and stay in the liquid near-dated region. The instrument is the [ICE Brent Crude American-style Options](https://www.ice.com/products/218/Brent-Crude-American-style-Options), and two features of its contract make the pricing clean enough to keep the focus on surface construction.

First, the option is written on a future, so the future serves directly as the forward the Black model takes as its input. Second, the option is futures-style margined. As shown in a [previous article]({{< relref "future_style_margining_options.md" >}}), the American option carries no early-exercise premium, so it can be priced as European. Together these collapse the pricing to the Black model, Black-Scholes written on the future in place of spot, with no discount factor on either term:

$$C = F\,N(d_1) - K\,N(d_2), \qquad P = K\,N(-d_2) - F\,N(-d_1), \qquad d_{1,2} = \frac{\ln(F/K) \pm \tfrac12 \sigma^2 \tau}{\sigma \sqrt{\tau}}.$$

Crude has fat tails on both sides, so the Brent smile rises in both wings, and its skew tilts with the market's fears: toward the calls when supply is threatened, as after Russia's invasion of Ukraine in 2022, and toward the puts when demand collapses, as during COVID in 2020.

## Building the Volatility Surface

### 1. Choosing the Anchor Grid

A surface is a function of two variables: across strike, where along the smile each volatility sits, and across maturity, which expiries the surface stores.

**The strike axis.** Where along the smile should the surface store its volatilities? Fixed strikes are one option, but as the future moves a fixed strike drifts to a different moneyness, so a common choice for Brent is to anchor at fixed deltas instead, a stable and portable coordinate. In practice the smile only stays roughly fixed in delta, not perfectly, but as a coordinate to store and interpolate volatilities it does the job.

With delta chosen as the coordinate, the next question is which deltas, and how many. For Brent there is no standard set, so the grid is left to the desk to decide. A denser grid captures more skew and curvature, but only where a quote or a wing extrapolation supports each point. Spacing the anchors too finely also magnifies the noise in their values into curvature, so the fitted smile turns wiggly and the greeks jump between recalibrations. A sparser grid is more stable but can miss curvature between its points. Where the grid lands is largely a modeling decision, though practical factors such as downstream risk or pricing system setup can also shape the choice.

For this article I will use a representative grid: closely spaced near the money and extending into each wing. The anchors are quoted as out-of-the-money options, using calls for strikes above the forward and puts for strikes below it:

$$10\Delta C,\; 15\Delta C,\; 25\Delta C,\; 35\Delta C,\; 50\Delta,\; 35\Delta P,\; 25\Delta P,\; 15\Delta P,\; 10\Delta P.$$

Each anchor can also be written as a put delta, which sweeps the whole smile monotonically:

$$-0.90,\; -0.85,\; -0.75,\; -0.65,\; -0.50,\; -0.35,\; -0.25,\; -0.15,\; -0.10.$$

**The maturity axis.** Which expiries should the surface store? Here the product decides for us. Brent options are listed on monthly futures, each expiring against its delivery month's future, so the natural maturity anchor is the contract month. For this article we take the first several contract months. In practice, we need to store enough contract months to cover a trading book's risk.

### 2. Selecting the Quotes

A volatility surface is only as good as the data that feeds it, and not every market quote should be used automatically. A quote earns its place only if it carries clean information about where the option trades. The criteria below are a reasonable starting set rather than an exhaustive list:

- **Use the out-of-the-money strike.** At any strike there is both a put and a call. The OTM one is the better quote to use: it is the more actively traded and tightly quoted of the two. At the strike equal to the future the two coincide: put-call parity makes the call and put equal there, so either serves, or their implied vols can be averaged.
- **Drop quotes below intrinsic.** No option can trade below its intrinsic value without implying an arbitrage.
- **Check monotonicity in strike.** Call prices must not rise as the strike rises, and puts must not fall. A strict reversal is an arbitrage, and the offending quotes should be dropped.
- **Drop quotes at the minimum tick.** A deep out-of-the-money option can be worth less than a single tick, so its price is usually floored at the minimum increment and carries almost no volatility information. The vol backed out of it is unreliable.
- **Filter on liquidity, if the feed carries it.** Open interest and volume indicate whether a strike is genuinely active; a strike with negligible activity can be dropped as stale.

### 3. Constructing Each Strip

We build one maturity strip at a time. The market gives prices at strikes, but the strip stores a volatility at each delta anchor, so the construction has to bridge the two. Once the strip is built, we can interpolate a volatility at any strike, not just the quoted ones. There are different ways to interpolate, and the method we use here is a cubic spline in log-moneyness, since volatility is measured on log returns.

**Invert each quote.** The Black model cannot be solved for $\sigma$ in closed form, so each quote, taken as the mid of bid and ask, is inverted numerically. This is exactly the single-option problem solved with [Newton's method]({{< relref "newton_vs_brent_vol_solver.md" >}}) in a previous article. Each quote goes in, its implied volatility comes out.

**Order into three vectors.** Suppose $n$ quotes survive selection for this maturity strip. With a volatility for each, compute two more quantities per strike: the log-moneyness $x = \ln(F/K)$, and the put delta

$$\Delta = N(d_1) - 1, \qquad d_1 = \frac{x + \tfrac12 \sigma^2 \tau}{\sigma \sqrt{\tau}}.$$

Sort all three by log-moneyness, giving three aligned vectors that describe the raw market smile:

$$\mathbf{x} = [x_1, x_2, \ldots, x_n], \qquad \boldsymbol{\sigma} = [\sigma_1, \sigma_2, \ldots, \sigma_n], \qquad \boldsymbol{\Delta} = [\Delta_1, \Delta_2, \ldots, \Delta_n].$$

The $i$-th entries $(x_i, \sigma_i, \Delta_i)$ belong to the same quote: its log-moneyness, its implied volatility, and the delta implied by that volatility and strike, with $x_1 < x_2 < \cdots < x_n$.

The anchor volatilities are what the strip stores, one at each delta anchor. Before solving for them, two things about the data need checking.

> **Do the quotes cover the delta anchors?** An anchor's volatility can be interpolated only if it sits within the deltas the quotes span. Suppose we want the $-0.90$ anchor but the quotes reach only $-0.85$. There is nothing to interpolate from, so the anchor has to be extrapolated past the data, which is the wing problem the later section handles. For Brent, front-month options are deeply liquid, so coverage there is rarely a problem; it is the back of the curve, where trading thins, that is more likely to leave the deep wings uncovered.
>
> **Is delta monotone in strike?** Matching a delta anchor to a strike only makes sense if each delta corresponds to a single strike. If the same delta were reached at two different log-moneyness values, that delta would have no unique implied volatility. On a well-behaved smile this holds automatically. To see when it can fail, differentiate $\Delta = N(d_1) - 1$ along the smile, with $\sigma$ depending on $x$:
>
> $$\frac{d\Delta}{dx} = n(d_1)\left(\frac{\partial d_1}{\partial x} + \frac{\partial d_1}{\partial \sigma}\,\sigma'(x)\right) = n(d_1)\left(\frac{1}{\sigma\sqrt{\tau}} - \frac{d_2}{\sigma}\,\sigma'(x)\right).$$
>
> The first term, $\partial d_1/\partial x$, is the moneyness effect: how much $d_1$ shifts for a step $dx$ in log-moneyness with the volatility held fixed. The second term, $\frac{\partial d_1}{\partial \sigma}\,\sigma'(x)$, is the skew effect: how much $d_1$ shifts because the volatility itself changes along the smile. Since $n(d_1) > 0$, delta is monotone as long as the skew effect stays below the moneyness effect, which rearranges to $\sigma'(x)\sqrt{\tau}\,d_2 < 1$. Near the money $d_2$ is close to zero and the condition holds easily; it is in a steep wing, or at a long maturity, that the skew effect can win and reverse delta. So the monotonicity is worth examining before proceeding.

**Solve the volatility at each anchor delta.** Take a target delta, say $-0.25$. The market rarely lists a strike whose delta is exactly $-0.25$, so we have to locate that strike ourselves, and here is the difficulty: to locate it we need the volatility there, but to read the volatility off the smile we need the strike. Each depends on the other. One way out of the circularity is a fixed-point iteration: guess a volatility, derive the strike it implies through the link between delta and log-moneyness, interpolate the smile based on the strike to get a new volatility, and repeat until the guess stops moving. Below are the steps for one iteration.

a. **Warm start.** Take the volatility of the market quote whose delta is nearest the target as the initial guess $\sigma^{(0)}$. Interpolating the volatilities of the two quotes that bracket the target delta gives a slightly better start, though the iteration converges quickly enough that the choice matters little.

b. **Solve for the strike ratio.** With the target delta and the current volatility $\sigma^{(k)}$ fixed, the delta equation inverts in closed form. With $d_1^\ast = N^{-1}(\Delta_{\text{target}} + 1)$,

$$x^{(k)} = d_1^\ast\,\sigma^{(k)}\sqrt{\tau} - \tfrac12 \big(\sigma^{(k)}\big)^2 \tau.$$

c. **Interpolate the smile.** Using the cubic spline $\mathrm{cs}(\cdot)$ fitted to the vectors $\mathbf{x}$ and $\boldsymbol{\sigma}$, interpolate the volatility at the new log-moneyness $x^{(k)}$ to get the next guess, $\sigma^{(k+1)} = \mathrm{cs}\big(x^{(k)}\big)$.

d. **Check consistency.** Recompute the delta from $x^{(k)}$ and $\sigma^{(k+1)}$. If the volatility has stopped moving and the delta is within tolerance of the target, stop; the anchor volatility is $\sigma^{(k+1)}$. Otherwise take $\sigma^{(k+1)}$ as the new guess and return to step b.

The iteration converges quickly in the covered range, where the quotes are dense and the smile well behaved. In code:

```python
from scipy.interpolate import CubicSpline
from scipy.stats import norm
import numpy as np

def solve_anchor_vol(target_delta, x, sigma, F, tau,
                     tol_vol=1e-6, tol_delta=1e-6, max_iter=200):
    smile = CubicSpline(x, sigma)            # vol as a function of x = ln(F/K)
    d1_star = norm.ppf(target_delta + 1.0)   # put delta = N(d1) - 1

    # warm start: vol of the quote whose delta is nearest the target
    delta_mkt = norm.cdf((x + 0.5*sigma**2*tau) / (sigma*np.sqrt(tau))) - 1.0
    s = sigma[np.argmin(np.abs(delta_mkt - target_delta))]

    for _ in range(max_iter):
        x_n = d1_star*s*np.sqrt(tau) - 0.5*s**2*tau   # strike ratio at fixed delta, vol
        if not (x[0] <= x_n <= x[-1]):                # outside the quotes: a wing problem
            raise ValueError("anchor lies beyond the quoted range; extrapolate the wing")
        s_new = float(smile(x_n))                     # vol the smile assigns there
        d1 = (x_n + 0.5*s_new**2*tau) / (s_new*np.sqrt(tau))
        delta = norm.cdf(d1) - 1.0
        if abs(s_new - s) < tol_vol and abs(delta - target_delta) < tol_delta:
            return s_new
        s = s_new
    raise RuntimeError("anchor solve did not converge")  # fall back to bracketed solve
```

This iteration is fast but not unconditionally convergent: it settles when the smile is gently sloped at the anchor, which usually holds across the covered range in practice, and it can slow where the wing steepens. When the iteration stalls or converges slowly, the fallback is a bracketed root-find on $x$: find two market points whose deltas straddle the target and solve with Brent's method.

<details style="margin-bottom: 1.5em;">
<summary>When does the anchor solve converge?</summary>

The anchor solve iterates on the volatility: guess a value, derive the strike it implies, read the smile there, repeat. Written as a map, $\sigma^{(k+1)} = \Phi(\sigma^{(k)})$ with $\Phi(\sigma) = \mathrm{cs}\big(x(\sigma)\big)$, where $\mathrm{cs}(\cdot)$ is the cubic spline through the market smile and $x(\sigma)$ is the closed-form strike solve at the fixed target delta.

Let $\sigma^\ast$ be the answer, the fixed point $\sigma^\ast = \Phi(\sigma^\ast)$, and $e_k = \sigma^{(k)} - \sigma^\ast$ the error. Subtracting the fixed-point equation from the iteration gives $e_{k+1} = \Phi(\sigma^{(k)}) - \Phi(\sigma^\ast)$, and expanding $\Phi$ to first order around $\sigma^\ast$ leaves

$$e_{k+1} \approx \Phi'(\sigma^\ast)\,e_k.$$

The error is multiplied by $\Phi'(\sigma^\ast)$ each step, so the iteration converges when $|\Phi'(\sigma^\ast)| < 1$, and faster the smaller that number.

What is $\Phi'$? It is a composition, so the chain rule gives $\Phi'(\sigma) = \mathrm{cs}'\big(x(\sigma)\big)\,x'(\sigma)$. The strike solve is $x(\sigma) = d_1^\ast\,\sigma\sqrt{\tau} - \tfrac12\sigma^2\tau$ with $d_1^\ast = N^{-1}(\Delta_{\text{target}}+1)$ constant, so $x'(\sigma) = d_1^\ast\sqrt{\tau} - \sigma\tau$. At the fixed point,

$$\Phi'(\sigma^\ast) = \mathrm{cs}'(x^\ast)\,\big(d_1^\ast\sqrt{\tau} - \sigma^\ast\tau\big) = \mathrm{cs}'(x^\ast)\,\sqrt{\tau}\,d_2^\ast.$$

Near the money $d_2^\ast$ is close to zero, which keeps $|\Phi'|$ small, so the iteration settles in a few steps even under a steep skew. In the wings both $|d_2^\ast|$ and the slope grow, so the iteration slows, and once $|\Phi'|$ reaches one it no longer converges. Note that this convergence condition is closely linked to the delta monotone condition discussed earlier, $\sigma'(x)\sqrt{\tau}\,d_2 < 1$: both depend on the same quantity, since $\sigma'$ is exactly $\mathrm{cs}'$.

</details>

### 4. Assembling the Surface

Repeating step 3 for each stored maturity gives one strip per contract month, and together these strips form the surface.

Reading a volatility off the surface takes two lookups:

- **Within a strip (one maturity).** A strike's volatility comes from translating it to log-moneyness $x = \ln(F/K)$, then evaluating that strip's cubic spline at that $x$.
- **Across strips.** Each strip belongs to a specific contract month's future, so what matters is which future an option settles against, not where its expiry sits in calendar time. An option on the September future reads its volatility off the September strip. That treats August and September as distinct underlyings rather than blending two distributions by calendar distance, which is common practice in commodities. The contrast is a single underlying seen at different horizons, like a stock or an equity index, where the strips share one distribution stretched over time, so a real maturity axis exists to interpolate across. There a common choice is to interpolate total variance linearly in time at fixed moneyness, since it is total variance, not volatility, that drives the dispersion of the price at each maturity.

## Spot the Hidden Issues

Suppose after going through construction steps 1 through 4, we end up with a surface that looks like this.

{{< vol_surface_spline >}}

It looks fine. Smooth strips, no kinks, every quote reproduced. But reproducing the quotes is not the same as being arbitrage-free, and the question we actually need to answer is whether the prices this surface implies are prices a rational market could support.

### Negative Density and Butterfly Arbitrage

The most basic property any set of option prices must satisfy is that the risk-neutral density it implies stays non-negative everywhere. So before trusting a strip, recover its implied density and check the sign.

But how do we translate the volatility surface into a density? Not through the lognormal density: the smile exists precisely because the distribution is not lognormal, with fatter tails on both sides. There is instead a model-free way, reading the density straight out of the call prices. It is the second derivative of the call price in strike,

$$q(K) = \frac{\partial^2 C}{\partial K^2}.$$

To see why, take the first derivative first. Move the strike from $K$ to $K + dK$. In every scenario where the option finishes out of the money it still pays zero, so those scenarios are unchanged; in every scenario where it finishes in the money it now pays $dK$ less. The drop in the call value is $dK$ weighted by the probability of landing in the money, $dC = -\,\mathbb{P}(F > K)\,dK$, so

$$\frac{\partial C}{\partial K} = -\,\mathbb{P}(F > K).$$

Now differentiate again. Moving $K$ to $K + dK$ changes the in-the-money probability by the mass in that band, and that mass over $dK$ is the density at $K$:[^density_df]

$$\frac{\partial^2 C}{\partial K^2} = -\,\frac{\partial}{\partial K}\mathbb{P}(F > K) = \frac{\mathbb{P}(K < F < K + dK)}{dK} = q(K).$$

This identity is the <strong id="breeden-litzenberger">Breeden-Litzenberger</strong> formula (1978). In practice this second derivative can be approximated numerically, on a grid of strikes spaced $\Delta K$ apart, by the central difference

$$q(K) \approx \frac{C(K - \Delta K) - 2\,C(K) + C(K + \Delta K)}{(\Delta K)^2}.$$

The numerator is the cost of a butterfly: buy one call at $K - \Delta K$, buy one at $K + \Delta K$, sell two at $K$. Its payoff is a tent peaking at $K$ and zero beyond the wings, never negative, so its cost can never be negative either, which is why the negative-density violation is called a butterfly arbitrage.

{{< butterfly_payoff >}}

The density check runs strip by strip. We still use the Black model to get the call prices, but here it imposes no distribution: fed a different vol at each strike, it only translates each quote from vol into price, and the shape of the density is set entirely by the smile.

```python
import numpy as np
from scipy.stats import norm

def implied_density(strike_grid, vol_of_strike, F, tau):
    """Risk-neutral density q(K) = d2C/dK2 from a fitted smile, on a futures-style option."""
    K = strike_grid
    sigma = vol_of_strike(K)                       # volatility the strip assigns each strike
    d1 = (np.log(F/K) + 0.5*sigma**2*tau) / (sigma*np.sqrt(tau))
    d2 = d1 - sigma*np.sqrt(tau)
    C = F*norm.cdf(d1) - K*norm.cdf(d2)            # undiscounted Black call price
    dCdK = np.gradient(C, K)                        # -P(F > K)
    return np.gradient(dCdK, K)                     # d2C/dK2 = q(K)
```

Run this over each strip in the surface and every strip but Dec 2026 comes back clean, a density that stays non-negative across every strike. Dec 2026 does not.

{{< vol_density >}}

Its density crosses below zero in the deep downside wing, near a strike of 53, bottoming around $-0.0014$. On the vol plot this is the $-0.10$-delta anchor, the far end of the Dec strip where the smile bends up into its wing. The bend is a touch too sharp there, and past a point extra curvature in the smile turns the call price concave in strike, which is a negative density. The negative dip is likely too small and too deep in the wing to trade against once bid-ask and fees are counted, so it is not necessarily free money. It is still a defect to fix: the surface is not arbitrage-free as constructed. The first step is to recheck the quotes that shape the smile near that anchor and drop or adjust any that are stale or noisy; if the quotes hold up, the problem lies in the interpolation rather than the data, since the cubic spline is forced to bend too sharply through them. So instead of relying on the spline from the existing quotes, we can leave the $-0.10$ delta anchor to the wing treatment discussed later.

The same strip also shows a choppy density closer to the money, a peak near 62, a dip near 66, another peak near 71. These sit around the $-0.25$-delta anchor, where the Dec smile bends up most sharply. The smile there looks smooth, but the density reads its curvature, not its level, and a cubic spline's curvature is only piecewise linear, with a kink at every node, even where the curve itself looks smooth. Every value stays positive, so no butterfly is mispriced, but the choppiness makes the Greeks and hedge ratios unstable. Desks that care more about exact repricing than smooth Greeks can tolerate this.

### Calendar Arbitrage Across Strips

The butterfly check validates each strip on its own. The other classic check looks across strips: calendar arbitrage. At a fixed moneyness, total variance $w = \sigma^2\tau$ must not decrease with maturity, because a longer-dated option is the same price process given more time to diffuse and cannot carry less variance than a shorter one. If it did, a calendar spread across the two would lock in a riskless profit. The comparison has to be at the same moneyness on each strip, not the same delta, since a fixed delta maps to a different strike at each maturity.

This assumes one underlying seen at two horizons, which holds for a stock or an equity index but not strictly for Brent, where each strip is its own future and the Sep and Oct futures are distinct underlyings. Correlation is what makes the check still worth running: front-month Brent futures move nearly one-for-one, so across the near months they behave close enough to a single underlying that the condition applies. For commodities it is a heuristic, not a strict law, and it can fail on volatility seasonality: a January natural gas contract expiring into peak winter demand carries higher volatility than a following April contract, so total variance falls from the winter month to the spring one with no arbitrage, since the two are different bets on different weather regimes.

For this surface the strips are the near months, where the check applies. Comparing total variance at fixed moneyness across maturities, it rises monotonically from Sep 2026 through Jan 2027.

{{< calendar_total_variance >}}

[^density_df]: These clean forms rely on the futures-style margining set up earlier: the quote is the undiscounted option value, so $C(K) = \mathbb{E}[(F-K)^+]$ with no discount factor. For a premium-paid option, $C(K) = e^{-r\tau}\,\mathbb{E}[(F-K)^+]$, and the discount factor carries through both derivatives, $\partial C/\partial K = -e^{-r\tau}\,\mathbb{P}(F > K)$ and $\partial^2 C/\partial K^2 = e^{-r\tau}\,q(K)$. Recovering a unit-mass density then means dividing by $e^{-r\tau}$.

## What Happens to the Wings

So far we have built the surface assuming the market provides liquid quotes covering the delta anchors. That is not always the case. A position may need vol at a strike so far out of the money the exchange lists no option there at all, and the deep quotes that do exist are often stale, so they fail the selection in the first place. Either way the strip has to say something about a region the data does not cover, and the cubic spline is the wrong tool for it: a spline interpolates, and past the last anchor it simply continues its end cubic, a curve with no reason to resemble a volatility wing.

So how do we extrapolate the wings? What is a reasonable shape, one that at least keeps the density non-negative, the butterfly condition we checked earlier?

### How Variance May Grow in the Wing

Rather than study volatility against log-moneyness $x = \ln(F/K)$, work in total variance $w(x) = \sigma^2(x)\,\tau$, since variance is what drives the dispersion of returns. The question is how fast $w$ may grow as $x$ runs out into the wing.

The butterfly condition from the density section can be written entirely in these coordinates. Non-negative density is equivalent to the function $g(x)$ below, built from $x$, $w$ and the first two derivatives of $w$, staying non-negative at every $x$,

$$g(x) = \left(1 - \frac{x\,w'}{2w}\right)^2 - \frac{w'^2}{4}\left(\frac{1}{w} + \frac14\right) + \frac{w''}{2} \;\ge\; 0,$$

where $w' = dw/dx$ and $w'' = d^2w/dx^2$. This form follows Gatheral and Jacquier (2014); a sketch of its derivation is in [Appendix A](#appendix-a-sketch-of-the-total-variance-butterfly-condition). To test how fast $w$ may grow in the wing, set $w = Cx^\alpha$ with $\alpha \ge 0$ and read the sign of $g$ as $x \to \infty$, the put wing; the call wing, $x \to -\infty$, follows a similar analysis after the sign change $u = -x$. A negative $\alpha$ would make total variance fall toward zero in the wing, while the Brent smile rises in both wings.

**Faster than linear, $\alpha > 1$.** With $w = Cx^\alpha$, $w' = C\alpha x^{\alpha-1}$, $w'' = C\alpha(\alpha-1)x^{\alpha-2}$, the three terms of $g$ scale as

$$\left(1 - \frac{x w'}{2w}\right)^2 = \left(1 - \frac{\alpha}{2}\right)^2, \qquad -\frac{w'^2}{4w} + \frac{w''}{2} = \frac{C\alpha(\alpha-2)}{4}\,x^{\alpha-2}, \qquad -\frac{w'^2}{16} = -\frac{C^2\alpha^2}{16}\,x^{2\alpha-2}.$$

The first is constant. The last dominates the middle one, since $2\alpha - 2 > \alpha - 2$, and its exponent is positive, so it diverges: $g \to -\infty$. Superlinear growth leads to butterfly arbitrages in the wing.

**Linear, $\alpha = 1$.** For $w = Cx$, $g \to \tfrac14 - C^2/16$ as $x \to \infty$, non-negative only when $C \le 2$. A linear wing is arbitrage-free up to that slope; steeper reintroduces the violation. This is the bound from Lee's (2004) moment formula, and the construction has to respect it. Note this is an asymptotic property: a slope within the bound keeps the density non-negative far out in the wing, but does not by itself guarantee it at every $x$.

**Slower than linear, $0 \le \alpha < 1$.** Here $g$ stays positive, so this wing is also arbitrage-free. But it assigns the deep wing a lower vol than the linear case, and hence a lower option price, which makes it the less conservative extrapolation. It can still be a legitimate choice: a common practice on desks is flat vol extrapolation, which is the case $\alpha = 0$.

**Linear is the sensible choice.** Superlinear growth introduces butterfly arbitrage in the deep wings. A linear wing with a slope of at most 2 rises at the maximum rate that no-arbitrage allows, so it gives the most conservative price for deep out-of-the-money options. We extrapolate the wing by continuing total variance linearly in log-moneyness.

$$w(x) = w(x^\ast) + w'(x^\ast)\,(x - x^\ast), \qquad x > x^\ast,$$

This is how we address the negative density on the Dec strip: rather than read the $-0.10$ delta anchor off the spline, we stop at the last anchor before it, $-0.15$ delta, which becomes $x^\ast$, take the value and slope of the spline there, and extend that line out to $-0.10$ delta. The line inherits its slope from the spline at $x^\ast$, so we should check that $w'(x^\ast) \le 2$, and if it is steeper, use a slope of 2 instead. Once the wing is in place, complete the butterfly arbitrage check at the join point and over the extrapolated part.

## Reflection

In this article we constructed the implied volatility surface through a cubic spline. It is straightforward and reprices the market quotes fairly accurately, but it is not always stable: in our example the Dec strip's density came back choppy, and a rippled density unsettles the Greeks. The spline is not the only way to build a surface, and it is not always the right one. Where a desk cares more about stable Greeks and hedges than about the last basis point of repricing, the usual choice is a parametric model. One classic choice is SSVI, a smooth parametric form whose parameters admit simple conditions that rule out butterfly and calendar arbitrage. In practice, richer parametric models may be used to obtain both the smooth smile of a parametric form and repricing accuracy.

The other idea worth revisiting is the density function. We recovered it as the Breeden-Litzenberger second derivative, $q(K) = \partial^2 C/\partial K^2$, and used it to check for arbitrage, but the formula is more than a diagnostic:

- **It connects market quotes to models.** Option prices across strikes pin down the risk-neutral density at each maturity, which gives models a concrete target to match. For example, [Dupire's formula]({{< relref "local_vol.md#local-vol-calibration-and-the-dupire-formula" >}}) for local volatility uses this link directly.
- **It explains the cubic spline.** The call price is a smooth function of volatility, so differentiating it twice in strike brings in the second derivative of the smile, and a continuous density needs a smile whose second derivative is continuous too. A cubic is the lowest-degree spline that delivers this.
- **It links to the BKM spanning result.** In an [earlier article]({{< relref "vol_skewness_kurtosis.md" >}}), the BKM result writes the expectation of any twice-differentiable payoff as a continuum of calls and puts, and we proved it with a Taylor expansion. The density offers an alternative route, since the expectation of a payoff can be computed by integrating it against the density function.


## Appendix A: Sketch of the Total-Variance Butterfly Condition

The main text uses the butterfly condition written in total variance and log-moneyness,

$$g(x) = \left(1 - \frac{x\,w'}{2w}\right)^2 - \frac{w'^2}{4}\left(\frac{1}{w} + \frac14\right) + \frac{w''}{2} \;\ge\; 0.$$

This appendix sketches how it follows from the density condition $\partial^2 C/\partial K^2 \ge 0$ of the
density section.

Work in log-moneyness $x = \ln(F/K)$ and total variance $w(x) = \sigma^2(x)\,\tau$. The Black
call price on the future, written with these variables, is

$$C = F\left[N(d_1) - e^{-x} N(d_2)\right], \qquad d_{1,2} = \frac{x}{\sqrt{w}} \pm \frac{\sqrt{w}}{2},$$

using $\ln(F/K) = x$ and $\sigma\sqrt\tau = \sqrt{w}$.

The risk-neutral density in strike is $q(K) = \partial^2 C/\partial K^2$, and non-negativity
of $q$ is what we want to express through $w(x)$. Write $p(x)$ for the same density in
log-moneyness. Changing variable only introduces a factor of $1/K$,

$$q(K) = \frac{p(x)}{K},$$

which is strictly positive, so the sign is unaffected and we can work with $p(x)$ throughout.

Because $w$ depends on $x$, each differentiation splits into a direct part and a smile part.
Carrying both derivatives through and using $F\,n(d_1) = K\,n(d_2)$, every term collects onto
the same Black factor, leaving

$$p(x) = \frac{g(x)}{\sqrt{2\pi w}}\,\exp\!\left(-\frac{d_2^2}{2}\right),$$

where $g(x)$ gathers the terms in $w$, $w'$, and $w''$ produced by differentiating $d_1$ and
$d_2$. The factors multiplying $g(x)$ are strictly positive, so the density is non-negative
exactly when

$$g(x) = \left(1 - \frac{x\,w'}{2w}\right)^2 - \frac{w'^2}{4}\left(\frac{1}{w} + \frac14\right) + \frac{w''}{2} \;\ge\; 0.$$

Note that $g$ is a pointwise comparison, not one against a single lognormal. At each $x$,
it is the ratio of the density to the Black density built from the total variance $w(x)$ at
that same point, and that variance changes along the smile. Setting $w' = w'' = 0$ gives $g = 1$ and recovers the Black density:

$$p(x) = \frac{1}{\sqrt{2\pi w}}\,\exp\!\left(-\frac{d_2^2}{2}\right).$$
