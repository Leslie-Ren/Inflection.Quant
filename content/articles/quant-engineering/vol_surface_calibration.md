---
title: "Constructing the Implied Volatility Surface: From Market Quotes to an Arbitrage-Free Fit"
date: 2026-06-26
draft: false
math: true
tags: ["volatility-surface", "implied-volatility", "black-scholes", "svi", "no-arbitrage", "commodities"]
---

## Why This Matters

For vanilla options, the simple models are usually sufficient. A plain call or put can be priced off Black-Scholes directly; you often do not need to reach for local volatility or a stochastic volatility model. Those heavier models earn their place with exotics, where the payoff depends on how the smile behaves rather than just its level today. For a vanilla, you take the market's implied volatility at the relevant strike and maturity and feed it into Black-Scholes. But that assumes a volatility surface already exists: before Black-Scholes can price anything, the surface it reads from has to be built, and building it is less straightforward than it appears.

Fitting the smile to quotes is routine, but that is not the whole story. A surface can look smooth and still imply an arbitrage, so the fit needs to be checked rather than trusted. And the deep wings, where the market quotes thin out or vanish, need a construction of their own. Both are part of building the surface, and this article works through the full construction.

## Crude Oil as a Working Example

It helps to fix a concrete product and carry it through every step. WTI, traded on CME, is the most liquid crude benchmark and its options surface is the more richly quoted. Brent trades on ICE with good liquidity in the near-dated contracts, thinning gradually in the back months and deep wings. There is more than one way to build its surface. One approach calibrates directly from Brent's own quotes, which reflects Brent's own skew and supply-demand. Another builds Brent relative to WTI through a Brent-WTI volatility spread. The two benchmarks are highly correlated, so a spread can be more stable than an independent fit, and it leans on WTI's deeper quotes where Brent runs sparse. The cost is that the spread itself has to be modeled, and its assumptions carry their own risk. Which approach is more suitable depends on the book, and on whether the risk reaches into the regions where Brent's quotes run out.

For this article I calibrate the Brent surface directly and stay in the liquid near-dated region, where the quotes support it. The instrument is the [ICE Brent Crude option](https://www.ice.com/products/218/Brent-Crude-American-style-Options), and two features of its contract make the pricing clean enough to keep the focus on surface construction.

First, the option is written on a future, so the future serves directly as the forward the Black model takes as its input, with no spot to carry. Second, the option is futures-style margined. The premium is not paid up front; it is posted as margin and marked to market daily, like the future itself. As shown in a previous article on [options under futures-style margining]({{< relref "future_style_margining_options.md" >}}), this has two consequences that simplify the pricing: the exchange quote is already the undiscounted option value, so no discount factor appears, and the American option carries no early-exercise premium, so it can be priced as European.

Together these collapse the pricing to the Black model on the future, with no discount factor on either term:

$$C = F\,N(d_1) - K\,N(d_2), \qquad P = K\,N(-d_2) - F\,N(-d_1), \qquad d_{1,2} = \frac{\ln(F/K) \pm \tfrac12 \sigma^2 \tau}{\sigma \sqrt{\tau}}.$$

The Brent smile is worth picturing before we start. Crude has fat tails on both wings, and the skew can run either way depending on conditions: whichever tail the market fears more lifts that wing, a supply shock steepening the calls, a demand collapse steepening the puts.

## Building the Volatility Surface

### 1. Choosing the Anchor Grid

A surface is a function of two variables: across strike, where along the smile each volatility sits, and across maturity, which expiries the surface stores.

**The strike axis.** Where along the smile should the surface store its volatilities? Not at fixed strikes: as the future moves, a fixed strike drifts to a different moneyness. A common choice is to anchor at fixed deltas instead, a stable and portable coordinate. In practice the smile only stays roughly fixed in delta, not perfectly, but as a coordinate to store and interpolate volatilities it does the job.

With delta chosen as the coordinate, the next question is which deltas, and how many. FX settles this by convention with a standard set of quoted deltas; Brent has no such standard, so the grid is left to the modeler as a tradeoff. More anchors capture more skew and curvature but are only as good as what supports them (a quote or a wing extrapolation), and spacing them too finely makes the fitted smile wigglier, which shows up as greeks that jump between recalibrations. Fewer anchors are more stable but can miss curvature between them. Where the grid lands is largely a modeling decision, driven mainly by the book: one carrying deep-wing risk wants wider anchors to price and hedge it, even past where the market quotes. Other factors, such as downstream risk or pricing system setup, can also impact the choice.

For this article I will use a representative grid: closely spaced near the money and extending into each wing. The anchors are quoted as out-of-the-money options, calls above the forward and puts below:

$$10\Delta C,\; 15\Delta C,\; 25\Delta C,\; 35\Delta C,\; 50\Delta,\; 35\Delta P,\; 25\Delta P,\; 15\Delta P,\; 10\Delta P.$$

Each anchor can also be written as a put delta, which sweeps the whole smile monotonically:

$$-0.90,\; -0.85,\; -0.75,\; -0.65,\; -0.50,\; -0.35,\; -0.25,\; -0.15,\; -0.10.$$

**The maturity axis.** Which expiries should the surface store? Here the product decides for us. Brent options are listed on monthly futures, each expiring against its delivery month's future, so the natural maturity anchor is the contract month. For this article we take the first several contract months. In practice, we need to store enough contract months to cover a trading book's risk.

### 2. Selecting the Quotes

A volatility surface is only as good as the data that feeds it, and not every market quote should be used automatically. A quote earns its place only if it carries clean information about where the option trades. The criteria below are a reasonable starting set rather than an exhaustive list:

- **Use the out-of-the-money strike.** At any strike there is both a put and a call. The OTM one is the better quote to use: it is the more actively traded and tightly quoted of the two. At the strike equal to the future the two coincide: put-call parity makes the call and put equal there, so either serves, or their implied vols can be averaged.
- **Drop quotes below intrinsic.** No option can trade below its intrinsic value without implying an arbitrage.
- **Check monotonicity in strike.** Call prices must not rise as the strike rises, and puts must not fall. Only a strict reversal is an arbitrage; equal adjacent prices do not violate it, though they may still be dropped for other reasons.
- **Drop quotes at the minimum tick.** A deep out-of-the-money option can be worth less than a single tick, so its price is usually floored at the minimum increment and carries almost no volatility information. The vol backed out of it is unreliable.
- **Filter on liquidity, if the feed carries it.** Open interest and volume indicate whether a strike is genuinely active; a strike with negligible activity can be dropped as stale.

### 3. Constructing Each Strip

We build one maturity strip at a time. The market gives prices at strikes, but the strip stores a volatility at each delta anchor, so the construction has to bridge from one to the other.

**Invert each quote.** The Black model cannot be solved for $\sigma$ in closed form, so each quote, taken as the mid of bid and ask, is inverted numerically. This is exactly the single-option problem from the previous article on [implied volatility]({{< relref "newton_vs_brent_vol_solver.md" >}}), where Newton's method and Brent's method were compared. Each quote goes in, its implied volatility comes out.

**Order into three vectors.** Suppose $n$ quotes survive selection for this maturity strip. With a volatility for each, compute two more quantities per strike: the log-moneyness $x = \ln(F/K)$, and the put delta

$$\Delta = N(d_1) - 1, \qquad d_1 = \frac{x + \tfrac12 \sigma^2 \tau}{\sigma \sqrt{\tau}}.$$

Sort all three by log-moneyness, giving three aligned vectors that describe the raw market smile:

$$\mathbf{x} = [x_1, x_2, \ldots, x_n], \qquad \boldsymbol{\sigma} = [\sigma_1, \sigma_2, \ldots, \sigma_n], \qquad \boldsymbol{\Delta} = [\Delta_1, \Delta_2, \ldots, \Delta_n].$$

Each index is one quote: its log-moneyness, its implied volatility, and the delta that volatility and strike imply, with $x_1 < x_2 < \cdots < x_n$.

The anchor volatilities are what the strip stores, one at each delta anchor. Before solving for them, two things about the data need checking.

> **Do the quotes cover the delta anchors?** An anchor's volatility can be interpolated only if it sits within the deltas the quotes span. Suppose we want the $-0.90$ anchor but the quotes reach only $-0.85$. There is nothing to interpolate from, so the anchor has to be extrapolated past the data, which is the wing problem the later section handles. For Brent, front-month options are deeply liquid, so coverage there is rarely a problem; it is the back of the curve, where trading thins, that is more likely to leave the deep wings uncovered.
>
> **Is delta monotone in strike?** Matching a delta anchor to a strike only makes sense if each delta corresponds to a single strike. On a well-behaved smile this holds automatically, and it breaks only if the vol rises steeply enough between two adjacent quotes that the volatility term in $d_1$ overwhelms the moneyness term, reversing delta so one delta maps to several strikes. A reversal is either a bad quote or a wing too steep to be arbitrage-free; either way the anchor has no unique strike, so the data is worth examining before proceeding rather than solving through it.

**Solve the volatility at each anchor delta.** Take a target delta, say $-0.25$. We want the volatility at the strike whose delta is $-0.25$, and here is the difficulty: to locate that strike we need the volatility there, but to read the volatility off the smile we need the strike. Neither comes first. To read a volatility at an arbitrary strike we first need the smile as a continuous curve, so interpolate the market quotes with a cubic spline of $\sigma$ against $x$. Log-moneyness is the natural axis to spline in. Volatility is measured on the log return, so absolute strike is the wrong coordinate. Delta sweeps the smile monotonically but is itself a function of the volatility being fit, which makes it a less clean axis than $x$. With that continuous $\sigma(x)$ in hand, the way out of the circularity is a fixed-point loop: guess a volatility, derive the strike it implies, read the spline there for a new volatility, and repeat until the guess stops moving. The steps below are one iteration.

a. **Warm start.** Take the volatility of the market quote whose delta is nearest the target as the initial guess $\sigma^{(0)}$. Interpolating the volatilities of the two quotes that bracket the target delta gives a slightly better start, though the iteration converges quickly enough that the choice matters little.

b. **Solve for the strike ratio.** With the target delta and the current volatility $\sigma^{(k)}$ fixed, the delta equation inverts in closed form. With $d_1^\ast = N^{-1}(\Delta_{\text{target}} + 1)$,

$$x^{(k)} = d_1^\ast\,\sigma^{(k)}\sqrt{\tau} - \tfrac12 \big(\sigma^{(k)}\big)^2 \tau.$$

c. **Read the smile.** Evaluate the spline at $x^{(k)}$ to get a new volatility $\sigma^{(k+1)} = \sigma\big(x^{(k)}\big)$.

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
        s_new = float(smile(x_n))                     # vol the smile assigns there
        d1 = (x_n + 0.5*s_new**2*tau) / (s_new*np.sqrt(tau))
        delta = norm.cdf(d1) - 1.0
        if abs(s_new - s) < tol_vol and abs(delta - target_delta) < tol_delta:
            return s_new
        s = s_new
    raise RuntimeError("anchor solve did not converge")  # fall back to bracketed solve
```

This iteration is fast but not unconditionally convergent: it settles when the smile is gently sloped at the anchor, which usually holds across the covered range in practice, and it can slow where the wing steepens. The convergence condition is derived in [Appendix A](#appendix-a-convergence-of-the-anchor-solve). When the iteration does stall, the fallback is a bracketed root-find on $x$: find two market points whose deltas straddle the target and solve with Brent's method.

### 4. Assembling the Surface

Repeat step 3 for each stored maturity and lay the strips side by side; that grid of strips is the surface.

Reading a price off it takes two lookups:

- **Within a strip (one maturity).** A strike's volatility comes from translating it to log-moneyness $x = \ln(F/K)$, then evaluating that strip's cubic spline at that $x$, the same spline built in step 3.
- **Across strips.** Each strip belongs to a specific contract month's future, so what matters is which future an option settles against, not where its expiry sits in calendar time. An OTC option that settles against the September future reads its volatility off the September strip, with the option's own $\tau$ in the pricing, even if it expires a little before the exchange-traded September options do. That treats August and September as distinct underlyings, which they are, rather than blending two distributions by calendar distance. The contrast is a single underlying seen at different horizons, like a spot equity index, where the strips share one distribution stretched over time, so a real maturity axis exists to interpolate across, and there the right variable is total variance, since it is total variance, not volatility, that drives the dispersion of the price at each maturity.

## Spot the Hidden Issues

Suppose that, after going through construction steps 1 through 4, we end up with a surface that looks like this.

{{< vol_surface_spline >}}

It looks fine. Smooth strips, no kinks, every quote reproduced. But reproducing the quotes is not the same as being arbitrage-free, and the question we actually need to answer is whether the prices this surface implies are prices a rational market could support.

### Negative Density and Butterfly Arbitrage

The most basic property any set of option prices must satisfy is that the risk-neutral density it implies stays non-negative everywhere. So before trusting a strip, recover its implied density and check the sign.

But how do we translate the volatility surface into a density? Not through the lognormal density: the smile exists precisely because the distribution is not lognormal, with fatter tails on both sides. There is instead a model-free way, reading the density straight out of the call prices. It is the second derivative of the call price in strike,

$$q(K) = \frac{\partial^2 C}{\partial K^2}.$$

To see why, take the first derivative first. Move the strike from $K$ to $K + dK$. Every option that finishes out of the money still pays zero, so those scenarios are unchanged; every option that finishes in the money now pays $dK$ less. The drop in the call value is $dK$ weighted by the probability of landing in the money, $dC = -\,\mathbb{P}(F > K)\,dK$, so

$$\frac{\partial C}{\partial K} = -\,\mathbb{P}(F > K).$$

Now differentiate again. Moving $K$ to $K + dK$ changes the in-the-money probability by the mass in that band, and that mass over $dK$ is the density at $K$:[^density_df]

$$\frac{\partial^2 C}{\partial K^2} = -\,\frac{\partial}{\partial K}\mathbb{P}(F > K) = \frac{\mathbb{P}(K < F < K + dK)}{dK} = q(K).$$

This identity is the Breeden-Litzenberger formula. In practice this second derivative can be approximated numerically, on a grid of strikes spaced $\Delta K$ apart, by the central difference

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

Its density crosses below zero in the deep downside wing, near a strike of 53, bottoming around $-0.0014$. On the vol plot this is the $-0.10$-delta anchor, the far end of the Dec strip where the smile bends up into its wing. The bend is a touch too sharp there, and past a point extra curvature in the smile turns the call price concave in strike, which is a negative density. The excursion is likely too small and too deep in the wing to trade against once bid-ask and fees are counted, so it is not free money. It is still a defect to fix: the surface is not arbitrage-free as constructed.

The same strip also shows a choppy density closer to the money, a peak near 62, a dip near 66, another peak near 71. These sit around the $-0.25$-delta anchor, where the Dec smile bends up most sharply. The smile there looks smooth, but the density reads its curvature, not its level, and a cubic spline's curvature is far less smooth than the curve. Every value stays positive, so no butterfly is mispriced. However, the choppiness makes the Greeks and hedge ratios unstable.

### Calendar Arbitrage Across Strips

The butterfly check validates each strip on its own. The other classic check looks across strips: calendar arbitrage. At a fixed moneyness, total variance $w = \sigma^2\tau$ must not decrease with maturity, because a longer-dated option is the same price process given more time to diffuse and cannot carry less variance than a shorter one. If it did, a calendar spread across the two would lock in a riskless profit. The comparison has to be at the same moneyness on each strip, not the same delta, since a fixed delta maps to a different strike at each maturity.

This assumes one underlying seen at two horizons, which holds for an equity index but not strictly for Brent, where each strip is its own future and the Sep and Oct futures are distinct underlyings. Correlation is what makes the check still worth running: front-month Brent futures move nearly one-for-one, so across the near months they behave close enough to a single underlying that the condition applies. For commodities it is a heuristic, not a strict law, and it can fail on volatility seasonality: a January natural gas contract expiring into peak winter demand carries higher volatility than a following April contract, so total variance falls from the winter month to the spring one with no arbitrage, since the two are different bets on different weather regimes.

For this surface the strips are the near months, where the check applies. Comparing total variance at fixed moneyness across maturities, it rises monotonically from Sep 2026 through Jan 2027.

{{< calendar_total_variance >}}

The surface is clean across strips as well as within them.

[^density_df]: These clean forms rely on the futures-style margining set up earlier: the quote is the undiscounted option value, so $C(K) = \mathbb{E}[(F-K)^+]$ with no discount factor. For a premium-paid option, $C(K) = e^{-r\tau}\,\mathbb{E}[(F-K)^+]$, and the discount factor carries through both derivatives, $\partial C/\partial K = -e^{-r\tau}\,\mathbb{P}(F > K)$ and $\partial^2 C/\partial K^2 = e^{-r\tau}\,q(K)$. Recovering a unit-mass density then means dividing by $e^{-r\tau}$.

## What Happens to the Wings

So far we have built the surface assuming the market provides liquid quotes covering the delta anchors. That is not always the case. A position may need vol at a strike so far out of the money the exchange lists no option there at all, and the deep quotes that do exist are often stale, so they fail the selection in the first place. Either way the strip has to say something about a region the data does not cover, and the cubic spline is the wrong tool for it: a spline interpolates, and past the last anchor it simply continues its end cubic, a curve with no reason to resemble a volatility wing and every tendency to run off.

So how do we extrapolate the wings? What is a reasonable shape, one that at least keeps the density non-negative, the butterfly condition we checked earlier?

### How Variance May Grow in the Wing

Rather than study volatility against log-moneyness $x = \ln(F/K)$, work in total variance $w(x) = \sigma^2(x)\,\tau$, since variance is what drives the dispersion of returns. The question is how fast $w$ may grow as $x$ runs out into the wing.

The butterfly condition from the density section can be written entirely in these coordinates. Non-negative density is equivalent to a single function of $w$ and its first two derivatives staying non-negative at every $x$,

$$g(x) = \left(1 - \frac{x\,w'}{2w}\right)^2 - \frac{w'^2}{4}\left(\frac{1}{w} + \frac14\right) + \frac{w''}{2} \;\ge\; 0,$$

where $w' = dw/dx$ and $w'' = d^2w/dx^2$; its derivation is in [Appendix B](#appendix-b-deriving-the-total-variance-butterfly-condition). Test the growth by setting $w = Cx^\alpha$ and reading the sign of $g$ as $x \to \infty$.

**Faster than linear, $\alpha > 1$.** With $w = Cx^\alpha$, $w' = C\alpha x^{\alpha-1}$, $w'' = C\alpha(\alpha-1)x^{\alpha-2}$, the two curvature terms give

$$-\frac{w'^2}{4w} + \frac{w''}{2} = \frac{C\alpha(\alpha-2)}{4}\,x^{\alpha-2},$$

which for $\alpha > 1$ grows without bound and drives $g \to -\infty$. Superlinear growth is a butterfly arbitrage in the wing.

**Linear, $\alpha = 1$.** For $w = a + bx$, $g \to \tfrac14 - b^2/16$ as $x \to \infty$, non-negative only when $b \le 2$. A linear wing is arbitrage-free up to that slope; steeper reintroduces the violation. This is Lee's bound, and the construction has to respect it.

**Slower than linear, $\alpha < 1$.** Here $g$ stays non-negative, so no arbitrage, but $w'(x) = C\alpha x^{\alpha-1} \to \infty$ as $x \to 0$: an infinite slope at the money that no real smile has. A sublinear wing can still extrapolate without admitting arbitrage, but it is not a natural variance shape and not a standard choice.

**Linear is the sensible choice.** It is the boundary between the two: fast enough to keep a fat, arbitrage-free tail, slow enough to stay smooth through the money. We extrapolate the wing by continuing total variance linearly in log-moneyness.

$$w(x) = w(x^\ast) + w'(x^\ast)\,(x - x^\ast), \qquad x > x^\ast,$$

past a point $x^\ast$ on each side where we stop trusting the spline and let the line take over. For example, to extrapolate a 1-delta anchor, we take the value and slope of the spline at the last trusted anchor, the 10-delta, and extend that line out to where 1-delta falls.

## Reflection

In this article we constructed the implied volatility surface through a cubic spline. It is straightforward and reprices the market quotes fairly accurately, but it is not always stable: in our example the Dec strip's density came back choppy, and a rippled density unsettles the Greeks. The spline is not the only way to build a surface, and it is not always the right one. Where a desk cares more about stable Greeks and hedges than about the last basis point of repricing, the usual choice is a parametric model. One classic choice is SSVI, whose parametric form is smooth and arbitrage-free by construction. In practice, richer parametric models may be used to obtain both the smooth smile of a parametric form and repricing accuracy.

The other thread worth pulling runs back to the density function. We recovered it as the Breeden-Litzenberger second derivative, $q(K) = \partial^2 C/\partial K^2$, and used it to check for arbitrage. This reminds me of another identity we discussed in an [earlier article ]({{< relref "vol_skewness_kurtosis.md" >}}), where the BKM result writes any expectation of a payoff as a continuum of calls and puts. There we proved it with a Taylor expansion; an alternative route to reach the same identity uses the density function discussed here. The same second derivative that gives the density also underlies that spanning; the two are connected, and the connection is worth working out.


## Appendix A: Convergence of the Anchor Solve

The anchor solve iterates on the volatility: guess a value, derive the strike it implies, read the smile there, repeat. Written as a map, $\sigma^{(k+1)} = \Phi(\sigma^{(k)})$ with $\Phi(\sigma) = \mathrm{cs}\big(x(\sigma)\big)$, where $\mathrm{cs}(\cdot)$ is the cubic spline through the market smile and $x(\sigma)$ is the closed-form strike solve at the fixed target delta. When does this settle?

Let $\sigma^\ast$ be the answer, the fixed point $\sigma^\ast = \Phi(\sigma^\ast)$, and $e_k = \sigma^{(k)} - \sigma^\ast$ the error. Subtracting the fixed-point equation from the iteration gives $e_{k+1} = \Phi(\sigma^{(k)}) - \Phi(\sigma^\ast)$, and linearising $\Phi$ near $\sigma^\ast$ leaves

$$e_{k+1} \approx \Phi'(\sigma^\ast)\,e_k.$$

The error is multiplied by $\Phi'(\sigma^\ast)$ each step, so the iteration converges when $|\Phi'(\sigma^\ast)| < 1$, and faster the smaller that number.

What is $\Phi'$? It is a composition, so the chain rule gives $\Phi'(\sigma) = \mathrm{cs}'\big(x(\sigma)\big)\,x'(\sigma)$. The strike solve is $x(\sigma) = d_1^\ast\,\sigma\sqrt{\tau} - \tfrac12\sigma^2\tau$ with $d_1^\ast = N^{-1}(\Delta_{\text{target}}+1)$ constant, so $x'(\sigma) = d_1^\ast\sqrt{\tau} - \sigma\tau$. At the fixed point,

$$\Phi'(\sigma^\ast) = \mathrm{cs}'(x^\ast)\,\big(d_1^\ast\sqrt{\tau} - \sigma^\ast\tau\big).$$

The controlling factor is $\mathrm{cs}'(x^\ast)$, the slope of the smile at the anchor; the second factor is set by the maturity and the target delta. Near the money the slope is small and the iteration converges in a few steps; in the wings it steepens and the iteration slows, and once $|\Phi'|$ reaches one it no longer converges.

## Appendix B: Deriving the Total-Variance Butterfly Condition

The main text uses the butterfly condition written in total variance and log-moneyness,

$$g(x) = \left(1 - \frac{x\,w'}{2w}\right)^2 - \frac{w'^2}{4}\left(\frac{1}{w} + \frac14\right) + \frac{w''}{2} \;\ge\; 0.$$

This appendix derives it from the density condition $\partial^2 C/\partial K^2 \ge 0$ of the density section.

Work in log-moneyness $x = \ln(F/K)$ and total variance $w(x) = \sigma^2(x)\,\tau$. The Black call price on the future, written with these variables, is

$$C = F\left[N(d_1) - e^{-x} N(d_2)\right], \qquad d_{1,2} = -\frac{x}{\sqrt{w}} \pm \frac{\sqrt{w}}{2},$$

using $\ln(F/K) = x$ and $\sigma\sqrt\tau = \sqrt{w}$. The risk-neutral density is $q(K) = \partial^2 C/\partial K^2$, and non-negativity of $q$ is what we want to express through $w(x)$.

Changing variables from $K$ to $x$ (with $dx/dK = -1/K$) and carrying the two derivatives through, the density can be written as

$$q(K) = \frac{g(x)}{\sqrt{8\pi\,w}\;K^2}\,\exp\!\left(-\frac{d_2^2}{2}\right),$$

where $g(x)$ collects the terms in $w$, $w'$, and $w''$ that come out of differentiating $d_1$ and $d_2$ twice. The prefactor $\big(\sqrt{8\pi w}\,K^2\big)^{-1}\exp(-d_2^2/2)$ is strictly positive, so the sign of the density is the sign of $g(x)$. Collecting those terms gives

$$g(x) = \left(1 - \frac{x\,w'}{2w}\right)^2 - \frac{w'^2}{4}\left(\frac{1}{w} + \frac14\right) + \frac{w''}{2}.$$


[^density_df]: These clean forms rely on the futures-style margining set up earlier: the quote is the undiscounted option value, so $C(K) = \mathbb{E}[(F-K)^+]$ with no discount factor. For a premium-paid option, $C(K) = e^{-r\tau}\,\mathbb{E}[(F-K)^+]$, and the discount factor carries through both derivatives, $\partial C/\partial K = -e^{-r\tau}\,\mathbb{P}(F > K)$ and $\partial^2 C/\partial K^2 = e^{-r\tau}\,q(K)$. Recovering a unit-mass density then means dividing by $e^{-r\tau}$.