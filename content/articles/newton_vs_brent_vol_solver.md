---
title: "Solving for Implied Volatility: Making Newton's Method Production-Ready"
date: 2026-04-09
draft: false
math: true
tags: ["Root Finding", "Newton-Raphson", "Brent's Method", "Secant Method", "Bisection Method", "Inverse Quadratic Interpolation", "Convergence Rate"]
---

## Why This Matters

I plan to write about calibrating a full volatility surface, and that work rests on a clear understanding of a simpler problem: solving for implied vol from a single option price. At its core, this is a root-finding problem. Given a market quote, we look for the volatility at which the model price matches the market price. A common algorithm for this is Newton's method. I learned it in class, saw it in textbooks, and implemented it in the traders' tools at work. Most of the time it just works. I knew convergence wasn't guaranteed, but that stayed abstract until I hit an edge case.

This article is about enhancing the basic Newton's method and making it robust. The iteration itself takes only a few lines of code. A production solver needs more safeguards to guarantee convergence. We also discuss what to use when Newton is not an option. The method needs the derivative of the price with respect to volatility. Under Black-Scholes this is simply vega, but many pricers cannot compute it analytically. For those cases, Brent's method is a reliable choice, so it's worth knowing how it works as well.

## Newton's Method

Newton's method linearizes the function around the current guess and steps to the root of that linear approximation. For implied volatility, define the pricing error:

$$f(\sigma) = C_\text{BS}(\sigma) - C_\text{mkt}$$

The first-order Taylor expansion

$$f(\sigma) \approx f(\sigma_n) + f'(\sigma_n)(\sigma - \sigma_n)$$

set to zero gives the update

$$\sigma_{n+1} = \sigma_n - \frac{f(\sigma_n)}{f'(\sigma_n)} = \sigma_n + \frac{C_\text{mkt} - C_\text{BS}(\sigma_n)}{\mathcal{V}(\sigma_n)}$$

where vega $\mathcal{V} = \partial C_\text{BS}/\partial \sigma$ plays the role of $f'$. Geometrically, each step moves to where the tangent line crosses zero.

### Quadratic Convergence

Let $\hat{\sigma}$ be the implied vol, so $f(\hat{\sigma}) = 0$, and define the error $e_n = \sigma_n - \hat{\sigma}$. Because vega is strictly positive, $f'(\hat{\sigma}) \neq 0$, and near $\hat{\sigma}$ Newton converges quadratically. Expanding $f(\hat{\sigma})$ around $\sigma_n$:

$$0 = f(\hat{\sigma}) = f(\sigma_n) - f'(\sigma_n)e_n + \frac{1}{2} f''(\xi_n)e_n^2$$

for some $\xi_n$ between $\sigma_n$ and $\hat{\sigma}$. Dividing by $f'(\sigma_n)$ and substituting into the update gives

$$e_{n+1} = \frac{f''(\xi_n)}{2f'(\sigma_n)}e_n^2$$

Once near the root, the number of correct digits roughly doubles each iteration.

## Newton's Failure Modes and Their Fixes

The Newton step makes two implicit assumptions: vega is large enough that dividing by it gives a reasonable step, and the current guess is close enough to the root that the first-order Taylor expansion is a good approximation. Breaking either assumption can stop Newton from converging.

### Vega Near Zero: Step Size Instability

Vega sits in the denominator of the step. When $\mathcal{V} \approx 0$, a small pricing error produces a huge volatility update, throwing the iterate far from the solution. Vega is small when the option is deep in or out of the money, or when time to maturity is very short. In these regimes the price is close to intrinsic value and barely responds to volatility.

One way to fix this is to switch from Newton to bisection when vega drops below a threshold. Bisection is simple. Start with an interval that surely contains the root, such as $[10^{-4}, 5]$. If the price at the current guess is above the market price, the root is below the guess, so we keep $[10^{-4}, \sigma_n]$ and its midpoint becomes the next guess. If the price is below the market price, the root is above the guess, so we keep $[\sigma_n, 5]$ and take its midpoint instead. Each step halves the interval, so bisection always converges, and the solver returns to Newton as soon as vega recovers above the threshold.

### Poor Initial Guess: Curvature Breakdown

Newton truncates the Taylor expansion

$$C_{\text{BS}}(\sigma + \Delta\sigma) \approx C_{\text{BS}}(\sigma) + \mathcal{V}\Delta\sigma + \frac{1}{2}\frac{\partial^2 C_{\text{BS}}}{\partial \sigma^2}(\Delta\sigma)^2 + \cdots$$

after the linear term. The omitted second-order term grows with the square of the step $\Delta\sigma$. When the guess is close to the true implied vol, the step is small and this term is negligible. When the guess is far off, the step is large, the curvature term dominates, and the tangent line no longer predicts where the price curve crosses the market price. The step can overshoot the root, or even land at a negative volatility, where the iteration breaks down.

The interactive chart below lets you drag the true implied vol and the initial guess to see this behaviour. It uses a European call with $S = 100$, $K = 130$, $r = 5\%$, $T = 0.5$ years, and no dividends. For example, set the true implied vol to 30% and the initial guess to 12%. The first step jumps to about 259%, and the second lands at a negative volatility, so the method fails to converge.

{{< newton_method_initial_guess>}}

A natural fix is to come up with a smart initial guess that is closer to the true implied vol. The idea is to approximate the option price with a simple function of volatility that can be inverted in closed form. Two common methods are the Brenner-Subrahmanyam and Corrado-Miller approximations.

The Brenner-Subrahmanyam approximation ([discussed in this article]({{< relref "understanding_brownian_motion.md#connection-to-brennersubrahmanyam" >}})) assumes an at-the-money option and zero rates. Expanding $N(d_1)$ and $N(d_2)$ to first order around zero makes the call price linear in volatility:

$$C \approx \frac{S\sigma\sqrt{T}}{\sqrt{2\pi}} \quad\Longrightarrow\quad \sigma_0 = \sqrt{\frac{2\pi}{T}}\,\frac{C_\text{mkt}}{S}$$

It works well for near-the-money options with short maturities, where $\sigma\sqrt{T}$ is small enough for the linear approximation to hold. It ignores moneyness, and it degrades for long-dated options, where $\sigma\sqrt{T}$ grows and rates can no longer be ignored.

The Corrado-Miller approximation keeps more terms of the expansion and brings in moneyness through the discounted strike $X = Ke^{-rT}$:

$$C \approx \frac{S - X}{2} + \frac{(S + X)\,\sigma\sqrt{T}}{2\sqrt{2\pi}} + \frac{(S - X)^2}{\sqrt{2\pi}\,(S + X)\,\sigma\sqrt{T}}$$

Multiplying through by $\sigma\sqrt{T}$ turns this into a quadratic in $\sigma\sqrt{T}$, and taking the larger root gives

$$\sigma_0\sqrt{T} = \frac{\sqrt{2\pi}}{S + X}\left[C_\text{mkt} - \frac{S - X}{2} + \sqrt{\left(C_\text{mkt} - \frac{S - X}{2}\right)^2 - \frac{(S - X)^2}{\pi}}\,\right]$$

When the strike equals the forward, $K = Se^{rT}$, the moneyness terms vanish and Corrado-Miller reduces exactly to Brenner-Subrahmanyam. Away from that point, the extra terms keep it accurate over a much wider range of strikes and maturities. The catch is the square root: for deep out-of-the-money options its argument can turn negative. When that happens, fall back to the Brenner-Subrahmanyam estimate.

The solver should also check each Newton step against the bisection interval from the previous fix. A step that lands outside the admissible range is rejected and replaced by a bisection step. The solver can never go negative or wander off, and near the root it takes pure Newton steps at full quadratic speed.

## A Production-Ready Newton Solver

With both failure modes handled, two practical details complete the solver: screening the quote before solving, and deciding when to stop.

### Check the Price Before Solving

An implied vol exists only if the market price lies strictly inside the no-arbitrage bounds. For a European call:

$$\max\left(S - Ke^{-rT},\ 0\right) < C_\text{mkt} < S$$

Quotes at or very close to the lower bound are not uncommon: the market often quotes deep in-the-money options close to expiry within a tick of their intrinsic value. The remaining time value is effectively zero, so no meaningful implied vol can be solved, and such quotes are usually discarded when building a vol surface. Without this check, Newton oscillates until it hits the iteration cap and returns the final iterated value. Screen the quote up front and flag it, so a bad quote doesn't become a bad volatility downstream.

### Stopping on Price vs Vol

A root solver needs a condition that says when the result is good enough. The natural choice is the price, for example stop when $|C_\text{BS}(\sigma_n) - C_\text{mkt}| < 10^{-6}$. But a small price error is not necessarily a small vol error. The two are linked through vega: expanding the price to first order around the implied vol $\hat\sigma$ gives

$$f(\sigma_n) = C_\text{BS}(\sigma_n) - C_\text{mkt} \approx \mathcal{V}\,(\sigma_n - \hat\sigma) = \mathcal{V}\,e_n \quad\Longrightarrow\quad e_n \approx \frac{f(\sigma_n)}{\mathcal{V}}$$

so the vol error is the price error divided by vega. Deep in the wing, where vega is around $10^{-4}$, a price error of $10^{-6}$ means a vol error of $10^{-2}$, a full vol point.

A price match is enough if repricing the vanilla option is the only purpose. But in many cases the vol surface is not the final product: it is used to price other instruments, or as an input to a local vol model. A vol error in the wings then carries into these downstream uses.

So a threshold on vol is generally preferred: stop when the step in $\sigma$ (or the interval width) is below a vol tolerance such as $10^{-6}$. Keep a price check as a secondary condition and an iteration cap as a last resort.

### Putting It Together

```python
def implied_vol(c_mkt, S, K, r, T, tol=1e-6, max_iter=100):
    # Validate: no implied vol outside the no-arbitrage bounds
    lower = max(S - K * exp(-r * T), 0.0)
    if not (lower < c_mkt < S):
        raise ValueError("price outside no-arbitrage bounds")

    # Thresholds, interval and initial guess
    vega_min = 1e-6 * S       # below this, vega is too small and Newton can oscillate
    price_tol = 1e-4 * S      # check that the result reprices the quote
    lo, hi = 1e-4, 5.0
    sigma = initial_guess(c_mkt, S, K, r, T)   # Brenner-Subrahmanyam / Corrado-Miller
    sigma = min(max(sigma, lo), hi)

    for _ in range(max_iter):
        f = bs_call(S, K, r, T, sigma) - c_mkt
        if f == 0.0:
            return sigma          # exact hit
        if f > 0:
            hi = sigma            # root is below the guess
        else:
            lo = sigma            # root is above the guess

        # Newton step, only if vega is usable and the step stays inside the interval
        vega = bs_vega(S, K, r, T, sigma)
        new = sigma - f / vega if vega > vega_min else None
        if new is None or not (lo < new < hi):
            new = 0.5 * (lo + hi) # bisection fallback

        # Stop on volatility
        done = abs(new - sigma) < tol or (hi - lo) < tol
        sigma = new
        if done:
            break
    else:
        raise RuntimeError("did not converge")

    # Check on price: the vol must reprice the quote
    if abs(bs_call(S, K, r, T, sigma) - c_mkt) > price_tol:
        raise RuntimeError("vol converged but does not reprice the quote")
    return sigma
```

With a good initial guess, a liquid near-the-money option converges in a few iterations and never needs the bracket method. If an iterate lands where vega is small, or a Newton step leaves the interval, the bisection method kicks in. This is most common on the wings.

## When the Analytical Derivative Is Not Available

The Newton solver above needs an analytical vega, which Black-Scholes conveniently provides. In practice we often need an implied vol where the price comes from a model with no analytical derivative. The most common case is an American option priced on a binomial tree or a PDE grid. Since most listed single-stock options are American, this case is routine. A second case is a vendor pricing library that returns a price but has no closed-form vega.

The obvious workaround is to approximate the derivative by finite differences:

$$f'(\sigma) \approx \frac{f(\sigma + h) - f(\sigma - h)}{2h}$$

This works, but it comes with weaknesses. Each Newton iteration now needs three pricer calls, $f(\sigma_n)$ and $f(\sigma_n \pm h)$, instead of one price plus a closed-form vega. And the step size $h$ has to be tuned: too large and the derivative is biased; too small and the pricer's numerical noise dominates the difference, so the derivative can be wildly wrong or even have the wrong sign.

Brent's method avoids both problems. It uses only function values, so each iteration needs a single pricer call and there is no $h$ to tune. It also keeps the root bracketed at every step, so a bad step caused by a noisy price cannot make it diverge.

## Brent's Method

Brent's method is an adaptive combination of three root finders: bisection, the secant method, and inverse quadratic interpolation. At each iteration it tries the fastest step that passes its safety checks and falls back to bisection when they fail. It needs one function evaluation per iteration and no derivative, it guarantees convergence, and when the function is well behaved it converges faster than linearly.

The method starts from a bracket $[\sigma_a, \sigma_b]$ with

$$f(\sigma_a) \cdot f(\sigma_b) < 0$$

so the pricing error changes sign across the interval and a root lies inside. Brent maintains this bracket at every iteration, in the same spirit as the interval kept by the Newton solver.

### The Three Building Blocks

#### Bisection: Slow but Always Converges

This is the same bisection used in the Newton solver's fallback for small vega. It evaluates $f$ at the midpoint of the bracket,

$$\sigma_\text{mid} = \frac{\sigma_a + \sigma_b}{2}$$

and replaces the endpoint whose $f$ has the same sign as $f(\sigma_\text{mid})$. The bracket halves every step, regardless of the shape of $f$.

**Convergence:** Bisection uses only the sign of $f$, and it discards all the other available information. Convergence is linear at rate $\frac{1}{2}$: reducing a bracket of width $W$ to tolerance $\varepsilon$ takes $\lceil \log_2(W/\varepsilon) \rceil$ iterations. For a bracket $[10^{-4}, 5]$ and tolerance $10^{-6}$, that is 23 iterations.

#### Secant Method: Newton Without the Derivative

The secant method draws a straight line through the two most recent iterates $(\sigma_{n-1}, f_{n-1})$ and $(\sigma_n, f_n)$ and steps to where it crosses zero:

$$\sigma_{n+1} = \sigma_n - f_n \frac{\sigma_n - \sigma_{n-1}}{f_n - f_{n-1}}$$

This is Newton with the derivative replaced by the slope between the last two points. The key difference from a finite difference Newton is that the slope reuses evaluations you already have, so each iteration costs one pricer call. Early on the points are also far apart, so pricer noise distorts the slope much less than with a small $h$.

**Convergence:** Near the root, a local analysis gives

$$e_{n+1} \approx C e_n e_{n-1}$$

where $C = f''(\hat\sigma)/\bigl(2f'(\hat\sigma)\bigr)$, the same constant as in Newton's error formula. This form is not surprising: the next iterate is built from the current and previous vols, so its error depends on both of their errors. This gives superlinear convergence of order $\varphi \approx 1.618$, the golden ratio. With an analytical derivative, Newton is faster. When the derivative has to be approximated by finite differences, the comparison changes. Newton still roughly doubles the number of correct digits per iteration, but each iteration now costs three pricer calls. The secant method multiplies the correct digits by about 1.6 per iteration at the cost of one call. For the same three calls, Newton doubles the digits once, while the secant method multiplies them by about $1.6^3 \approx 4$. Without a bracket, the secant step can overshoot when $f$ is strongly curved, so convergence is not guaranteed.

#### Inverse Quadratic Interpolation: Fast but Fragile

Inverse quadratic interpolation (IQI) fits a quadratic through the three most recent points $(\sigma_a, f_a)$, $(\sigma_b, f_b)$, $(\sigma_c, f_c)$. The key detail is that it fits $\sigma$ as a function of $f$, so evaluating at $f = 0$ gives the next iterate directly. By Lagrange interpolation:

$$\sigma_\text{new} = \sigma_a \frac{f_b f_c}{(f_a - f_b)(f_a - f_c)} + \sigma_b \frac{f_a f_c}{(f_b - f_a)(f_b - f_c)} + \sigma_c \frac{f_a f_b}{(f_c - f_a)(f_c - f_b)}$$

**Convergence:** Near the root, a local analysis gives

$$e_{n+1} \approx K\, e_n e_{n-1} e_{n-2}$$

where $K$ depends on the first three derivatives of $f$ at the root. As with the secant method, the next iterate is built from the last three vols, so its error depends on all three of their errors. The third point lets IQI capture local curvature, which the secant method ignores, giving convergence of order $q \approx 1.839$. However, when the function values are close together, the fit is ill-conditioned, and because IQI extrapolates to $f = 0$, the step can land far outside the bracket. That is why IQI is rarely used alone.

### Combining the Three Methods

Bisection always converges but is slow. The secant method and IQI are fast near the root but can fail when the function is badly behaved. Brent's idea is to let each method do what it does best: try a fast step first, check whether it can be trusted, and fall back to bisection when it cannot. Since the root stays inside a bracket at every step, convergence is guaranteed.

Each iteration runs through the following steps:

1. **Propose a fast step.** Brent keeps three points: the best estimate so far $\sigma_\text{best}$ (smallest $|f|$), the contrapoint $\sigma_\text{contra}$, with $f(\sigma_\text{contra})$ opposite in sign to $f(\sigma_\text{best})$ so the root lies between them, and the previous best estimate $\sigma_\text{prev}$. Use IQI through $\sigma_\text{contra}$, $\sigma_\text{best}$ and $\sigma_\text{prev}$ when the three function values are available and distinct. Otherwise use a secant step through $\sigma_\text{contra}$ and $\sigma_\text{best}$.
2. **Check the step.** Accept it only if both conditions hold:
    - It lands in the safe part of the bracket: close to $\sigma_\text{best}$, away from $\sigma_\text{contra}$. A step that jumps nearly across the bracket to $\sigma_\text{contra}$ usually means the interpolation is misbehaving, and it would barely shrink the bracket. In Brent's original algorithm, the safe part is the three quarters of the bracket nearest $\sigma_\text{best}$, that is, between $(3\sigma_\text{contra} + \sigma_\text{best})/4$ and $\sigma_\text{best}$.
    - It shows real progress: near the root, a converging method takes smaller and smaller steps. Brent requires them to shrink at least about as fast as bisection would, while tolerating an occasional larger step, since the fast methods do not shrink their steps evenly. In Brent's original algorithm, the step must be smaller than half the step taken two iterations earlier. If interpolation stalls, this check fails and forces bisection.

    If either fails, take the bisection midpoint $\sigma_\text{mid} = (\sigma_\text{contra} + \sigma_\text{best})/2$ instead.
3. **Evaluate the pricing error** $f(\sigma_\text{new}) = C(\sigma_\text{new}) - C_\text{mkt}$, the single pricer call of the iteration.
4. **Update the points.** Save the current $\sigma_\text{best}$ as $\sigma_\text{prev}$. Then the new point replaces one end of the bracket, chosen so the root stays inside: if $f(\sigma_\text{new})$ has the opposite sign to $f(\sigma_\text{contra})$, it replaces $\sigma_\text{best}$; otherwise it replaces $\sigma_\text{contra}$. Finally, if $|f(\sigma_\text{contra})| < |f(\sigma_\text{best})|$, the contrapoint is now the better estimate, so the two exchange roles: $\sigma_\text{contra}$ becomes $\sigma_\text{best}$, and the old $\sigma_\text{best}$ becomes $\sigma_\text{contra}$.
5. **Check for convergence.** Stop if $f(\sigma_\text{best}) = 0$ or the half-width $|\sigma_\text{best} - \sigma_\text{contra}|/2$ is below the tolerance. Otherwise return to step 1.

The bracket never loses the root, and the safeguards bound the worst case: Brent can never be much slower than bisection overall, and in practice it is usually far faster.

### Brent in Production

Brent's method shares the production checks of the Newton solver. The no-arbitrage check still comes first. If the quote passes it, an implied vol exists, and since the price increases with volatility, the interval $[10^{-4}, 5]$ used for the Newton fallback is also a valid starting bracket for Brent, with $f(10^{-4}) < 0 < f(5)$. For an American option, the lower bound in the no-arbitrage check is the immediate exercise value $S - K$ rather than the European bound.

The stopping rule carries over too: set the tolerance on $\sigma$, not on the price error, since a close price match can still hide a large vol error in the wings. With a tree or PDE pricer there is one more limit: the price is only known to within the pricer's own numerical error. Brent stays robust because in the worst case it only needs the sign of $f$, but the implied vol it returns inherits the pricer's error: a price error $\varepsilon$ becomes a vol error of about $\varepsilon / \mathcal{V}$.

Finally, Brent is slower than Newton near the root. When a closed-form derivative is available, Newton's quadratic rate wins, so Brent is the tool for when Newton is not an option, not a replacement for it.

## Conclusion

Safeguarded Newton's method and Brent's method share more than they differ. Both keep a bracket so the solver always converges, and both fall back to bisection when the fast step is rejected. The choice between them comes down to vega: use Newton when an analytical vega is available, and use Brent otherwise.

Newton's update uses only the first derivative. For Black-Scholes implied vol, Peter Jäckel's "Let's Be Rational" (2015) pushes the same recipe further. A carefully constructed approximation gives an accurate starting point across all strikes, and a Householder iteration, a generalization of Newton's method that also uses the second and third derivatives, reaches the maximum attainable precision in as few as two iterations. It is worth reading for anyone building a production-grade Black-Scholes solver.

## References

Brenner, M., and Subrahmanyam, M. G. (1988). A simple formula to compute the implied standard deviation. *Financial Analysts Journal*, 44(5), 80-83.

Brent, R. P. (1973). *Algorithms for Minimization without Derivatives*. Prentice-Hall.

Corrado, C. J., and Miller, T. W. (1996). A note on a simple, accurate formula to compute implied standard deviations. *Journal of Banking and Finance*, 20(3), 595-603.

Jäckel, P. (2015). Let's be rational. *Wilmott*, 2015(75), 40-53.
