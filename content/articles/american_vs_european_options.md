---
title: "Early Exercise of American Options: Call Equivalence and the Put Premium"
date: 2026-04-03
draft: false
math: true
tags: ["American Option", "European Option", "Put-Call Parity", "Black-Scholes PDE", "Dividend"]
---

## Why This Matters

One of the first results I learned in my derivatives pricing course is that an American call on a non-dividend-paying stock is worth the same as its European counterpart, while an American put can be worth more. The result is easy to remember, but what actually produces the asymmetry did not register with me at the time, and that is what I wanted to pin down in this article.

For both the call and the put, the decision comes down to the same comparison, what exercising delivers today against what holding is worth. Interest rates and dividends both impact that comparison, and the standard derivation makes simplifying assumptions about them, taking rates as positive and the stock as paying no dividends. I also wanted to work out what happens when the rate is negative, and when the stock does pay a dividend.

## Why Early Exercise of an American Call Is Suboptimal

Exercising an American call early delivers $S_t - K$ at time $t < T$. Whether that is ever worth doing depends on what the option is worth if the holder does not exercise. For European options on a non-dividend-paying stock, put-call parity rearranged to isolate the call says that a European call amounts to owning the stock, holding a European put for downside protection, and deferring payment of the strike until maturity:

$$C_t = S_t + P_t - Ke^{-r(T-t)}$$

An American call is worth at least as much as the European call, since its holder can always choose to wait, so this value is a lower bound on the American call's worth. Comparing it against the $S_t - K$ from exercising:

$$C_t - (S_t - K) = P_t + K\left(1 - e^{-r(T-t)}\right)$$

The two terms on the right are exactly what early exercise gives up. The first is the option's time value, the optionality it still carries. Exercising swaps the kinked payoff for the straight line $S - K$. Below $K$ that line runs negative, where the option's payoff would have stopped at zero, and that floor at zero is the protection $P_t$ measures. The second is the interest on $K$, earned by holding the strike until maturity instead of paying it today.

Under $r > 0$ and no dividends both terms are positive, so the American call is never exercised early, and American and European calls have the same price.

### What Happens When the Assumptions Fail

**When the rate is not positive.** At $r = 0$ the term $K(1 - e^{-r(T-t)})$ vanishes. There is no interest to earn on $K$, so that half of the argument is gone. The first term survives, since the protection below $K$ is worth something as long as volatility and time remain, and early exercise stays suboptimal on the strength of the optionality alone.

At $r < 0$ the term $K(1 - e^{-r(T-t)})$ turns negative. Paying $K$ later costs more than paying it today, so deferring the strike is now a cost rather than a benefit. The comparison becomes $P_t$ against a negative second term, and that sum can come out negative, so early exercise of a call can become optimal. That leaves a genuine trade-off between taking intrinsic value now and keeping the optionality, and a boundary where the two balance. It is the same trade-off that makes early exercise optimal for the put under positive rates, which we work out in the [American puts section](#why-american-puts-are-worth-more-than-european-puts).

**When the stock pays a dividend.** Keep $r > 0$ and let the stock pay a dividend. On the ex-dividend date the stock trades without the right to that dividend, and its price drops. In practice the drop is usually smaller than the dividend itself.[^tax] We assume the dividend $D$ equals the drop and is deterministic. This is a simplification, as dividends further out are unknown and need to be forecasted. Note that the option holder does not own the stock, only a claim on it, and is therefore not entitled to the dividend unless they exercise before the ex-dividend date to own the shares.

Today's stock price already contains the dividends the stock will pay before $T$. Buying it now gets you two things, the dividend stream and the stock's worth at expiry. Strip out the dividend, and the remaining $S_t - \text{PV}_t(D)$ is what you are paying for the terminal price $S_T$ alone. Here $\text{PV}_t(D)$ is the present value at $t$ of the dividends, each discounted from its own ex-dividend date:

$$\text{PV}_t(D) = \sum_i D_i e^{-r(t_i - t)}$$

where $t_i$ is the ex-dividend date of $D_i$. An option's payoff depends only on $S_T$, so put-call parity with dividends becomes:

$$C_t = S_t - \text{PV}_t(D) + P_t - Ke^{-r(T-t)}$$

Comparing against exercising, as before:

$$C_t - (S_t - K) = P_t + K\left(1 - e^{-r(T-t)}\right) - \text{PV}_t(D)$$

The dividend enters with a minus sign, since it is the one thing exercising gains rather than gives up, and early exercise is attractive only when it outweighs the first two terms.

The stock only drops on ex-dividend dates, so exercise can only be optimal immediately before one. Exercising any earlier captures the same dividend but pays the strike sooner and gives up more optionality, so waiting until the last moment dominates. Take a single dividend $D$ with ex-dividend date $t_d$, and write $t_d^-$ for the moment just before the drop and $t_d^+$ for the moment just after it, so that

$$S(t_d^+) = S(t_d^-) - D$$

The decision is made at $t_d^-$, the last point at which exercising still captures the dividend. Comparing the two choices just after the date:

| | Exercise at $t_d^-$ | Hold |
|---|---|---|
| Shares held | $S(t_d^+)$ | none |
| Dividend received | $D$ | none |
| Strike paid | $K$ | not yet |
| Option retained | none | a call on a stock with no dividends left |
| Value at $t_d^+$ | $S(t_d^-) - K$ | $C(t_d^+)$ |

Exercising nets to $S(t_d^-) - K$ because the shares are worth $D$ less after the drop and the dividend pays exactly that back. Holding keeps an option on a stock with no dividends left before expiry, so parity applies in its original form at the ex-dividend price:

$$C(t_d^+) = S(t_d^+) + P(t_d^+) - Ke^{-r(T-t_d)}$$

Dropping $P(t_d^+)$ gives a lower bound on what holding is worth, so holding beats exercising whenever

$$S(t_d^-) - D - Ke^{-r(T-t_d)} \geq S(t_d^-) - K$$

which gives $D \leq K(1 - e^{-r(T-t_d)})$. Over a short remaining life this is approximately $D/K \leq r(T-t_d)$, and with $K$ close to $S$ the left side is the dividend yield. When the dividend yield falls short of the interest, early exercise is ruled out. We simplified by taking a single dividend. With multiple dividends spread through the life of the option, the same test applies at each ex-dividend date, and $T - t_d$ is smallest at the last one. So with $D$ assumed constant, the condition $D/K \leq r(T - t_d)$ is hardest to satisfy at the final ex-dividend date. If early exercise is ruled out at the final ex-dividend date, it is also ruled out at the prior dates.

When the dividend yield exceeds the interest, early exercise can be worthwhile, though not automatically. Comparing the value of holding against the value of exercising:

$$C(t_d^+) - \left(S(t_d^-) - K\right) = P(t_d^+) + K\left(1 - e^{-r(T-t_d)}\right) - D$$

Exercise is optimal precisely when $D$ exceeds $P(t_d^+) + K(1 - e^{-r(T-t_d)})$, which happens when the put value is small, that is when the stock price is high.

## Why American Puts Are Worth More Than European Puts

With $r > 0$ and no dividends, exercising calls early gives up two things at once, the optionality and the interest on $K$. Nothing is received in return, so there is never a reason to exercise early. Put options are a different case. Exercising still gives up the optionality, but it collects $K - S$ in cash today rather than at expiry. Early exercise of the put therefore carries both a cost and a benefit, and early exercise can be optimal. Keeping $r > 0$ and no dividends for now, the Black-Scholes PDE is the clearest way to see where the balance tips.

Under the Black-Scholes assumptions, the value $P(S,t)$ of a put on a non-dividend-paying stock satisfies the following PDE as long as the option is held rather than exercised:

$$\underbrace{\frac{\partial P}{\partial t}}_{\text{Time Decay}} + \underbrace{\frac{1}{2}\sigma^2 S^2 \frac{\partial^2 P}{\partial S^2}}_{\text{Convexity Gain (Gamma)}} + \underbrace{rS\frac{\partial P}{\partial S}}_{\text{Drift of Underlying}} - \underbrace{rP}_{\text{Carry Cost}} = 0$$

This PDE is derived by constructing a delta-hedged portfolio and requiring that its value grows at the risk-free rate under no-arbitrage. Each term has a financial meaning:

- **Time Decay** $\frac{\partial P}{\partial t}$: the rate at which the option loses value as expiry approaches, holding $S$ fixed
- **Convexity Gain** $\frac{1}{2}\sigma^2 S^2 \frac{\partial^2 P}{\partial S^2}$: the gain from being long gamma, because the put is convex in $S$, the holder benefits on average from large moves in either direction
- **Drift of Underlying** $rS\frac{\partial P}{\partial S}$: since $\frac{\partial P}{\partial S} < 0$ for a put, the risk-neutral upward drift of the stock works against the put holder
- **Carry Cost** $-rP$: the opportunity cost of holding the option rather than investing its value at the risk-free rate

The PDE is a balance condition: the convexity gain from being long gamma exactly offsets the combined drag from time decay, adverse drift, and carry cost.

### The Hold Region and the Exercise Region

For a European put, the PDE holds everywhere. The holder has no choice but to wait. For an American put, the holder has agency, and the $(S, t)$ plane splits into two regions.

In the **hold (continuation) region**, where $P(S,t) > K - S$, the option is worth more alive than dead. The convexity gain is sufficient to justify the carry cost and the drag from drift. The PDE holds with **equality**:

$$\frac{\partial P}{\partial t} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 P}{\partial S^2} + rS\frac{\partial P}{\partial S} - rP = 0$$

In the **exercise (stopping) region**, where $P(S,t) = K - S$, the balance breaks down. Substituting $P = K - S$:

- $\frac{\partial P}{\partial t} = 0$, intrinsic value has no time decay
- $\frac{1}{2}\sigma^2 S^2 \frac{\partial^2 P}{\partial S^2} = 0$, intrinsic value is linear in $S$, so gamma is zero
- $rS\frac{\partial P}{\partial S} = -rS$, the delta of $K - S$ is $-1$
- $-rP = -r(K-S)$

The PDE evaluates to:

$$0 + 0 + (-rS) - r(K - S) = -rK < 0$$

The equality strictly fails. Gamma is zero, so nothing offsets the carry cost and the adverse drift, which together come to $rK$. With $r > 0$ the residual is strictly negative, and the PDE inequality signals that continuation is dominated by immediate exercise.

### Where the Exercise Boundary Sits

We have shown that an exercise region exists, but not where it is. At a given time $t$ before maturity, when $S$ is close to $0$ the payoff $K - S$ is already close to $K$, the largest amount the put can ever pay. The value of optionality, which comes from the convexity of the payoff, is close to zero, so the holder is better off exercising to collect $K - S$ now and earn interest on it. As $S$ increases, the payoff moves further away from its ceiling and the value of optionality increases. Somewhere in between, exercising and holding are exactly balanced, and that level of $S$ is what we want to pin down. Repeating this at every $t$ traces out the exercise boundary $S^*(t)$. For $S \le S^*(t)$ the put is exercised, and for $S > S^*(t)$ it is held.

At each $t$ the time to maturity $T - t$ is different, so the balance point moves with $t$. Let us first simplify the problem by assuming a perpetual put, with $T = \infty$.

**The perpetual put.** With no expiry the option looks the same at every date, so $\partial P/\partial t = 0$ and the PDE in the hold region becomes an ODE:

$$\frac{1}{2}\sigma^2 S^2 P'' + rSP' - rP = 0$$

Notice that in each term the power of $S$ matches the order of the derivative: $S^2$ multiplies $P''$, $S$ multiplies $P'$, and $S^0$ multiplies $P$. For a power $P = S^\beta$, each derivative lowers the power by one and each factor of $S$ raises it back, so every term becomes a constant times $S^\beta$. This makes $S^\beta$ the natural guess. Substituting it and dividing by $S^\beta$ turns the ODE into a quadratic in $\beta$:

$$\frac{1}{2}\sigma^2\beta(\beta - 1) + r\beta - r = (\beta - 1)\left(\frac{1}{2}\sigma^2\beta + r\right) = 0$$

The roots are $\beta = 1$ and $\beta = -\alpha$, where $\alpha = 2r/\sigma^2$. The ODE is linear and second order, so these two independent solutions give every solution:

$$P = AS^{-\alpha} + BS$$

The term $BS$ is simply $B$ shares of stock, which grows without bound. The put must vanish as $S \to \infty$, which forces $B = 0$, so $P = AS^{-\alpha}$ in the hold region. This leaves two unknowns, the constant $A$ and the exercise boundary $S^*$, so we need two conditions at the boundary.

1. **Value matching.** At the boundary the held option is worth exactly its intrinsic value, so the value is continuous as $S$ crosses from the hold region into the exercise region:

   $$AS^{*-\alpha} = K - S^*$$

2. **Delta matching.** At the boundary the delta of the option must be $-1$, a condition usually called smooth pasting:

   $$-\alpha A S^{*-\alpha-1} = -1$$

A delta of $-1$ means the option moves one for one against the stock, exactly like the payoff $K - S$ it is about to become. If the delta at the boundary were above $-1$, the value curve would meet the payoff line at a kink, and waiting would beat exercising there.

<details>
<summary>Why the delta has to be -1</summary>

A put's delta can never be below $-1$, so the only question is why it cannot be above $-1$ at the boundary. Suppose it were, with delta $\Delta > -1$ just above $S^*$, say $\Delta = -0.8$. Just below $S^*$ we fall into the exercise region, so the value is $K - S$ with slope $-1$. Just above it we are in the hold region, where the value curve has slope $\Delta$, so the two pieces meet at a kink. Stand exactly at $S^*$ and compare the two choices over a short time $dt$.

**Exercise now.** The holder receives $K - S^*$ in cash and invests it at the risk-free rate, gaining

$$r(K - S^*)\,dt$$

**Wait.** The holder keeps the option while the stock moves by

$$dS = rS^*\,dt + \sigma S^*\,dW, \qquad dW = \sqrt{dt}\,Z, \quad Z \sim N(0,1)$$

If the stock falls, the option value increases by $-dS$. If it rises, it changes by $\Delta\,dS$. Therefore

$$dP = -dS + (1 + \Delta)\max(dS, 0)$$

Since $E[dS] = rS^*\,dt$ and, to leading order, $E[\max(dS, 0)] = \sigma S^*\sqrt{dt}\,E[\max(Z, 0)] = \sigma S^*\sqrt{dt/2\pi}$, the expected change in the option value is

$$E[dP] = (1 + \Delta)\,\sigma S^*\sqrt{\frac{dt}{2\pi}} - rS^*\,dt$$

**Compare.** Subtracting what exercising earns from what waiting earns:

$$E[dP] - r(K - S^*)\,dt = \underbrace{(1 + \Delta)\,\sigma S^*\sqrt{\frac{dt}{2\pi}}}_{\text{benefit from the kink}} \;-\; \underbrace{rK\,dt}_{\text{interest on } K}$$

The benefit is of order $\sqrt{dt}$ and the cost is of order $dt$, and for small $dt$ the square root is far larger. With $\Delta = -0.8$, $\sigma = 30\%$, $S^* = 50$, $K = 100$, $r = 5\%$ and $dt = 0.0001$, the benefit is about $0.012$ against a cost of $0.0005$, and shrinking $dt$ only widens the gap. So for any $\Delta > -1$, waiting is strictly better than exercising at $S^*$, which contradicts $S^*$ being the boundary.

</details>

<br>

Solving the two conditions together gives the boundary:

$$S^*_\infty = \frac{\alpha}{1+\alpha}K = \frac{2r}{2r + \sigma^2}K$$

When $r = 0$ the boundary is at $0$ and the exercise region disappears, because there is no interest on $K$ to collect. As $r$ increases, the interest earned by receiving $K$ early grows, so exercising becomes more attractive and the boundary rises. When $\sigma = 0$ the boundary is at $K$, because there is no optionality left to give up. As $\sigma$ increases, the value of the optionality grows, so holding becomes more attractive and the boundary falls. The boundary never exceeds the strike, since exercising above $K$ has no payoff.

**Finite maturity.** Once the expiry is brought back, $\partial P/\partial t$ no longer drops out, the boundary has to be solved together with the PDE, and there is no closed form for $S^*(t)$. We can, however, bound it. A put with less time left has less optionality to give up than the perpetual put, so it is exercised sooner. The boundary therefore sits between the perpetual level and the strike:

$$S^*_\infty \le S^*(t) \le K$$

The boundary can be computed with a binomial tree that takes the larger of the continuation value and $K - S$ at every node. The figure shows it for $K = 100$ and $\sigma = 30\%$.

{{< exercise-boundary >}}

Far from expiry, each boundary drifts slowly down toward its perpetual level, and as expiry approaches it climbs quickly to the strike. A higher rate lifts the boundary at every maturity, and the gap between the two curves widens with time to expiry.

**The call under negative rates.** The same boundary analysis applies to the call when $r < 0$, the case where we found early exercise can become optimal. With $r < 0$, paying $K$ later costs more than paying it today, so exercising early saves that cost, and the saving is weighed against the optionality given up. The optionality shrinks as the call moves deeper in the money, so exercising wins when $S$ is high, which places the exercise region above the early exercise boundary rather than below it. As expiry approaches, the boundary falls quickly toward $K$ from above, just as the put's boundary rises toward $K$ from below in the figure above.

### What Happens When the Assumptions Fail

**When the rate is not positive.** The only benefit of exercising the put early is receiving $K - S$ today and earning interest on it until expiry. At $r = 0$ that benefit is gone, so exercising early would give up the optionality for nothing, and the American put is worth the same as the European put. At $r < 0$ it becomes a cost, since cash received early shrinks rather than grows, so exercising early loses on both counts. This is the mirror of the call under positive rates, where exercising early gives up the optionality and pays the strike too soon.

**When the stock pays a dividend.** Keep $r > 0$ and let the stock pay a single dividend $D$ on the ex-dividend date $t_d$. The stock drops on that date, which raises the put's payoff, so unlike the call, the put benefits from the dividend by holding rather than by exercising. Holding keeps the optionality and gains from the drop, while exercising collects the interest on $K$ until the ex-dividend date. If the dividend alone exceeds the interest, holding wins regardless of the optionality:

$$D \geq K\left(e^{r(t_d - t)} - 1\right)$$

Over a short interval this is approximately $t_d - t \leq D/(rK)$. This means close to the ex-dividend date the put should not be exercised. With a longer time to the ex-dividend date, the interest grows past the dividend, and exercise can be optimal. There the optionality has to be weighed alongside the dividend against the interest, which is the same trade-off as in the no-dividend case, now with an additional dividend term on the holding side.

The figure below shows how dividends change the early exercise boundary, for a sample contract with $K = 100$, $\sigma = 30\%$, $r = 5\%$ and $T = 1$, paying $0.5$ each quarter.

{{< dividend-boundary >}}

The shape looks very different from the no-dividend boundary, which rises steadily. This is because there are now two factors at play. Without dividends the only clock is expiry: the optionality decays as it approaches, so exercising becomes increasingly attractive, and the boundary rises. A dividend adds a second clock, running toward the ex-dividend date. As that date comes closer, exercising becomes less attractive, and the boundary falls. Near each ex-dividend date the second clock dominates the first, which is why the boundary runs the other way and eventually disappears: no stock price is low enough to justify exercising with the drop coming soon. So the dividend postpones exercise rather than ruling it out.

[^tax]: This is commonly attributed to the differing tax treatment of dividends and capital gains.

## Summary

Every case is the same comparison: what holding gains against what exercising gains.

| Call | Holding gains | Exercising gains | Optimal to exercise early? |
|---|---|---|---|
| $r > 0$, no dividends | optionality, and the interest on $K$ while payment is deferred | nothing | **No**, worth the same as European |
| $r = 0$ | optionality | nothing | **No** |
| $r < 0$ | optionality | the interest on $K$, which is now a cost to defer | **Yes**, once $S$ rises to the call's exercise boundary, which sits above $K$ |
| Dividends, $r > 0$ | optionality, and the interest on $K$ | the dividend $D$ | **Yes**, but only just before an ex-dividend date, and only if $D$ is large enough |

<br>

| Put | Holding gains | Exercising gains | Optimal to exercise early? |
|---|---|---|---|
| $r > 0$, no dividends | optionality | the interest on $K$, received today rather than at expiry | **Yes**, once $S$ falls to the exercise boundary, worth more than European |
| $r = 0$ | optionality | nothing | **No**, worth the same as European |
| $r < 0$ | optionality, and $K$ is worth more received later than today | nothing | **No**, worth the same as European |
| Dividends, $r > 0$ | optionality, and the dividend $D$ through the drop in $S$ | the interest on $K$ until the ex-dividend date | **Yes**, but not when less than about $D/(rK)$ remains until an ex-dividend date |

**Bonus question.** What happens when $r < 0$ and the stock pays a dividend? Follow the same trade-off analysis and work out whether the interest rate and the dividend reinforce each other or work against each other. I will leave that one to the readers to figure out.