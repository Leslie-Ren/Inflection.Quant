---
title: "Brownian Motion: From Random Walks to Option Prices"
date: 2026-03-26
draft: false
math: true
tags: ["Brownian Motion","Random Walk", "Option Price", "Brenner–Subrahmanyam approximation"]
---
## Why This Matters

Brownian motion, the mathematical model underlying everything from stock prices to heat diffusion, has one of its most elegant properties: the variance of its position at time $t$ grows linearly with time. Not $t^2$, not $\sqrt{t}$, but exactly $t$. This seemingly abstract fact has a concrete consequence in financial markets: under the idealised conditions of an at-the-money option with zero rates, it is precisely why option prices scale with $\sqrt{T}$ rather than $T$, a direct fingerprint of Brownian motion inside Black-Scholes. Understanding why requires looking at both **physical observations** and the **mathematical construction** of Brownian motion.

---

## Physical Motivation

Brownian motion was first observed by Robert Brown in 1827 as the erratic motion of pollen particles in water. Later, Einstein (1905) quantified it, showing that a particle's mean squared displacement grows proportionally to time:

$$
\mathbb{E}[(X(t) - X(0))^2] \sim t
$$

Einstein's key insight was that this $\sim t$ scaling follows directly from particles receiving many small, independent random kicks — exactly the structure that the Central Limit Theorem later formalises. Because each kick is independent, their variances add, and the total variance accumulates linearly with the number of kicks, and hence with time. This linear growth — where the spread of positions increases proportionally to elapsed time, not faster — is the hallmark of diffusive motion, and it is the same additivity argument that underpins the random-walk construction below.

---

## Definition of Standard Brownian Motion

Mathematically, a standard Brownian motion $B(t)$ is a continuous-time stochastic process defined by:

1. $B(0) = 0$.
2. Independent increments: The increment $B(t+s) - B(s)$ is independent of the past.
3. Normally distributed increments: $B(t+s) - B(s) \sim N(0, t)$, where the variance equals the length of the interval.

Notice that the last condition encodes linear variance growth: an increment over a time interval of length $\Delta t$ has variance $\Delta t$. But why is this the correct scaling rather than an arbitrary choice? The answer comes from modeling Brownian motion as a limit of a random walk.

---

## Brownian Motion as a Limit of a Random Walk

Consider a simple symmetric random walk:

$$
S_n = X_1 + X_2 + \dots + X_n
$$

where each $X_i$ is $\pm 1$ with equal probability. Because the steps are independent:

$$
\text{Var}(S_n) = \text{Var}(X_1) + \dots + \text{Var}(X_n) = n
$$

Variance grows linearly with the number of steps.

### Scaling Step Size to Match Time

To approximate Brownian motion over a fixed time horizon $t$, divide it into $n$ small steps of length:

$$
\Delta t = \frac{t}{n}
$$

Now, assign a step size $\delta$ to each increment such that the total variance matches the observed linear growth in $t$. Let:

$$
\delta = \sqrt{\Delta t} = \sqrt{\frac{t}{n}}
$$

Then, the Brownian motion approximation is:

$$
B(t) \approx \delta (X_1 + X_2 + \dots + X_n)
$$

and the total variance becomes:

$$
\text{Var}[B(t)] = n \cdot \delta^2 = n \cdot \frac{t}{n} = t
$$

This shows that the step size must scale as $\sqrt{\Delta t}$ to reproduce the observed linear growth of variance in continuous time. Moreover, as $n \to \infty$, the Central Limit Theorem guarantees that the normalised sum $\delta(X_1 + \dots + X_n)$ converges in distribution to a Gaussian — which is why the limit is not merely variance-correct, but fully normally distributed, justifying the $N(0, t)$ increments in the formal definition above.

### Why Other Step Sizes Fail

| Step size | Resulting variance | Problem |
|---|---|---|
| Constant $c$ | $c^2 n$ | Grows without bound as $n \to \infty$; no well-defined limit. |
| $1/n$ | $(1/n)^2 n = 1/n$ | Variance → 0; the process becomes deterministic. |
| $t/\sqrt{n}$ | $(t/\sqrt{n})^2 n = t^2$ | Quadratic growth in time; inconsistent with physical observation. |

Only $\sqrt{\Delta t}$ gives variance proportional to $t$, consistent with both physical observation and the continuous-time limit.

---

## Application: Option Price Scaling with Expiry

A concrete place where the $\sqrt{T}$ scaling shows up in practice is in the price of at-the-money options. To isolate the effect of Brownian motion cleanly, we work under two idealising assumptions: the option is struck exactly at the current price ($S = K$), and the risk-free rate is zero ($r = 0$). In practice, implied volatility varies across expiries and rates are non-zero, both of which distort the pure $\sqrt{T}$ relationship — but stripping these away lets the Brownian motion signature come through clearly.

Under these assumptions, consider a European call option on a stock with no drift. Intuitively, the option's value should grow with time-to-expiry $T$ — more time means more opportunity for the stock to move in your favour. But by how much? If variance grew as $T^2$, prices would scale linearly with $T$; if it were constant, prices would not change with expiry at all. Because Brownian motion gives $\text{Var}(B_T) = T$, the standard deviation of the stock's position scales as $\sqrt{T}$, and since an ATM option's value is essentially compensating the seller for the expected absolute deviation of the terminal stock price, the price scales accordingly.

The table below shows Black-Scholes ATM call prices for $S = K = 100$, $\sigma = 20\%$, $r = 0$:

| Expiry $T$ | Call Price | Ratio to 1-month price |
|---|---|---|
| 1 month | 2.31 | 1.00 |
| 4 months | 4.62 | 2.00 |
| 9 months | 6.93 | 3.00 |
| 16 months | 9.24 | 4.00 |

The expiry quadruples from 1 to 4 months, yet the price only doubles — exactly the $\sqrt{T}$ signature of Brownian motion. This is not an approximation artifact.

### Why the ATM Black-Scholes Formula Reduces to $\sigma\sqrt{T}$

Under the two conditions $S = K$ and $r = 0$, the Black-Scholes $d_1$ and $d_2$ terms simplify considerably. Recall:

$$
d_1 = \frac{\ln(S/K) + (r + \frac{1}{2}\sigma^2)T}{\sigma\sqrt{T}}, \qquad d_2 = d_1 - \sigma\sqrt{T}
$$

Setting $S = K$ (so $\ln(S/K) = 0$) and $r = 0$:

$$
d_1 = \frac{\frac{1}{2}\sigma^2 T}{\sigma\sqrt{T}} = \frac{\sigma\sqrt{T}}{2}, \qquad d_2 = -\frac{\sigma\sqrt{T}}{2}
$$

The call price formula $C = S \cdot N(d_1) - K \cdot e^{-rT} N(d_2)$ then becomes:

$$
C = S \left( N\left(\tfrac{\sigma\sqrt{T}}{2}\right) - N\left(-\tfrac{\sigma\sqrt{T}}{2}\right) \right) = S \left( 2N\left(\tfrac{\sigma\sqrt{T}}{2}\right) - 1 \right)
$$

The price is an exact function of $\sigma\sqrt{T}$ alone — there is no separate dependence on $\sigma$ or $T$ individually, only through their product $\sigma\sqrt{T}$. Doubling $T$ is equivalent to doubling $\sigma^2$. To make the $\sqrt{T}$ dependence fully explicit, apply a first-order Taylor expansion of $N(x)$ around $x = 0$. Since $N'(x) = \phi(x)$ and $\phi(0) = \frac{1}{\sqrt{2\pi}}$, we have $N(x) \approx \frac{1}{2} + \frac{1}{\sqrt{2\pi}} x$ for small $x$, and therefore:

$$
2N(x) - 1 \approx \sqrt{\frac{2}{\pi}}x
$$

Substituting $x = \frac{\sigma\sqrt{T}}{2}$ gives the leading-order ATM call price:

$$
C \approx S \cdot \sqrt{\frac{2}{\pi}} \cdot \frac{\sigma\sqrt{T}}{2} = \frac{S\sigma\sqrt{T}}{\sqrt{2\pi}}
$$

with $\frac{\sigma}{\sqrt{2\pi}}$ acting as a simple proportionality constant. The approximation is accurate for $\sigma\sqrt{T} \leq 0.3$ (e.g. 20% vol with under roughly two years to expiry), and deteriorates for long-dated or high-vol options where the cubic correction term in the Taylor expansion becomes material.

Any model that replaced Brownian motion with a process whose variance scaled differently would produce option prices inconsistent with this pattern, which is one reason the continuous-time Brownian framework is so deeply embedded in derivatives pricing.

### Connection to Brenner–Subrahmanyam

The relationship between ATM call price and $\sigma\sqrt{T}$ is exactly what motivates the Brenner–Subrahmanyam initial guess for implied volatility:

$$
\sigma_0 \approx \frac{C_\text{mkt}}{S} \cdot \sqrt{\frac{2\pi}{T}}
$$

By taking the observed market price and inverting this relationship, we get a smart starting point for **Newton-Raphson** iterations when solving for implied volatility. For a deeper dive into solving for implied vol, including **Newton vs Brent methods**, see my article: [Solving for Implied Volatility: Newton's Method vs Brent's Method]({{< relref "newton_vs_brent_vol_solver.md" >}})


---

## Conclusion

The variance of Brownian motion grows linearly with time because it models the diffusive behaviour of real particles, and the $\sqrt{\Delta t}$ scaling in the random-walk limit ensures this property holds in continuous time. This is not merely a mathematical convenience — it is the property that ties the abstract construction directly to observable phenomena, from the spread of particle positions to the pricing of financial derivatives. The $\sqrt{T}$ signature appears wherever diffusion governs the dynamics, and recognising it is one of the more transferable intuitions in quantitative finance.

---

