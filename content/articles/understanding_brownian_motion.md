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

Brownian motion was first observed by Robert Brown in 1827 as the erratic motion of microscopic particles in water. Einstein (1905) quantified it, showing that a particle's mean squared displacement grows proportionally to time: 

$$ \mathbb{E}[(X(t) - X(0))^2] \propto t $$

Einstein's key insight was that this linear scaling follows from particles receiving many small, independent random kicks. Because the kicks are independent, their variances add, so the total variance grows linearly with the number of kicks, and hence with time. The Central Limit Theorem also tells us that the sum of many small independent kicks is approximately Gaussian, which is why the displacement of a Brownian particle is normally distributed. This linear growth in variance is the hallmark of diffusive motion, and it is the same additivity argument that drives the random-walk construction below.

---

## Definition of Standard Brownian Motion
A standard Brownian motion $B(t)$ is a stochastic process indexed by $t \in [0, \infty)$ satisfying:

1. Starts at the origin: $B(0) = 0$.
2. Continuous paths: With probability 1, the map $t \mapsto B(t)$ is
   continuous.
3. Independent increments: For any $0 \le t_0 < t_1 < \cdots < t_n$, the
   increments $B(t_1) - B(t_0),\ \ldots,\ B(t_n) - B(t_{n-1})$ are mutually
   independent.
4. Normally distributed increments: For $0 \le s < t$,
   $$B(t) - B(s) \sim N(0,\, t - s),$$

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

A concrete place where the $\sqrt{T}$ scaling shows up in practice is in the price of
at-the-money options. To isolate the effect of Brownian motion cleanly, we work under two
idealising assumptions: the option is struck exactly at the current price ($S = K$), and the
risk-free rate is zero ($r = 0$). In practice, implied volatility varies across expiries and
rates are non-zero, both of which distort the pure $\sqrt{T}$ relationship. Stripping them
away lets the Brownian motion signature come through clearly.

Under these assumptions, consider a European call on a stock whose *log-price* is Brownian
with variance $\sigma^2 T$. The option's value should grow with time-to-expiry $T$, since more time means
more opportunity for the stock to move in your favour. But by how much? The premium tracks the
*standard deviation* of the terminal price rather than its variance, because an ATM call is
compensating the seller for an expected absolute deviation. If variance grew as $T^2$, the standard deviation would grow as $T$ and
prices would scale linearly with expiry; if variance were constant, the standard deviation
would be too, and prices would not change with expiry at all. Because Brownian motion gives
$\text{Var}(B_T) = T$, the proportional dispersion of the terminal price scales as
$\sigma\sqrt{T}$, and the premium scales the same way.

The table below shows Black-Scholes ATM call prices for $S = K = 100$, $\sigma = 20\%$,
$r = 0$:

| Expiry $T$ | Call Price | Ratio to 1-month price |
|---|---|---|
| 1 month | 2.303 | 1.000 |
| 4 months | 4.604 | 1.999 |
| 9 months | 6.901 | 2.997 |
| 16 months | 9.193 | 3.992 |

The expiry quadruples from 1 to 4 months, yet the price only doubles. That is the $\sqrt{T}$
signature of Brownian motion. The ratios track $\sqrt{T}$ closely,
and the small shortfall reflects curvature in the normal CDF, explained below.

### Why the ATM Black-Scholes Formula Reduces to a Function of $\sigma\sqrt{T}$

Under the two conditions $S = K$ and $r = 0$, the Black-Scholes $d_1$ and $d_2$ terms simplify
considerably. Recall:

$$
d_1 = \frac{\ln(S/K) + (r + \frac{1}{2}\sigma^2)T}{\sigma\sqrt{T}}, \qquad d_2 = d_1 - \sigma\sqrt{T}
$$

Setting $S = K$ and $r = 0$:

$$
d_1 = \frac{\frac{1}{2}\sigma^2 T}{\sigma\sqrt{T}} = \frac{\sigma\sqrt{T}}{2}, \qquad d_2 = -\frac{\sigma\sqrt{T}}{2}
$$

The call price formula $C = S \cdot N(d_1) - K \cdot e^{-rT} N(d_2)$ then becomes:

$$
C = S \left( N\left(\tfrac{\sigma\sqrt{T}}{2}\right) - N\left(-\tfrac{\sigma\sqrt{T}}{2}\right) \right) = S \left( 2N\left(\tfrac{\sigma\sqrt{T}}{2}\right) - 1 \right)
$$

This is exact, and it is the formula behind the table above. Note that the price depends on
$\sigma$ and $T$ only through the combination $\sigma\sqrt{T}$, with no separate dependence on
either.

To make the $\sqrt{T}$ dependence fully explicit, expand the normal CDF about the origin.
Since $N'= \phi$ and $\phi(0) = \frac{1}{\sqrt{2\pi}}$, for small arguments

$$
2N\left(\tfrac{\sigma\sqrt{T}}{2}\right) - 1 = \sqrt{\frac{2}{\pi}}\left(\frac{\sigma\sqrt{T}}{2} - \frac{1}{6}\left(\frac{\sigma\sqrt{T}}{2}\right)^{3} + O\!\left((\sigma\sqrt{T})^{5}\right)\right)
$$

Keeping only the leading term:

$$
C \approx S \cdot \sqrt{\frac{2}{\pi}} \cdot \frac{\sigma\sqrt{T}}{2} = \frac{S\sigma\sqrt{T}}{\sqrt{2\pi}} \approx 0.4\, S \sigma \sqrt{T}
$$

which is the familiar desk rule of thumb. The cubic term is always negative, so the linear
approximation always *overstates* the true premium, by a relative amount of roughly
$\sigma^2 T / 24$. At 20% vol that is under 0.1% inside six months and about 0.2% at the
16-month point in the table. The same curvature is what pulls the observed ratios slightly
below 2.000, 3.000 and 4.000. Any model that replaced Brownian motion with a process whose variance scaled differently would
produce option prices inconsistent with this pattern, which is one reason the continuous-time
Brownian framework is so deeply embedded in derivatives pricing.

### Connection to Brenner–Subrahmanyam

The relationship between ATM call price and $\sigma\sqrt{T}$ is exactly what motivates the
Brenner–Subrahmanyam initial guess for implied volatility:

$$
\sigma_0 \approx \frac{C_\text{mkt}}{S} \cdot \sqrt{\frac{2\pi}{T}}
$$

Inverting the leading-order relationship gives a smart starting point for **Newton-Raphson**
iterations when solving for implied volatility. For a deeper dive into solving for implied vol, including
**Newton vs Brent methods**, see my article:
[Solving for Implied Volatility: Newton's Method vs Brent's Method]({{< relref "newton_vs_brent_vol_solver.md" >}}).

---

## Conclusion

The variance of Brownian motion grows linearly with time because it models the diffusive behaviour of real particles, and the $\sqrt{\Delta t}$ scaling in the random-walk limit ensures this property holds in continuous time. This is not merely a mathematical convenience — it is the property that ties the abstract construction directly to observable phenomena, from the spread of particle positions to the pricing of financial derivatives. The $\sqrt{T}$ signature appears wherever diffusion governs the dynamics, and recognising it is one of the more transferable intuitions in quantitative finance.

---

