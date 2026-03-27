---
title: "Understanding Brownian Motion: Why Variance Grows Linearly in Time"
date: 2026-03-26
draft: false
math: true
tags: ["Quant", "Brownian Motion"]
---

Brownian motion — the mathematical model behind everything from stock prices to the way heat spreads — has a surprisingly elegant property: the variance of a particle’s position at time $t$ grows **linearly** with time. Not $t^2$, not $\sqrt{t}$, but exactly $t$. At first glance, this might seem like a purely mathematical choice, but it actually comes from both physical observations and mathematical construction.

---

## Physical Motivation

The story starts with Robert Brown, who in 1827 observed pollen grains jittering randomly in water. The motion seemed erratic, almost chaotic. Decades later, Einstein (1905) made this quantitative, showing that a particle’s **mean squared displacement** (MSD) grows in proportion to time:

$$
\mathbb{E}[(X(t) - X(0))^2] \sim t
$$

What does this mean in plain terms? A few key points:

* At each instant, the particle moves in a random direction.  
* Over longer periods, the spread of particle positions increases — but importantly, it grows **linearly** with elapsed time, not faster.  

This linear behavior is a hallmark of diffusive motion and one of the first hints that something simple underlies the apparent randomness.

---

## Definition of Standard Brownian Motion

Mathematically, a **standard Brownian motion** $B(t)$ is a continuous-time stochastic process defined by:

1. $B(0) = 0$  
2. Independent increments: Each increment $B(t+s) - B(s)$ is independent of the past.  
3. Normally distributed increments: $B(t+s) - B(s) \sim N(0, t)$, meaning the variance equals the length of the time interval.  

Notice how the last condition encodes **linear variance growth**: an increment over a time interval of length $\Delta t$ has variance exactly $\Delta t$. But why this scaling, and not something else? To answer that, we can think of Brownian motion as the limit of a simple random walk.

---

## Brownian Motion as a Limit of a Random Walk

Imagine a simple symmetric random walk:

$$
S_n = X_1 + X_2 + \dots + X_n
$$

where each $X_i$ is either $+1$ or $-1$ with equal probability. Because the steps are independent:

$$
\text{Var}(S_n) = \text{Var}(X_1) + \dots + \text{Var}(X_n) = n
$$

Notice how the variance grows linearly with the **number of steps** — already a hint of the behavior we want to capture in continuous time.

---

#### Scaling Step Size to Match Time

To stretch this random walk over a fixed time horizon $t$, we divide time into $n$ small intervals:

$$
\Delta t = \frac{t}{n}
$$

Now, the trick is to pick a step size $\delta$ for each increment so that the total variance matches the linear growth observed in real particles. Setting:

$$
\delta = \sqrt{\Delta t} = \sqrt{\frac{t}{n}}
$$

gives:

$$
B(t) \approx \delta (X_1 + X_2 + \dots + X_n)
$$

and the total variance becomes:

$$
\text{Var}[B(t)] = n \cdot (\delta)^2 = n \cdot \frac{t}{n} = t
$$

In other words, the **$\sqrt{\Delta{t}}$** scaling is exactly what’s needed to reproduce the observed linear growth in variance.

---

#### Why Other Step Sizes Don’t Work

What happens if we pick a different scaling for the steps?  

| Step size | Resulting variance | Problem |
|-----------|-----------------|---------|
| Constant $c$ | $c^2 n$ | Variance blows up as $n \to \infty$; the process has no well-defined limit. |
| $1/n$ | $(1/n)^2 n = 1/n$ | Variance → 0; the process becomes almost deterministic. |
| $t/\sqrt{n}$ | $(t/\sqrt{n})^2 n = t^2$ | Quadratic growth; the particle spreads faster than what experiments show. |

Only the **$\sqrt{\Delta{t}}$** choice gives variance proportional to $t$, consistent with what we actually observe in physical diffusion.

---

## Conclusion

The variance of Brownian motion grows linearly with time because it reflects the diffusive spread of real particles, and the $\sqrt{\Delta{t}}$ scaling in the random-walk limit is what makes this property hold mathematically — a simple rule hiding behind apparently random motion.