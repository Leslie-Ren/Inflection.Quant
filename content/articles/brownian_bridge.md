---
title: "The Brownian Bridge: What Brownian Motion Looks Like When You Know the Endpoints"
date: 2026-05-28
draft: false
math: true
tags: ["Brownian Motion", "Brownian Bridge", "Conditional Expectation", "Survival Probability"]
---

## Why This Matters

In my earlier article on [Brownian motion]({{< relref "understanding_brownian_motion.md" >}}), I worked through the forward view: a process starting at a known value, diffusing into an uncertain future. Sometimes we know more than just the starting point. We also know where the process ended up, and we want to characterize the path in between. The object that answers this is the Brownian bridge: a Brownian motion conditioned on its terminal value.

My motivation for writing this is a follow-up article on [Monte Carlo variance reduction]({{< relref "mc_variance_reduction.md" >}}), where the bridge supports two of the techniques: conditional Monte Carlo uses it to accelerate the pricing of barrier options and remove the monitoring bias, and quasi-Monte Carlo uses it to reduce the effective dimension of the simulation for faster convergence. To keep that article focused on variance reduction, I am introducing the Brownian bridge here as a standalone reference. This article is meant as a brief introduction rather than a deep dive: the aim is to define the bridge, construct its form, derive the survival probability behind conditional Monte Carlo, and show how the bridge reduces the effective dimension for quasi-Monte Carlo.

## Definition

Formally, the Brownian bridge from $a$ to $b$ on $[0, T]$ is the conditional process

$$
\{W_t \mid W_0 = a,\ W_T = b\},\quad 0 \leq t \leq T.
$$

The definition is abstract. To work with the bridge in practice we need a concrete formula that, given a sample of an ordinary Brownian motion, produces a sample of the bridge.

## Decomposing $W_t$
 
To get a concrete formula we split $W_t$ into a piece that is predictable from $W_T$ and a piece that is independent of it. The familiar version of this idea, for two standard normals $X$ and $Y$ with correlation $\rho$, is
 
$$
Y = \rho X + \sqrt{1 - \rho^2}\, Z,
$$
 
where $Z$ is a standard normal independent of $X$. The first term is the predictable piece (a multiple of $X$); the second is the independent residual.
 
We want the same kind of split for $W_t$ in terms of $W_T$, but $W_t$ and $W_T$ have variances $t$ and $T$, not 1. The role of $\rho X$ in the familiar form is to be the best linear predictor of $Y$ given $X$. For variables with arbitrary variance, the best linear predictor of $W_t$ given $W_T$ is $\text{Cov}(W_t, W_T) / \text{Var}(W_T) \cdot W_T$[^1]:
 
$$
\mathbb{E}[W_t \mid W_T] = \frac{\text{Cov}(W_t, W_T)}{\text{Var}(W_T)} W_T = \frac{t}{T} W_T.
$$
 
For Brownian motion $\text{Cov}(W_t, W_T) = t$ and $\text{Var}(W_T) = T$, giving the slope $t/T$.
 
[^1]: This is the least-squares regression slope. Differentiating $\mathbb{E}[(Y - \beta X)^2]$ in $\beta$ and setting to zero gives $\beta = \mathbb{E}[XY] / \mathbb{E}[X^2]$, which for zero-mean variables is $\text{Cov}(X, Y) / \text{Var}(X)$.
 
Subtracting the projection from $W_t$ leaves the residual
 
$$
R_t = W_t - \frac{t}{T} W_T,
$$
 
which is independent of $W_T$ (zero covariance, both Gaussian). So we have the decomposition
 
$$
W_t = \underbrace{\frac{t}{T} W_T}_{\text{predictable from } W_T} \;+\; \underbrace{R_t}_{\text{independent of } W_T}.
$$

## Conditioning on $W_T = b$

With the split in hand, conditioning on the endpoint is straightforward, because it acts on each piece separately. The first piece is a function of $W_T$, so fixing $W_T = b$ turns it into the deterministic value $\frac{t}{T} b$. The second piece is independent of $W_T$, so its distribution is unchanged. Putting them back together:

$$
\{W_t \mid W_T = b\} = \frac{t}{T} b + \left(W_t - \frac{t}{T} W_T\right).
$$

This is a bridge from $0$ to $b$. To start at $a$ instead, add the deterministic term $a(1 - t/T)$, which equals $a$ at $t = 0$ and zero at $t = T$:

$$
B_t = a + \frac{t}{T}(b - a) + \left(W_t - \frac{t}{T} W_T\right),\quad 0 \leq t \leq T.
$$

To simulate $B_t$, note that the $W_T$ inside the bracket is the random terminal value of the unconditioned Brownian motion we started with, not the fixed endpoint $b$. The formula is a recipe: simulate an ordinary Brownian motion, read off $W_t$ and $W_T$, and apply the transformation. The two are correlated, since $W_T$ contains everything $W_t$ did, so they are generated in order from two independent standard normals $Z_1$ and $Z_2$:

$$
W_t = \sqrt{t}\,Z_1, \qquad W_T = W_t + \sqrt{T - t}\,Z_2.
$$

## Mean and Covariance

Before putting the bridge to use, it helps to know what a typical bridge path looks like. The construction expresses $B_t$ as a linear combination of Gaussian variables, so $B_t$ is itself Gaussian at every time, and the joint distribution of $(B_{t_1}, \ldots, B_{t_n})$ at any finite set of times is multivariate Gaussian. The process is fully characterized by its mean and covariance:

$$
\mathbb{E}[B_t] = a + \frac{t}{T}(b - a),
$$

$$
\text{Cov}(B_s, B_t) = \frac{s(T - t)}{T},\quad 0 \leq s \leq t \leq T,
$$

with variance $\text{Var}(B_t) = t(T - t)/T$, a parabola vanishing at both endpoints and peaking at $T/2$. The uncertainty starts at zero, grows toward the middle of the interval, then shrinks back to zero as the process approaches the fixed endpoint. The bridge is therefore "pulled" toward the terminal value as maturity approaches.

Knowing that $B_t$ is Gaussian with this mean and variance also gives a more direct way to simulate it than the recipe in the previous section. Instead of simulating $W_t$ and $W_T$ and transforming them, a single standard normal $Z$ is enough:

$$
B_t = a + \frac{t}{T}(b - a) + \sqrt{\frac{t(T - t)}{T}}\;Z.
$$

This samples the bridge at one time. Points at several times are correlated, so they cannot be drawn independently this way; the last section shows how to build a whole path.

## Survival Probability via the Reflection Principle
 
Now we can turn to the question that makes the bridge so useful in pricing: given where a path started and where it ended, how likely is it that it never touched a barrier in between? This survival probability has a closed-form expression. This is the foundation of analytical down-and-out barrier option pricing and of the [conditional Monte Carlo]({{< ref "mc_variance_reduction.md#technique-3-conditional-monte-carlo" >}}) technique I cover in the variance reduction article.
 
The bridge is the object whose survival probability we want, but the derivation does not use the bridge construction directly. The conditional probability is computed from unconditional Brownian motion via Bayes' rule, with the reflection principle supplying the joint density of the minimum and the terminal value.
 
I will derive the simple case here: standard Brownian motion starting at $W_0 = 0$ with unit volatility, conditioned on $W_T = b$, with a lower barrier $L < 0$ and $b > L$ (since the path starts at $0$, the lower barrier must be negative). The question is
 
$$
P\!\left(\min_{0 \leq s \leq T} W_s > L \;\middle|\; W_T = b\right) = ?
$$
 
### The Reflection Principle
 
The tool we need is the reflection principle for standard Brownian motion. For any barrier $L < 0$ and any terminal value $b > L$,
 
$$
P\!\left(\min_{0 \leq s \leq T} W_s \leq L,\ W_T = b\right) = P(W_T = 2L - b).
$$
 
This is a density identity: $P(W_T = b)$ is shorthand for the density $\phi_T(b)$, and the equation should be read as "the density of $W_T$ at $b$ on the event that the path crossed $L$ equals the unconditional density of $W_T$ at $2L - b$." The idea behind it: every path that touches $L$ before time $T$ and ends at $b$ can be reflected about $L$ from its first hitting time onward, producing a path that ends at $2L - b$. The reflection is a bijection between the two sets of paths, and Brownian motion's symmetry means it preserves the density.
 
{{< reflection_principle >}}
 
### From Reflection to Survival
 
Reading the reflection identity as a density statement: the density of $W_T$ at $b$ on the event $\{m_T \leq L\}$ is $\phi_T(2L - b)$, where $m_T = \min_{s \leq T} W_s$ and $\phi_T$ is the density of $W_T \sim \mathcal{N}(0, T)$. Dividing by the unconditional density $\phi_T(b)$ gives the conditional probability of having hit the barrier given the terminal value:
 
$$
P\!\left(m_T \leq L \;\middle|\; W_T = b\right) = \frac{\phi_T(2L - b)}{\phi_T(b)} = \exp\!\left(-\frac{2L(L - b)}{T}\right).
$$
 
The survival probability is the complement:
 
$$
\boxed{\;P\!\left(\min_{0 \leq s \leq T} W_s > L \;\middle|\; W_T = b\right) = 1 - \exp\!\left(-\frac{2L(L - b)}{T}\right).\;}
$$

The survival probability depends on how far each endpoint sits from the barrier. If either endpoint approaches the barrier, the survival probability falls toward zero. If both endpoints move far away from the barrier, the exponent becomes large and negative, pushing the survival probability toward one.
  
### General Case
 
For Brownian motion starting at $W_0 = a$ with volatility $\sigma$, conditioned on $W_T = b$, with lower barrier $L$ and $a, b > L$, the general result follows from the simple case by a change of variables. Define
 
$$
\tilde{W}_s = \frac{W_s - a}{\sigma}.
$$
 
Then $\tilde{W}$ is standard Brownian motion starting at $0$. The event $\{W_s \leq L\}$ becomes $\{\tilde{W}_s \leq (L - a)/\sigma\}$, and the conditioning $\{W_T = b\}$ becomes $\{\tilde{W}_T = (b - a)/\sigma\}$. Substituting into the simple-case formula with barrier $(L - a)/\sigma$ and terminal value $(b - a)/\sigma$:
 
$$
P\!\left(\min_{0 \leq s \leq T} W_s > L \;\middle|\; W_0 = a,\ W_T = b\right) = 1 - \exp\!\left(-\frac{2(a - L)(b - L)}{\sigma^2 T}\right).
$$
 
The symmetric form in $(a - L)$ and $(b - L)$ reflects the fact that the bridge is time-reversible: running the bridge backwards in time (mapping $s$ to $T - s$) gives another Brownian bridge with the start and end swapped. Reversing time does not change which values the path visits, only the order, so the minimum is the same in both directions. The survival probability must therefore give the same answer whether we call the endpoints $(a, b)$ or $(b, a)$, which forces the formula to be symmetric in the two.

## Building a Path from the Endpoint Inward

So far we have used the bridge to answer questions about a Brownian path whose endpoints are given. The same conditional distribution also lets us simulate Brownian motion paths in a different order from the usual one. The usual way builds a path forward in time, one increment after another. The bridge lets us fix the endpoint first, then fill in points between values we already know.

Suppose $W_u$ and $W_v$ are known at times $u < v$, and we want $W_s$ at some $u < s < v$. Given its two ends, the segment in between is a Brownian bridge, so from the mean and covariance above, $W_s$ is normal with

$$
\mathbb{E}[W_s \mid W_u, W_v] = \frac{(v-s)\,W_u + (s-u)\,W_v}{v-u}, \qquad \text{Var}(W_s \mid W_u, W_v) = \frac{(s-u)(v-s)}{v-u}.
$$

A path on a grid can then be built by bisection. Draw the endpoint $W_T = \sqrt{T}\,Z_1$, then the midpoint given both ends using $Z_2$, then the two quarter points using $Z_3$ and $Z_4$, and so on until every grid point is filled. Each new point consumes one standard normal, so the path uses the same number of normals as the forward construction and has exactly the same distribution.

What changes is how much each normal matters. Each new point is placed in the middle of a gap between two points we already know, and its uncertainty depends only on the size of that gap. Setting $s$ halfway between $u$ and $v$ in the conditional variance above shows that a point in the middle of a gap of length $g$ has variance $g/4$. The first normal sets the endpoint, with the full uncertainty of a one-year Brownian motion. The next sets the midpoint, which is already pinned down from both sides, so its variance is a quarter of the endpoint's. Every round halves the gaps, and the new points become steadily more constrained. With $T = 1$:

| Round | Normals used | Points filled | Gap between known points | Variance of each new point |
|---|---|---|---|---|
| 0 | $Z_1$ | endpoint, $t = 1$ | none | 1 |
| 1 | $Z_2$ | midpoint, $t = \tfrac{1}{2}$ | 1 | $\tfrac{1}{4}$ |
| 2 | $Z_3, Z_4$ | quarters, $t = \tfrac{1}{4}, \tfrac{3}{4}$ | $\tfrac{1}{2}$ | $\tfrac{1}{8}$ |
| 3 | $Z_5$ to $Z_8$ | eighths | $\tfrac{1}{4}$ | $\tfrac{1}{16}$ |
| … | … | … | … | … |
| 8 | the last ~120 | each step of a 252-step grid | $\tfrac{1}{128}$ | $\tfrac{1}{512}$ |
| … | … | … | … | … |

The first few normals therefore decide the overall shape of the path: where it ends, where it sits halfway, and so on. The later ones only add small wiggles between neighbors that are already fixed. For ordinary Monte Carlo this ordering makes no difference, since every normal is equally random. For quasi-Monte Carlo it matters a great deal, because the first few inputs of a low-discrepancy sequence are its most evenly spread. The [Monte Carlo variance reduction article]({{< ref "mc_variance_reduction.md#reducing-the-effective-dimension" >}}) shows how this reduces the effective dimension when pricing a barrier option.

<details style="margin-bottom: 1.5em;">
<summary>Brownian bridge path construction</summary>

```python
def bridge_order(M):
    # bisection order: endpoint first, then midpoints of each known interval
    order, left, right = [M - 1], [-1], [-1]
    queue = [(-1, M - 1)]
    while queue:
        nxt = []
        for a, b in queue:
            if b - a > 1:
                m = (a + b) // 2
                order.append(m); left.append(a); right.append(b)
                nxt += [(a, m), (m, b)]
        queue = nxt
    return order, left, right

def bridge_paths(Z, T, M):
    t = np.arange(1, M + 1) * T / M
    order, left, right = bridge_order(M)
    W = np.zeros_like(Z)
    W[:, M - 1] = np.sqrt(T) * Z[:, 0]
    for k in range(1, M):
        i, a, b = order[k], left[k], right[k]
        ta = 0.0 if a < 0 else t[a]
        Wa = 0.0 if a < 0 else W[:, a]
        tb, ti = t[b], t[i]
        mean = ((tb - ti) * Wa + (ti - ta) * W[:, b]) / (tb - ta)
        W[:, i] = mean + np.sqrt((ti - ta) * (tb - ti) / (tb - ta)) * Z[:, k]
    return W, t
```

Given an array `Z` of standard normals with one column per time step, `bridge_paths(Z, T, M)` returns the Brownian paths `W` on the grid `t`. A geometric Brownian motion path then follows as `S0 * np.exp((r - 0.5*sigma**2)*t + sigma*W)`.

</details>
