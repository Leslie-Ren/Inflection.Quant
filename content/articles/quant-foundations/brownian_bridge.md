---
title: "The Brownian Bridge: What Brownian Motion Looks Like When You Know the Endpoints"
date: 2026-05-28
draft: false
math: true
tags: ["Brownian Motion", "Brownian Bridge", "Conditional Expectation", "Survival Probability"]
---

## Why This Matters

In my earlier article on [Brownian motion]({{< relref "understanding_brownian_motion.md" >}}), I worked through the forward view: a process starting at a known value, diffusing into an uncertain future. Sometimes we know more than just the starting point. We also know where the process ended up, and we want to characterise the path in between. The object that answers this is the Brownian bridge: a Brownian motion conditioned on its terminal value.

My motivation for writing this is a follow-up article on Monte Carlo variance reduction, where one of the techniques, conditional Monte Carlo, uses the bridge to accelerate the pricing of barrier options. To keep that article focused on variance reduction, I am introducing the Brownian bridge here as a standalone reference. This article is meant as a brief introduction rather than a deep dive: the aim is to define the bridge, construct its form, and derive the survival probability that sits at the core of the variance reduction technique.

---

## Definition

Formally, the Brownian bridge from $a$ to $b$ on $[0, T]$ is the conditional process

$$
\{W_t \mid W_0 = a,\ W_T = b\},\quad 0 \leq t \leq T.
$$

The definition is abstract. To work with the bridge in practice we need a concrete formula that, given a sample of an ordinary Brownian motion, produces a sample of the bridge.

---

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

Conditioning acts on each piece separately. The first piece is a function of $W_T$, so fixing $W_T = b$ turns it into the deterministic value $\frac{t}{T} b$. The second piece is independent of $W_T$, so its distribution is unchanged. Putting them back together:

$$
\{W_t \mid W_T = b\} = \frac{t}{T} b + \left(W_t - \frac{t}{T} W_T\right).
$$

This is a bridge from $0$ to $b$. To start at $a$ instead, add the deterministic term $a(1 - t/T)$, which equals $a$ at $t = 0$ and zero at $t = T$:

$$
B_t = a + \frac{t}{T}(b - a) + \left(W_t - \frac{t}{T} W_T\right),\quad 0 \leq t \leq T.
$$

### A Word on Notation

The $W_T$ inside the bracket on the right is the random terminal value of the unconditioned Brownian motion we started with. It is not equal to $b$. The formula is a recipe: take any sample path of an unconditioned Brownian motion, look at its terminal value $W_T$, and apply the transformation. The subtraction removes the part of the sampled path that is correlated with its own terminal value; the linear interpolation then puts the desired endpoint $b$ back in. At $t = T$ the bracket vanishes by construction (the path's own terminal value cancels itself) and $B_T = b$.

---

## Mean and Covariance

The construction expresses $B_t$ as a linear combination of Gaussian variables, so $B_t$ is itself Gaussian at every time, and the joint distribution of $(B_{t_1}, \ldots, B_{t_n})$ at any finite set of times is multivariate Gaussian. The process is fully characterised by its mean and covariance:

$$
\mathbb{E}[B_t] = a + \frac{t}{T}(b - a),
$$

$$
\text{Cov}(B_s, B_t) = \frac{s(T - t)}{T},\quad 0 \leq s \leq t \leq T,
$$

with variance $\text{Var}(B_t) = t(T - t)/T$, a parabola vanishing at both endpoints and peaking at $T/2$. The uncertainty starts at zero, grows toward the middle of the interval, then shrinks back to zero as the process approaches the fixed endpoint. The bridge is therefore "pulled" toward the terminal value as maturity approaches.

---
## Survival Probability via the Reflection Principle
 
The survival probability of the bridge, the probability that the bridge stays above a barrier, has a closed-form expression. This is the foundation of analytical down-and-out barrier option pricing and of the conditional Monte Carlo technique I cover in the [variance reduction article]({{< relref "mc_variance_reduction.md" >}}).
 
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
 
 