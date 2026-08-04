---
title: "Monte Carlo Variance Reduction: What We Average, and How We Sample"
date: 2026-06-02
draft: false
math: true
tags: ["monte-carlo", "variance-reduction", "options", "numerical-methods"]
---

## Why This Matters

In the article on the [Feynman-Kac theorem]({{< relref "feynman_kac.md" >}}), we saw that the price of a derivative can be expressed equivalently as the solution to a deterministic PDE or as the expectation of a discounted payoff under the risk-neutral measure. This gives us two complementary numerical approaches to pricing. For low-dimensional problems with smooth payoffs, finite difference methods on the PDE side are efficient and accurate. For high-dimensional problems, path-dependent payoffs, or models where the PDE is hard to derive, Monte Carlo (MC) on the expectation side becomes the natural choice.

Anyone who has implemented MC in a production pricer notices a challenge immediately: it converges slowly. The standard error of the estimator scales in proportion to $1/\sqrt{N}$, so gaining one extra digit of precision requires a hundredfold increase in computational cost. For a risk system that needs to revalue thousands of trades overnight, or for an intraday hedging system that revalues frequently as the market moves, the computation time adds up quickly.

In this article we walk through five common variance reduction techniques used in practice. Before getting to them, we take a quick look at where the variance comes from and how it shapes the way we approach the problem.

---

## Where the Variance Comes From

Consider the MC estimator for a European-style payoff:

$$\hat{V}_N = \frac{e^{-rT}}{N}\sum_{i=1}^N g(S_T^{(i)})$$

where $g$ is the payoff function and $S_T^{(i)}$ are independent samples of the terminal price under the risk-neutral measure. The variance of this estimator is:

$$\text{Var}(\hat{V}_N) = \frac{e^{-2rT}\,\text{Var}(g(S_T))}{N}$$

The expression points to two distinct levers we can pull. We can **reshape what we average** by replacing $g(S_T)$ with a related quantity that has the same expectation but lower variance. Or we can **change how we sample** by drawing the paths in a smarter way, without touching the payoff itself.

Each of the five techniques in this article fits one of these two patterns:

> - **Reshape what we average:** control variates, conditional MC, importance sampling[^1].
> - **Change how we sample:** antithetic variates, quasi-MC.

The techniques look different on the surface and involve different machinery, but each is just a different way of applying one of these two ideas.

[^1]: Importance sampling sits at the boundary between the two: it changes the sampling distribution, but the reason it reduces variance is that the reweighted payoff has lower variance than the original. We place it on the payoff side because that is what drives its effectiveness.

---

## The Running Example: A Continuously-Monitored Down-and-Out Call

We use the same example throughout the article so the variance reduction from each technique can be compared on a like-for-like basis. The example is a continuously-monitored down-and-out call on a single underlying. The payoff is:

$$g(S) = \max(S_T - K, 0) \cdot \mathbb{1}\{\min_{0 \leq t \leq T} S_t > B\}$$

The option pays the standard call payoff at expiry if and only if the underlying stays strictly above the barrier $B$ at every instant of the path. In practice, barriers can be monitored continuously or at discrete observation dates such as daily or weekly fixings, and the choice is a contract specification. Most of the techniques in this article apply to both; conditional MC is the exception and is best suited to continuous monitoring, as discussed in its section.

The parameters for the example:

| Parameter | Value |
|---|---|
| Spot $S_0$ | 100 |
| Strike $K$ | 100 |
| Barrier $B$ | 85 |
| Volatility $\sigma$ | 30% |
| Risk-free rate $r$ | 5% |
| Maturity $T$ | 1 year |
| Time grid $M$ | 252 steps |
| Number of paths $N$ | 100,000 |

The barrier sits 15% below spot, close enough that a meaningful fraction of paths knock out but far enough that a meaningful fraction survive. This is the regime where the option payoff has meaningful variance and where variance reduction techniques can make a real difference.

Two caveats on the example. In practice, barrier options are typically priced under local vol or stochastic vol models rather than constant-vol GBM; the simplification here is to keep the focus on variance reduction. Under this simplified GBM the barrier option also has a closed-form price (around **11.87** for our parameters), which lets us check our Monte Carlo estimates against the exact value.

### Naive Monte Carlo Baseline

The baseline simulates $N$ paths of the underlying under Geometric Brownian motion on a discrete grid of $M$ steps, checks the barrier at each grid point, and averages the surviving payoffs.

```python
import numpy as np

def simulate_paths(S0, sigma, r, T, M, N):
    dt = T / M
    Z = np.random.standard_normal((N, M))
    increments = (r - 0.5 * sigma**2) * dt + sigma * np.sqrt(dt) * Z
    log_paths = np.log(S0) + np.cumsum(increments, axis=1)
    return np.exp(log_paths)  # shape (N, M)

def naive_mc(S0, K, B, sigma, r, T, M, N):
    paths = simulate_paths(S0, sigma, r, T, M, N)
    survived = np.min(paths, axis=1) > B
    payoffs = np.maximum(paths[:, -1] - K, 0) * survived
    price = np.exp(-r * T) * np.mean(payoffs)
    se = np.exp(-r * T) * np.std(payoffs) / np.sqrt(N)
    return price, se
```

Note that the naive simulation only checks the barrier at the $M$ grid points, but the contract specifies continuous monitoring. The estimator therefore carries a small upward bias from missing between-grid breaches: on our parameters the naive estimator converges to roughly **12.19** rather than the true value of **11.87**, an overprice of about 0.32. This article focuses on variance reduction so we set this bias aside, but we will see that one of our techniques removes the bias directly.

---

## Technique 1: Antithetic Variates

Where does the variance in the MC estimator come from? Each path is an independent realisation of the underlying stochastic process, and thus produces a random payoff. The estimator averages these independent draws, and the variance comes from the dispersion of the payoff across paths. If we could pair each path with a second one whose payoff tends to land on the opposite side of the mean, the dispersion within each pair would partially cancel and the estimator would be tighter.

The simplest way to construct such a pair is to negate the random draws. If a path is generated from $Z_1, \ldots, Z_M$, we generate a partner from $-Z_1, \ldots, -Z_M$. Both draws are valid samples from the standard normal, so each path on its own is a legitimate MC path. The two paths are mirror images of each other: when the original drifts up, the partner drifts down.

Does this actually reduce variance? Let $g(S^+)$ and $g(S^-)$ be the payoffs on the original and negated paths. For each pair we compute the average:

$$\bar{g}^{(i)} = \frac{g(S^+) + g(S^-)}{2}$$

The variance of a single pair average is:

$$\text{Var}(\bar{g}) = \frac{1}{4}\left[\text{Var}(g(S^+)) + \text{Var}(g(S^-)) + 2\,\text{Cov}(g(S^+), g(S^-))\right]$$

Since $Z$ and $-Z$ have the same distribution, the paths $S^+$ and $S^-$ are both valid samples of the same GBM dynamics, so $\text{Var}(g(S^+)) = \text{Var}(g(S^-)) = \text{Var}(g(S))$. The expression simplifies to:

$$\text{Var}(\bar{g}) = \frac{1}{2}\,\text{Var}(g(S)) + \frac{1}{2}\,\text{Cov}(g(S^+), g(S^-))$$

Compare this against two independent paths, whose average has variance $\frac{1}{2}\,\text{Var}(g(S))$. The antithetic average beats the independent average whenever the covariance is negative, i.e., when $g(S^+)$ and $g(S^-)$ move in opposite directions. This is exactly what happens when $g$ is **monotonic** in the underlying: a path that ends up high gives a large payoff and its negated partner gives a small one. Monotonically decreasing payoffs work by the same argument.

If $g$ is non-monotonic the argument breaks down. The covariance can be positive, and the technique can actually *increase* variance compared to independent sampling. This is the main caveat: antithetic variates is safe on calls and puts, useful on barrier options where the payoff is mostly monotonic away from the barrier, and potentially harmful on payoffs with strong non-linearities such as straddles or butterflies.

**Implementation.**

```python
def antithetic_mc(S0, K, B, sigma, r, T, M, N):
    dt = T / M
    Z = np.random.standard_normal((N // 2, M))
    Z_full = np.concatenate([Z, -Z], axis=0)
    increments = (r - 0.5 * sigma**2) * dt + sigma * np.sqrt(dt) * Z_full
    paths = np.exp(np.log(S0) + np.cumsum(increments, axis=1))
    survived = np.min(paths, axis=1) > B
    payoffs = np.maximum(paths[:, -1] - K, 0) * survived
    # pair averages
    pair_avg = (payoffs[:N//2] + payoffs[N//2:]) / 2
    price = np.exp(-r * T) * np.mean(pair_avg)
    se = np.exp(-r * T) * np.std(pair_avg) / np.sqrt(N // 2)
    return price, se
```

**Result.** On our example antithetic variates gives a per-path variance reduction of around 1.5x. There is also a wall-clock benefit: because we only need to draw $N/2$ independent normal vectors and negate them to produce $N$ paths, the random number generation cost is roughly halved. The wall-clock efficiency (tighter estimate in less time) comes out around 1.8x. The technique is essentially free to implement and combines well with the others, which makes it a sensible default to layer on top of any other variance reduction approach.

---

## Technique 2: Control Variates

The down-and-out call shares most of its structure with a vanilla European call. They have the same strike, the same maturity, the same underlying. The only difference is the knock-out feature, which kills some paths that would otherwise have paid off. On every path where the barrier is not breached, the two options pay exactly the same amount.

That suggests a natural question: since we know the vanilla call price analytically from Black-Scholes, can we use it to help price the barrier option? The simplest thing to try is a linear approximation:

$$g_{barrier}(S) = \alpha + \beta \, g_{vanilla}(S) + \varepsilon(S)$$

where $\alpha$ and $\beta$ are constants we want to find and $\varepsilon$ is the residual error on each path.

This is exactly the setup of a linear regression of $g_{barrier}$ on $g_{vanilla}$, with $\alpha$ as the intercept and $\beta$ as the slope. The optimal $\beta$ is the standard regression coefficient:

$$\beta = \frac{\text{Cov}(g_{barrier}, g_{vanilla})}{\text{Var}(g_{vanilla})}$$

It measures how much the barrier payoff moves, on average, per unit movement in the vanilla payoff.

Taking expectations of the linear approximation under $\mathbb{Q}$:

$$\mathbb{E}[g_{barrier}] = \alpha + \beta \, \mathbb{E}[g_{vanilla}]$$

We know $\mathbb{E}[g_{vanilla}]$ from Black-Scholes (it is $e^{rT} C_{BS}$, the undiscounted vanilla price). So if we can estimate $\alpha$, we recover the barrier price directly.

The advantage of estimating $\alpha$ rather than $\mathbb{E}[g_{barrier}]$ comes from its variance. The intuition is that $g_{barrier}$ is noisy across paths, but a large portion of that noise is shared with $g_{vanilla}$ (the two payoffs move together on most paths). By subtracting the scaled vanilla payoff $\beta\,g_{vanilla}$, we cancel out the shared noise and are left with the residual that is specific to the barrier feature itself. The more correlated the two payoffs are, the more noise is cancelled. Working out the variance of $g_{barrier} - \beta\,g_{vanilla}$ at the optimal $\beta$ gives:

$$\text{Var}(\alpha) = \text{Var}(g_{barrier})\,(1 - \rho^2)$$

where $\rho$ is the correlation between the barrier and vanilla payoffs. The closer $\rho$ is to 1, the larger the variance reduction.

We estimate $\alpha$ as the sample mean of $g_{barrier} - \beta\,g_{vanilla}$ over $N$ paths:

$$\hat{\alpha} = \frac{1}{N}\sum_{i=1}^N \left[ g_{barrier}(S^{(i)}) - \beta\,g_{vanilla}(S^{(i)}) \right]$$

and then recover the discounted barrier price by adding back $\beta\,\mathbb{E}[g_{vanilla}] = \beta\,e^{rT}C_{BS}$ and discounting:

$$\hat{V}_{CV} = e^{-rT}\left(\hat{\alpha} + \beta\,e^{rT}C_{BS}\right)$$

where $C_{BS}$ is the analytical Black-Scholes price. In practice $\beta$ is estimated from the simulation itself by computing the sample covariance and variance.

**Implementation.**

```python
from scipy.stats import norm

def bs_call(S0, K, sigma, r, T):
    d1 = (np.log(S0/K) + (r + 0.5*sigma**2)*T) / (sigma*np.sqrt(T))
    d2 = d1 - sigma*np.sqrt(T)
    return S0*norm.cdf(d1) - K*np.exp(-r*T)*norm.cdf(d2)

def control_variate_mc(S0, K, B, sigma, r, T, M, N):
    paths = simulate_paths(S0, sigma, r, T, M, N)
    survived = np.min(paths, axis=1) > B
    barrier_payoffs = np.maximum(paths[:, -1] - K, 0) * survived
    vanilla_payoffs = np.maximum(paths[:, -1] - K, 0)
    # estimate beta
    beta = np.cov(barrier_payoffs, vanilla_payoffs)[0, 1] / np.var(vanilla_payoffs)
    # known mean of vanilla payoff under Q
    vanilla_mean = np.exp(r*T) * bs_call(S0, K, sigma, r, T)
    adjusted = barrier_payoffs - beta * (vanilla_payoffs - vanilla_mean)
    price = np.exp(-r * T) * np.mean(adjusted)
    se = np.exp(-r * T) * np.std(adjusted) / np.sqrt(N)
    return price, se
```

**Result.** On our example the correlation between the barrier and vanilla payoffs is high, giving a variance reduction factor of around 10x. This is a substantial improvement over antithetic variates and reflects how much of the barrier option's variance is "explained" by the vanilla call.

**When to reach for it.** Control variates is most powerful when the correlation between the payoff and the control is high. For the barrier option this works well when the barrier is reasonably far from spot, so most paths survive and the two payoffs move together. When the barrier is very close to spot and most paths knock out, the correlation drops and the technique loses much of its power.

---

## Technique 3: Conditional Monte Carlo

The variance of the naive per-path payoff $Y = \max(S_T - K, 0) \cdot \mathbb{1}\{\text{survived}\}$ comes from two distinct sources. Some of the variability comes from $S_T$ varying across paths: the call payoff lands at different places, and the survival probability is different at different endpoints. But there is also variation that exists *even when $S_T$ is held fixed*: two paths with the same terminal value can have very different outcomes if one touches the barrier and the other does not. If we could remove this second source of variation while keeping the first, the estimator would tighten.

A scatter of $Y$ against $S_T$ makes the two sources visible:

{{< conditional_mc_scatter >}}

For each value of $S_T$ above the strike, paths fall into two groups: surviving paths sit on the diagonal $Y = S_T - K$, and knocked-out paths sit at $Y = 0$. The vertical spread at any given $S_T$ is the within-$S_T$ variance: variation in $Y$ that has nothing to do with where $S_T$ landed. The orange curve is the conditional mean $\mathbb{E}[Y \mid S_T]$, which passes between the two groups in proportion to the survival probability. The slope and position of the orange curve as $S_T$ varies is the between-$S_T$ variance.

The law of total variance is the formal version of this split. For any random payoff $Y$ and any random variable $X$,

$$\text{Var}(Y) = \mathbb{E}[\text{Var}(Y \mid X)] + \text{Var}(\mathbb{E}[Y \mid X]).$$

The first term is the average within-$X$ variance (the vertical spread at each $S_T$ in the picture). The second term is the variance of the conditional mean (how much the orange curve moves as $S_T$ varies). Together they sum to the full variance.

If we can compute $\mathbb{E}[Y \mid X]$ in closed form, we can sample $X$ and average $\mathbb{E}[Y \mid X]$ instead of $Y$ itself. The estimator remains unbiased: $\mathbb{E}[Y] = \mathbb{E}[\mathbb{E}[Y \mid X]]$, so averaging the conditional mean gives the same answer as averaging $Y$. Its per-path variance drops from $\text{Var}(Y)$ to just $\text{Var}(\mathbb{E}[Y \mid X])$, the second term. The reduction equals the first term: exactly the within-$X$ spread we identified in the picture. This is conditional Monte Carlo.

For our barrier option a natural choice is $X = S_T$. The call payoff depends on $S_T$ directly and the survival indicator is correlated with $S_T$, so $S_T$ absorbs most of the variation in $Y$. The conditional expectation $\mathbb{E}[Y \mid S_T]$ is the call payoff times the probability that the path stayed above the barrier given the endpoints. For continuously-monitored barriers on Constant-coefficient GBM the Brownian bridge gives that survival probability in closed form. Applying the [general-case bridge survival formula]({{< relref "brownian_bridge.md" >}}#general-case) to the log-price $\ln S_t$:

$$P(\min_{0 \leq t \leq T} S_t > B \mid S_0, S_T) = 1 - \exp\left(-\frac{2 \ln(S_0/B)\,\ln(S_T/B)}{\sigma^2 T}\right)$$

The bridge integrates over all possible continuous paths connecting $S_0$ and $S_T$, returning the survival probability without us ever needing to simulate intermediate points.

Computationally this changes the simulation in a more fundamental way than the other techniques. The other four techniques all simulate full paths of $M$ steps and modify how the resulting payoffs are combined into the estimator. Conditional MC via the one-step bridge skips the intermediate simulation entirely: we sample only $S_T$, apply the bridge formula once, and combine with the call payoff. The simulation dimension collapses from $M = 252$ to $M = 1$. This is a side effect of the conditional expectation being available in closed form, not the source of the variance reduction itself.

The estimator is:

$$\hat{V}_{CMC} = \frac{e^{-rT}}{N}\sum_{i=1}^N \max(S_T^{(i)} - K, 0) \cdot P(\text{survival} \mid S_0, S_T^{(i)})$$

where $S_T^{(i)}$ is drawn from the lognormal terminal distribution and the survival probability is computed in closed form.

**Implementation.**

```python
def conditional_mc(S0, K, B, sigma, r, T, N):
    Z = np.random.standard_normal(N)
    ST = S0 * np.exp((r - 0.5*sigma**2)*T + sigma*np.sqrt(T)*Z)
    above = ST > B
    log_S0_B = np.log(S0 / B)
    log_ST_B = np.log(np.maximum(ST / B, 1e-12))
    survival_prob = np.where(above, 1 - np.exp(-2 * log_S0_B * log_ST_B / (sigma**2 * T)), 0.0)
    payoffs = np.maximum(ST - K, 0) * survival_prob
    price = np.exp(-r * T) * np.mean(payoffs)
    se = np.exp(-r * T) * np.std(payoffs, ddof=1) / np.sqrt(N)
    return price, se
```

Notice the implementation does not take $M$ as an argument. There are no intermediate steps to simulate.

**Result.** The per-path variance reduction is modest, around 1.2x. This ties back to control variates: the vanilla call (a function of $S_T$ alone) already explained roughly 9/10 of the barrier payoff's variance via linear projection. The full conditional expectation $\mathbb{E}[Y \mid S_T]$ that conditional MC averages is itself nearly linear in the vanilla payoff over the in-the-money region, so it captures almost the same variance the linear projection did. Most of $\text{Var}(Y)$ was already attributable to $S_T$, leaving only a small residual for conditional MC to remove. But because each path now costs roughly $1/M$ of a naive MC path, the wall-clock efficiency improvement is closer to **250x**. The right comparison metric for this technique is standard error per unit of computational time, not per path.

Conditional MC also removes the discretisation bias of the naive estimator. The naive simulation monitors the barrier only on the $M$ grid points and misses between-grid breaches, biasing the price upward. The naive, antithetic, control variate, importance sampling, and quasi-MC estimators all converge near **12.19**, while conditional MC (which handles continuous monitoring analytically through the bridge) converges to the true value of **11.87**.

**When to reach for it.** Conditional MC via the one-step bridge is most powerful when the conditional expectation has a clean closed form for the full interval. Continuously-monitored barrier options under constant-coefficient GBM are the canonical case. For more complex models (local vol, stochastic vol, time-dependent barriers) the one-step shortcut is not available, and one falls back to a multi-step version where the bridge is evaluated segment by segment. The multi-step version still removes bias and provides modest variance reduction, but the dramatic computational speedup is lost. Outside barriers, the same principle applies to any payoff where the path can be integrated out analytically conditional on a low-dimensional summary.

For discretely-monitored barriers, the one-step shortcut also does not apply: there is no discretisation bias to correct (the simulation can check the barrier on the actual fixing dates), the conditional expectation given the terminal value has no closed form, and the dimensional collapse is therefore gone.

---

## Technique 4: Importance Sampling
 
In the naive simulation, more than half of paths contribute zero to the payoff: some knock out at the barrier, others survive but finish out-of-the-money. The variance of the estimator comes from the productive paths, those that both survive and end in-the-money, where the payoff varies widely from path to path. We would like to draw more paths from this productive region. But the price we want is an expectation under the risk-neutral measure $\mathbb{Q}$, so any shift in the sampling distribution must be accompanied by a correction that recovers $\mathbb{Q}$-expectations.
 
This is exactly the kind of reweighting developed in the [Girsanov article]({{< relref "girsanov.md" >}}). There we showed that drift lives in the probability weights: shifting the drift of a process is the same as reweighting paths by the Radon-Nikodym derivative. We use the same machinery here but for a different purpose. In the pricing context Girsanov takes us from the real-world measure $\mathbb{P}$ to the risk-neutral measure $\mathbb{Q}$. In Monte Carlo we already start under $\mathbb{Q}$; we now reweight from $\mathbb{Q}$ to some shifted measure $\tilde{\mathbb{Q}}$ that makes productive paths more likely, then divide back by the Radon-Nikodym derivative to recover the original $\mathbb{Q}$-expectation.
 
Under $\mathbb{Q}$, the underlying follows $dS_t = rS_t\,dt + \sigma S_t\,dW_t$. Under the shifted measure $\tilde{\mathbb{Q}}$, it follows $dS_t = (r + \eta)S_t\,dt + \sigma S_t\,d\tilde{W}_t$, where $\eta$ is the drift shift we choose. The Radon-Nikodym derivative is:
 
$$\frac{d\mathbb{Q}}{d\tilde{\mathbb{Q}}} = \exp\left(-\frac{\eta}{\sigma}\tilde{W}_T - \frac{1}{2}\frac{\eta^2}{\sigma^2}T\right)$$
 
The importance sampling estimator is:
 
$$\hat{V}_{IS} = \frac{e^{-rT}}{N}\sum_{i=1}^N g(S^{(i)}) \cdot \frac{d\mathbb{Q}}{d\tilde{\mathbb{Q}}}(S^{(i)})$$
 
where the paths $S^{(i)}$ are now simulated under $\tilde{\mathbb{Q}}$.
 
**Choosing $\eta$.** The choice involves a tradeoff between two effects. A larger shift pushes more paths into the productive region (surviving and in-the-money), which lowers variance. But a larger shift also increases the variance of the likelihood ratio itself. Writing $L$ for the likelihood ratio $d\mathbb{Q}/d\tilde{\mathbb{Q}}$, we have $\log L = -(\eta/\sigma)\tilde W_T - \tfrac{1}{2}(\eta/\sigma)^2 T$, which has variance $\eta^2 T / \sigma^2$. Since $L$ is the exponential of this, the variance of $L$ grows even faster, and the importance-sampled estimator inherits that growth. The optimal $\eta$ sits where the two effects balance; pushing past it increases estimator variance and the technique can do worse than naive MC. A good rule of thumb for an out-of-the-money option is to shift the drift so that the new expected terminal price lands near the strike. For a barrier option with the barrier below spot, we shift the drift upward to keep paths away from the barrier and to push the terminal distribution toward the in-the-money region.
 
**Implementation.**
 
```python
def importance_sampling_mc(S0, K, B, sigma, r, T, M, N, eta):
    dt = T / M
    Z = np.random.standard_normal((N, M))
    # simulate under shifted measure with drift r + eta
    increments = (r + eta - 0.5*sigma**2)*dt + sigma*np.sqrt(dt)*Z
    paths = np.exp(np.log(S0) + np.cumsum(increments, axis=1))
    # likelihood ratio: depends on terminal Brownian motion under shifted measure
    W_T_tilde = np.sum(sigma*np.sqrt(dt)*Z, axis=1) / sigma
    likelihood = np.exp(-(eta/sigma)*W_T_tilde - 0.5*(eta/sigma)**2*T)
    survived = np.min(paths, axis=1) > B
    payoffs = np.maximum(paths[:, -1] - K, 0) * survived * likelihood
    price = np.exp(-r * T) * np.mean(payoffs)
    se = np.exp(-r * T) * np.std(payoffs) / np.sqrt(N)
    return price, se
```
 
For our example, $\eta$ around 0.2 gives a variance reduction factor of about 4.3x. The optimal shift could be found by a pilot simulation, but the result is moderately robust within a sensible range.
 
**When to reach for it.** Importance sampling shines when many paths contribute zero to the payoff, so concentrating samples in the productive region tightens the estimator: deep out-of-the-money options, rare event probabilities, and barrier options similar to our example. It requires more care than the other techniques because the choice of new measure can make or break the variance reduction. 
 
---
 
## Technique 5: Quasi-Monte Carlo
 
When we say "Monte Carlo," we usually mean simulating random paths and averaging the payoffs. But sampling and averaging are two separate steps. Random sampling is one way to draw paths from the underlying distribution, but it is not the only way. Do the samples themselves need to be random?
 
Random samples have an inefficiency. $N$ independent uniform draws on $[0, 1]$ do not produce $N$ evenly-spaced points; they produce clumps where the algorithm happened to draw close together and gaps where it did not. The clumps waste effort by oversampling the same region; the gaps leave information on the table. This dispersion is what makes the error shrink only as $1/\sqrt{N}$.
 
Under the hood, every random normal in our simulation comes from a transformation of uniform random numbers; the underlying randomness enters as uniforms in $[0, 1]^d$, where $d$ is the number of random inputs per path. Quasi-MC replaces these uniform draws with a deterministic sequence of points constructed to fill the unit cube evenly, then applies the inverse normal CDF to convert each one into a normal. The resulting normals concentrate where the probability density is high, in the same way random normals do, but without the clumps and gaps that come from independent random draws.
 
The deterministic sequences used in practice are called **low-discrepancy sequences**. Sobol is the standard choice in finance: each new point in the sequence is placed where the existing points have left the largest gap, so the first $N$ points are well-distributed in the unit cube for any $N$. Because the points are deterministic, convergence is no longer governed by the spread of a random estimator but by how evenly the points cover the unit cube; for smooth integrands this decays faster than $1/\sqrt{N}$.
 
The intuition behind the faster rate is that random points waste evaluations. With $N$ random points, two might land near each other (the second contributing little new information) while a gap elsewhere is left unsampled (no information about $f$ in that region). The fluctuating spacing means each new point's contribution to the average is itself noisy, and the noise averages out only as $1/\sqrt{N}$. Deterministic low-discrepancy points avoid both problems: each new point lands where the existing points have left the largest gap, so every evaluation contributes maximally. The estimate gets sharper as $N$ grows in a way the random average cannot match.
 
Mechanically the estimator looks identical to the naive one. The only change is in how the underlying uniform numbers are generated. Everything downstream is unchanged.
 
**Implementation.**
 
```python
from scipy.stats import norm, qmc
 
def quasi_mc(S0, K, B, sigma, r, T, M, N):
    dt = T / M
    sampler = qmc.Sobol(d=M, scramble=True)
    U = sampler.random(N)  # uniforms in (0,1)^M
    Z = norm.ppf(U)
    increments = (r - 0.5*sigma**2)*dt + sigma*np.sqrt(dt)*Z
    paths = np.exp(np.log(S0) + np.cumsum(increments, axis=1))
    survived = np.min(paths, axis=1) > B
    payoffs = np.maximum(paths[:, -1] - K, 0) * survived
    price = np.exp(-r * T) * np.mean(payoffs)
    se = np.exp(-r * T) * np.std(payoffs) / np.sqrt(N)
    return price, se
```
 
A note on the standard error. Pure Sobol is fully deterministic: same code, same $N$, same answer every time. There is no run-to-run variation, so the standard error formula does not apply. To recover an error estimate, we use scrambled Sobol: a randomising transformation applied to the Sobol sequence. Each scramble produces a different sequence, but all of them keep Sobol's uniform coverage property. Running several independently-scrambled versions and looking at the spread of their estimates gives the same kind of measurement we get for naive MC, which is what lets us compare quasi-MC against the other techniques on a like-for-like basis. If you only need a point estimate, unscrambled Sobol is a reasonable alternative.
 
**Result.** On our example, randomised Sobol gives a per-path variance reduction of around 8x. The wall-clock efficiency is lower, around 5x, because Sobol sequence generation is more expensive than the pseudorandom number generation used by the other techniques.
 
The mathematical machinery behind quasi-MC's convergence rate, including discrepancy bounds and the role of effective dimension, is beyond the scope of this article.
 
---
 
## Comparative Summary
 
| Technique | Lever | Per-path VR | Wall-clock efficiency | Key requirement |
|---|---|---|---|---|
| Antithetic variates | Change how we sample | ~1.5x | ~1.8x | Monotonicity of payoff |
| Control variates | Reshape what we average | ~10x | ~10x | Correlated quantity with known mean |
| Conditional MC (one-step bridge) | Reshape what we average | ~1.2x | ~250x | Analytical bridge for the full interval |
| Importance sampling | Reshape what we average | ~4.3x at $\eta = 0.2$ | ~4x | Good choice of shifted measure |
| Quasi-MC (randomised Sobol) | Change how we sample | ~8x | ~5x | Smooth payoff structure |
 
The numbers above are measured on the specific parameter set in this article. Wall-clock figures are measured on a single machine and depend on hardware and library versions; the per-path variance reduction is the more robust quantity to compare. The wall-clock efficiency column matters most for conditional MC, where the simulation dimension collapses from $M = 252$ to $M = 1$ and each path costs a tiny fraction of a naive MC path.
 
Most of these techniques can be layered: antithetic with control variates, antithetic with quasi-MC, control variates with quasi-MC. The combined variance reduction is typically less than the product of the individual factors because the techniques exploit overlapping structure.
 
---
 
## Conclusion
 
Our running example was deliberately simple, with an analytical price for verification, but the techniques apply much more broadly. Wherever Monte Carlo is used to price a derivative, the same questions arise about variance, bias, and wall-clock efficiency.
 
Two structural observations from the article I found interesting. Control variates and conditional MC are two applications of the same idea: both split the payoff into a piece explained by another variable and a residual, and average the explained piece while dealing with the residual separately. Control variates uses a linear function of a correlated payoff; conditional MC uses the full conditional expectation given some variable. And quasi-MC, despite the name, is not really a Monte Carlo method: pure Sobol is fully deterministic.
 
 