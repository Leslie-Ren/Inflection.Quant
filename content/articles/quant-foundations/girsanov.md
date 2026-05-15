---
title: "Drift Lives in the Measure: An Intuitive Look at Girsanov's Theorem"
date: 2026-05-06
draft: false
math: true
tags: ["Girsanov", "Change of Measure", "Radon-Nikodym derivative", "Risk-Neutral Pricing", "Stochastic Processes"]
---

## Why This Matters

We want to price a derivative. Under the real world measure $\mathbb{P}$, we face
two problems. First, we do not know the true drift $\mu$ of the underlying, and
historical estimates are notoriously unreliable. Second, even if we knew $\mu$,
taking the expected payoff under $\mathbb{P}$ would still not give the market price.
Risky cash flows must be discounted more heavily than guaranteed ones because
investors are risk averse. Pricing under $\mathbb{P}$ requires both the true
probabilities of outcomes and a model for how the market prices risk. Both are
fundamentally unobservable. So what can we do?

It turns out there is a way to bypass both problems entirely. Rather than estimating
$\mu$ and modelling risk aversion separately, we can work under a different
probability measure where pricing is simple. But this raises an immediate question:
is such a measure legitimate? And how does it relate to the real world?

This article is my attempt to answer both. We will take the First Fundamental
Theorem of Asset Pricing (FTAP) as given, use it to work out what the risk-neutral measure must look like, and then show where Girsanov's theorem enters and what it adds.

---

## What We Know Before Girsanov

### Taking FTAP as Given

The First Fundamental Theorem of Asset Pricing states that in an arbitrage-free market, there exists a measure $\mathbb{Q}$ under which the price of any asset, expressed in units of the numeraire, is a martingale. Equivalently, the price of any derivative equals the expected value of its payoff discounted by the numeraire. We will take FTAP as given here without proving it.

In principle, any positive traded asset can serve as a numeraire. Different choices lead to different probability measures, but they all give consistent prices for the same assets. For now we choose the simplest option: the risk-free bond $e^{rt}$. With this choice, the derivative price is:

$$V_0 = e^{-rT}\mathbb{E}^{\mathbb{Q}}[g(S_T)]$$

This is a remarkable result. The unknown drift $\mu$ and the unobservable risk
aversion have both disappeared. Pricing reduces to computing an expectation under
$\mathbb{Q}$. But what do we actually know about $\mathbb{Q}$?

### Two Properties That Follow Immediately

**Property 1: $\mathbb{Q}$ must price the stock correctly.**

The stock is already trading at an observed price $S_0$ today. The stock is itself a tradable payoff, in the same sense as a derivative, since it pays $S_T$ at a future date $T$. Under $\mathbb{Q}$, it must be priced consistently:

$$S_0 = e^{-rT}\mathbb{E}^{\mathbb{Q}}[S_T]$$

Equivalently, the discounted stock price $e^{-rt}S_t$ must be a martingale under
$\mathbb{Q}$.

**Property 2: The drift of the stock under $\mathbb{Q}$ must be $r$.**

This follows directly from Property 1. If the stock follows $S_t = S_0 e^{\mu t +
\sigma W_t}$ under $\mathbb{P}$, then for the discounted stock to be a martingale
under $\mathbb{Q}$, the drift must be exactly $r$. Any other drift would allow an
arbitrage between the stock and the risk-free bond: if the stock drifted faster than
$r$ under $\mathbb{Q}$, we could borrow at $r$ and buy the stock for a riskless
profit.

So the drift being $r$ under $\mathbb{Q}$ is not an assumption. It is forced on us
by the risk-free rate $r$, and no-arbitrage.

### What We Still Do Not Know

These properties tell us a great deal about $\mathbb{Q}$. But two questions
remain unanswered.

First, is $\mathbb{Q}$ actually a legitimate probability measure? FTAP guarantees
its existence abstractly, and we know its drift is $r$. But a valid probability measure must assign probabilities to all events consistently and integrate to one.

Second, how exactly does $\mathbb{Q}$ relate to $\mathbb{P}$? Knowing the drift is
$r$ under $\mathbb{Q}$ gives us one piece of information, but not the full structure of how the two measures are connected. If we want to convert expectations between $\mathbb{P}$ and $\mathbb{Q}$, or understand how the real world is reshaped into the pricing world, we need this relationship explicitly.

This is precisely what Girsanov's theorem provides.

---

## The Intuition: Reweighting Paths

### Drift Lives in the Probability Weights

Our goal is to construct a new measure $\mathbb{Q}$ under which the process has
a different drift (replacing $\mu$ with $r$). How to do that? Do we
change the paths to change the drift?

That instinct turns out to be wrong.

To see why, consider two processes. Under $\mathbb{P}$, a particle moves up with
probability 0.7 and down with probability 0.3. Under $\mathbb{Q}$, it moves up and
down with probability 0.5 each. The important point is that the set of possible paths is identical under both measures. A trajectory like up, up, down exists in both worlds unchanged. Nothing
about the paths themselves has been altered, only how likely each path is.

This is where a common intuition breaks. When people first see stochastic processes,
they may read drift directly off a single realised path. A steadily rising stock
chart is labelled as "positive drift," while a flat or noisy one is seen as "low
drift." But this is an inference from one sample path, not a property of the
process.

That interpretation fails because a single path carries no information about typical
behaviour. A zero-drift process can produce strong upward trends, and a
positive-drift process can still fall over finite horizons. What changes across
models is not the path, but the distribution over paths.

So drift is not something encoded in trajectories. It is encoded in how probability
is assigned to trajectories. If drift lives in the probability weights, then changing the drift from $\mu$ to $r$ cannot mean modifying paths. It must mean reweighting them. That reweighting of path probabilities is exactly what a change of measure does.

### The Per-Step Reweighting

Consider a random walk where the underlying has drift $\mu$ and volatility $\sigma$.
At each step, the particle moves up by $\sigma\sqrt{\Delta t}$ or down by
$\sigma\sqrt{\Delta t}$. Under $\mathbb{P}$, to match a process with drift $\mu$,
the probabilities are:

$$\mathbb{P}(\text{up}) = \frac{1}{2} + \frac{\mu\sqrt{\Delta t}}{2\sigma}, \qquad
\mathbb{P}(\text{down}) = \frac{1}{2} - \frac{\mu\sqrt{\Delta t}}{2\sigma}$$

Under $\mathbb{Q}$, we want drift $r$ instead of $\mu$, so the same step size
$\sigma\sqrt{\Delta t}$ must now be weighted to give a net drift of $r$:

$$\mathbb{Q}(\text{up}) = \frac{1}{2} + \frac{r\sqrt{\Delta t}}{2\sigma}, \qquad
\mathbb{Q}(\text{down}) = \frac{1}{2} - \frac{r\sqrt{\Delta t}}{2\sigma}$$

We now have two probability measures sitting side by side. Our goal is to compute the pricing expectation $\mathbb{E}^{\mathbb{Q}}[g(S_T)]$, but $\mathbb{Q}$ is not a distribution we can sample from directly. It is a mathematical construction whose existence FTAP guarantees, while its explicit form is still something we need to construct. What we can do is simulate paths under $\mathbb{P}$, since $\mathbb{P}$ corresponds to historically observed dynamics, or any baseline measure from which we can simulate. So we need a way to evaluate a $\mathbb{Q}$-expectation using $\mathbb{P}$-paths, by reweighting each path by how much more or less likely it is under $\mathbb{Q}$ than under $\mathbb{P}$. The formal object that does this
reweighting, the ratio $d\mathbb{Q}/d\mathbb{P}$, is called the
**Radon-Nikodym derivative**. At each step it is simply the ratio of
$\mathbb{Q}$-probability to $\mathbb{P}$-probability for the outcome that occurred.

The Radon-Nikodym derivative at each step is:

$$\frac{d\mathbb{Q}}{d\mathbb{P}}\bigg|_{\text{up}} =
\frac{1/2 + r\sqrt{\Delta t}/(2\sigma)}{1/2 + \mu\sqrt{\Delta t}/(2\sigma)}
\approx 1 - \frac{(\mu-r)\sqrt{\Delta t}}{\sigma}$$

$$\frac{d\mathbb{Q}}{d\mathbb{P}}\bigg|_{\text{down}} =
\frac{1/2 - r\sqrt{\Delta t}/(2\sigma)}{1/2 - \mu\sqrt{\Delta t}/(2\sigma)}
\approx 1 + \frac{(\mu-r)\sqrt{\Delta t}}{\sigma}$$

Since $r < \mu$ in a typical equity market, upward moves get downweighted and
downward moves get upweighted, just enough to cancel the excess drift $\mu - r$.
Note that the step size $\sigma$ appears naturally in the denominator: a larger
volatility means each step is larger, so a smaller probability adjustment is needed
to shift the drift by the same amount. Over a full path of $n = T/\Delta t$ steps:

$$\frac{d\mathbb{Q}}{d\mathbb{P}} = \prod_{i=1}^{n}
\left(1 - \frac{(\mu-r)}{\sigma}\xi_i\sqrt{\Delta t}\right)$$

where $\xi_i = \pm 1$ records whether step $i$ was up or down. Here we are treating
$\mu$, $r$, and $\sigma$ as constants, so $\theta = (\mu - r)/\sigma$ is constant
across all steps. This keeps the random walk tractable.

### Taking the Continuous Limit

Writing $\theta = (\mu - r)/\sigma$ for the constant excess drift per unit of volatility,
taking the logarithm and using $\log(1-x) \approx -x - x^2/2$ for small $x$:

$$\log\frac{d\mathbb{Q}}{d\mathbb{P}} \approx \sum_{i=1}^{n}
\left(-\theta\xi_i\sqrt{\Delta t} - \frac{\theta^2\Delta t}{2}\right)
= -\theta\sum_{i=1}^{n}\xi_i\sqrt{\Delta t} - \frac{\theta^2 T}{2}$$

As $\Delta t \to 0$, the sum $\sum_{i=1}^{n}\xi_i\sqrt{\Delta t}$ converges to
Brownian motion $W^{\mathbb{P}}_T$ by the same argument as in
[Brownian Motion: From Random Walks to Option Prices]({{< relref "understanding_brownian_motion.md" >}}).
This is a $\mathbb{P}$-Brownian motion specifically because the steps $\xi_i$ were drawn according to $\mathbb{P}$-probabilities. The product of per-step Radon-Nikodym factors becomes:

$$\frac{d\mathbb{Q}}{d\mathbb{P}} = \exp \left(-\theta W^{\mathbb{P}}_T -
\frac{\theta^2 T}{2}\right)$$

This is the **Girsanov exponential**, the continuous-time Radon-Nikodym derivative
of $\mathbb{Q}$ with respect to $\mathbb{P}$. 

The process $Z_t = \exp\left(-\theta W^{\mathbb{P}}_t - \frac{\theta^2 t}{2}\right)$
is a local martingale under $\mathbb{P}$. This follows from Itô's lemma: applying it to $Z_t$ shows that $dZ_t = -\theta Z_t dW^{\mathbb{P}}_t$, which has no $dt$ term and
is therefore a local martingale. Under standard integrability conditions, this process is in fact a true martingale[^1], so its expectation remains constant: $\mathbb{E}^{\mathbb{P}}[Z_T] = Z_0 = 1$.

[^1]: Under the Novikov condition
$\mathbb{E}^{\mathbb{P}}\left[\exp\left(\frac{1}{2}\int_0^T\theta_t^2dt\right)\right] < \infty$,
$Z_t$ is a true martingale rather than merely a local martingale. A local martingale
has zero drift locally but may fail to have constant expectation globally. The
Novikov condition rules this out and ensures $\mathbb{E}^{\mathbb{P}}[Z_T] = Z_0 = 1$
holds exactly, which is what the validity argument in the next section requires.

---
 
## Is $\mathbb{Q}$ a Valid Probability Measure?
 
We now have an explicit formula for $\mathbb{Q}$, but we should check that it is actually a legitimate probability measure. Two things are required: the probabilities must be non-negative, and they must sum to one.
 
Non-negativity is immediate. The Girsanov exponential is an exponential function, so it is strictly positive for every path.
 
Summing to one requires a short argument. The total probability assigned by $\mathbb{Q}$ to all events is:
 
$$\mathbb{Q}(\Omega) = \mathbb{E}^{\mathbb{Q}}[\mathbf{1}]$$
 
By the definition of the Radon-Nikodym derivative, any $\mathbb{Q}$-expectation can be converted to a $\mathbb{P}$-expectation by reweighting the probabilities, i.e. applying the conversion factor $Z_T = d\mathbb{Q}/d\mathbb{P}$:

$$\mathbb{E}^{\mathbb{Q}}[\mathbf{1}]= \mathbb{E}^{\mathbb{P}}[\mathbf{1}\cdot Z_T ]$$
 
We showed earlier that $Z_t$ is a martingale under $\mathbb{P}$, so $\mathbb{E}^{\mathbb{P}}[Z_T] = Z_0 = 1$. Therefore $\mathbb{Q}(\Omega) = 1$ and $\mathbb{Q}$ is a proper probability measure.
 
### Why $\mathbb{Q}$ and $\mathbb{P}$ Must Agree on What Is Possible
 
There is one further requirement that is easy to overlook. Since the Girsanov exponential is strictly positive, every path that has positive probability under $\mathbb{P}$ also has positive probability under $\mathbb{Q}$, and vice versa. The two measures agree on which events are possible and which are not. This property is called **equivalence of measures**, and it is not just a technicality.
 
To see why it matters, suppose $\mathbb{Q}$ assigned zero probability to some event that $\mathbb{P}$ considered possible, say a large downward move in the stock. Then a derivative that pays off only in that scenario would be priced at zero under $\mathbb{Q}$, even though it has a genuine chance of paying out in the real world. A trader who knew this could buy the derivative for nothing and collect a positive expected payoff under $\mathbb{P}$, which is a pure arbitrage. Equivalent measures rule this out by ensuring that anything that can happen in the real world is also priced as possible under $\mathbb{Q}$.
 
---

## Girsanov's Theorem

We now have all the pieces. Let $X_t$ be a process under $\mathbb{P}$ with drift
$\mu$ and diffusion $\sigma$:

$$dX_t = \mu_t dt + \sigma_t dW^{\mathbb{P}}_t$$

In the random walk above we used a constant $\theta = (\mu - r)/\sigma$. The continuous-time theorem makes no such restriction. Define the **market price of risk** $\theta_t$, which is now allowed to vary over time, as the excess drift removed per unit of volatility at each instant. Girsanov's theorem states that there exists a
measure $\mathbb{Q}$, whose Radon-Nikodym derivative with respect to $\mathbb{P}$
is:

$$\frac{d\mathbb{Q}}{d\mathbb{P}} = \exp\left(-\int_0^T \theta_t dW^{\mathbb{P}}_t - \frac{1}{2}\int_0^T \theta_t^2 dt\right)$$

under which $X_t$ has drift $r$ instead of $\mu$:

$$dX_t = r_tdt + \sigma_t dW^{\mathbb{Q}}_t$$

The diffusion $\sigma_t$ is unchanged between the two measures. The paths are the same. Only the drift and the probability weights have changed.

| Component | Meaning |
|---|---|
| $\theta_t = (\mu_t- r_t)/\sigma_t$ | Market price of risk: excess drift removed per unit of volatility at time $t$ |
| $\int_0^T \theta_tdW^{\mathbb{P}}_t$ | Stochastic integral capturing how each infinitesimal shock is reweighted over time|
| $\frac{1}{2}\int_0^T \theta_t^2dt$ | Deterministic quadratic-variation correction that naturally emerges from exponentiating Brownian motion and ensures the exponential remains normalized|
| $dW^{\mathbb{Q}}_t = dW^{\mathbb{P}}_t + \theta_tdt$ | Same Brownian increments, recentred to remove excess drift |

This answers both open questions from earlier. The Girsanov exponential is the
explicit Radon-Nikodym derivative relating $\mathbb{P}$ and $\mathbb{Q}$. And the
theorem guarantees that this reweighting always produces a valid probability measure.

---

## The Connection to the Feynman-Kac Article

We can now close the loop with the [Feynman-Kac article]({{< relref "feynman_kac.md" >}}). There we took $\mathbb{Q}$
as given and showed that expectations under it satisfy a PDE. Here we have shown
where $\mathbb{Q}$ comes from and why it is legitimate.

The full picture across the two articles is:

1. FTAP guarantees a measure $\mathbb{Q}$ exists under which pricing is simple.
2. Girsanov tells us $\mathbb{Q}$ is obtained from $\mathbb{P}$ by the explicit
Radon-Nikodym derivative above, and that the result is always a valid measure.
3. Feynman-Kac tells us expectations under $\mathbb{Q}$ satisfy a PDE, giving us
a second way to compute prices. The PDE approach is particularly valuable when early exercise is possible, cases where
computing the expectation directly becomes intractable.

Every time we write $\mathbb{E}^{\mathbb{Q}}$ in a pricing formula, both Girsanov
and FTAP are working quietly in the background.

---

## When Does Girsanov Get Used Explicitly?

| Situation | How Girsanov is Used |
|---|---|
| Switching from $\mathbb{P}$ to $\mathbb{Q}$ |At the SDE level, it changes the drift from $\mu$ to $r$. At the distribution level, the Radon–Nikodym derivative reweights entire paths, so that the market price of risk is absorbed into probabilities rather than dynamics |
| Change of numeraire | Switching numeraire corresponds to a drift change; Girsanov guarantees each switch produces a valid measure |
| Stochastic volatility models | Volatility risk cannot be hedged; its market price of risk is a free parameter and each choice gives a valid $\mathbb{Q}$ via Girsanov |
| Importance sampling in Monte Carlo | Replace the original sampling measure with a more convenient one that makes rare payoff events more likely, and correct expectations using the Radon–Nikodym derivative to preserve unbiased pricing while reducing variance|

## Looking Ahead
So far we have taken the risk-free bond as the numeraire, which led us to the risk-neutral measure $\mathbb{Q}$ and the familiar drift of $r$. But any strictly positive self-financing wealth process can play that role, and each choice produces a different valid measure via exactly the same Girsanov machinery. The [next  article]({{< relref "change_of_numeraire.md" >}}) develops this in full, exploring what happens when we switch to other natural numeraires and how each choice simplifies the pricing of a different class of derivatives.