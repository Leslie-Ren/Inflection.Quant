---
title: "The Measure We Choose: How Numéraires Simplify Pricing"
date: 2026-05-12
draft: false
math: true
tags: ["Girsanov", "Change of Measure", "Radon-Nikodym derivative", "Risk-Neutral Pricing", "Stochastic Processes"]
---

## Why This Matters

In the [article on Girsanov's Theorem]({{< relref "girsanov.md" >}}), we studied how the real-world measure $\mathbb{P}$ and the risk-neutral measure $\mathbb{Q}$ relate, and showed that switching between them amounts to reweighting paths via the Girsanov exponential. Throughout, the risk-free bond was the numéraire: the asset against which all prices were expressed. But this is a convenient choice, not a fundamental one.

Any strictly positive self-financing wealth process can serve as a numéraire, and each choice gives a different probability measure under which asset prices, expressed in units of that numéraire, become martingales. The price of a derivative is invariant; what changes is how the problem is represented. So instead of viewing pricing as a fixed-measure expectation problem, it is often more natural to think of it as choosing the numéraire that best matches the structure of the payoff.

We will see this concretely through the exchange option, pricing it under two different measures. One choice makes the calculation almost trivial, while the other leaves the true structure of the problem hidden in plain sight.

---

## A Motivating Example: The Exchange Option

An exchange option gives the holder the right to exchange asset $S^2$ for asset $S^1$ at time $T$. Its payoff is:

$$V_T = \max(S^1_T - S^2_T, 0)$$

Assume both assets follow geometric Brownian motion under $\mathbb{P}$:

$$dS^i_t = \mu_i S^i_t dt + \sigma_i S^i_t dW^i_t, \qquad i = 1, 2$$

with $d\langle W^1, W^2 \rangle_t = \rho dt$. 

Let us price it two ways: under the risk-neutral measure, and under the stock measure using $S^2$ as numéraire.

### Approach 1: Pricing Under the Risk-Neutral Measure $\mathbb{Q}$

Under the risk-neutral measure $\mathbb{Q}$, the drifts are replaced by the risk-free rate $r$, and the terminal values are:

$$S^i_T = S^i_0 \exp\left[\left(r - \tfrac{1}{2}\sigma_i^2\right)T + \sigma_i \sqrt{T} \xi^i\right]$$

where $(\xi^1, \xi^2)$ is a bivariate standard normal with correlation $\rho$.

The price is:

$$V_0 = e^{-rT} \mathbb{E}^{\mathbb{Q}}\left[\max(S^1_T - S^2_T, 0)\right]$$

Write $X = \ln S^1_T$ and $Y = \ln S^2_T$, and define the log-ratio $Z = X - Y = \ln(S^1_T / S^2_T)$. Under $\mathbb{Q}$, $Z$ is normally distributed:

$$Z \sim \mathcal{N}\left(\ln\frac{S^1_0}{S^2_0} + \left(-\tfrac{1}{2}\sigma_1^2 + \tfrac{1}{2}\sigma_2^2\right)T,\sigma^2 T\right)$$

where $\sigma^2 = \sigma_1^2 + \sigma_2^2 - 2\rho\sigma_1\sigma_2$. The payoff condition $S^1_T > S^2_T$ is simply $Z > 0$, so we can write:

$$V_0 = e^{-rT}\mathbb{E}^{\mathbb{Q}}\!\left[e^X \mathbf{1}_{Z > 0}\right] - e^{-rT}\mathbb{E}^{\mathbb{Q}}\!\left[e^Y \mathbf{1}_{Z > 0}\right]$$

Each term involves a log-normal multiplied by an indicator on a correlated normal, so the expectations do not factorise directly. They can be evaluated using the log-normal identity[^1], applied to the pairs $(X, Z)$ and $(Y, Z)$, giving:
 
$$\boxed{V_0 = S^1_0 \, \Phi(d_1) - S^2_0 \, \Phi(d_2)}$$
 
where $d_1 = \dfrac{\ln(S^1_0/S^2_0) + \frac{1}{2}\sigma^2 T}{\sigma\sqrt{T}}$ and $d_2 = d_1 - \sigma\sqrt{T}$. 

The result is correct, but getting there required peeling apart two correlated expectations, computing covariances, and standardising shifted events. And after all that work, it is still not clear why the price depends on two stocks in that particular way. The stock measure will make it obvious.
 
[^1]: For two jointly normal random variables $A$ and $B$: $\mathbb{E}\!\left[e^A f(B)\right] = e^{\mu_A + \frac{1}{2}\sigma_A^2} \, \mathbb{E}\!\left[f(B + \text{Cov}(A, B))\right]$. Applied to the first term with $A = X$, $B = Z$: $\text{Cov}(X, Z) = (\sigma_1^2 - \rho\sigma_1\sigma_2)T$, the scalar $e^{-rT} \cdot e^{\mu_X + \frac{1}{2}\sigma_1^2 T}$ simplifies to $S^1_0$, and the shifted event standardises to $\Phi(d_1)$. The second term follows identically with $\text{Cov}(Y, Z) = (\rho\sigma_1\sigma_2 - \sigma_2^2)T$.

---

### Approach 2: Pricing Under the Stock Measure $\mathbb{Q}^2$

Now use $S^2$ as the numéraire. Under the associated measure $\mathbb{Q}^2$, any asset price *divided by $S^2$* is a martingale. The pricing formula becomes:

$$V_0 = S^2_0 \, \mathbb{E}^{\mathbb{Q}^2}\!\left[\frac{V_T}{S^2_T}\right] = S^2_0 \, \mathbb{E}^{\mathbb{Q}^2}\!\left[\max\!\left(\frac{S^1_T}{S^2_T} - 1, \, 0\right)\right]$$

Define $R_t = S^1_t / S^2_t$. The option has become a standard call on $R$ with strike 1 under $\mathbb{Q}^2$.

**What does $R_t$ look like under $\mathbb{Q}^2$?** Since $S^1/S^2$ must be a martingale under $\mathbb{Q}^2$, by Itô's formula:

$$\frac{d R_t}{R_t} = (\ldots) \, dt + \sigma_1 \, d\widetilde{W}^1_t - \sigma_2 \, d\widetilde{W}^2_t$$

where $\widetilde{W}^i$ are Brownian motions under $\mathbb{Q}^2$. The drift must be zero. The volatility is unchanged between measures, so:

$$\frac{d R_t}{R_t} = \sigma \, d\widetilde{W}_t$$

where $\sigma = \sqrt{\sigma_1^2 + \sigma_2^2 - 2\rho\sigma_1\sigma_2}$ is the relative volatility and $\widetilde{W}$ is a single Brownian motion under $\mathbb{Q}^2$.

So $R_T$ is log-normal with no drift, volatility $\sigma$, and initial value $R_0 = S^1_0/S^2_0$. The option is a call on $R$ with strike 1 — exactly a Black-Scholes call:

$$V_0 = S^2_0 \left[R_0 \, \Phi(d_1) - 1 \cdot \Phi(d_2)\right] = S^1_0\,\Phi(d_1) - S^2_0\,\Phi(d_2)$$

Same answer, and this derivation is more direct. Choosing $S^2$ as numéraire turns the exchange option into a call on the ratio $S^1/S^2$, with no covariance bookkeeping and no shifted events. More importantly, the stock measure reveals what the risk-neutral measure keeps hidden: the exchange option is fundamentally a bet on the relative performance of two assets, and $\sigma = \sqrt{\sigma_1^2 + \sigma_2^2 - 2\rho\sigma_1\sigma_2}$ is the only risk that matters.

---

## The Stock Measure and Its Relation to $\mathbb{Q}$
 
What exactly happens to the measure when we switch numéraire from $B_t$ to $S^2_t$? The exchange option showed us that $R_t = S^1_t/S^2_t$ becomes driftless, but not why. Here we make that clear.
 
$$\mathbb{Q} \xrightarrow{\text{change of numéraire}} \mathbb{Q}^2$$
 
The numéraire for $\mathbb{Q}$ is the bond $B_t = e^{rt}$. The numéraire for $\mathbb{Q}^2$ is $S^2_t$. To find the Radon-Nikodym derivative between them, recall that any asset $V_t$, expressed in units of $B_t$, must be a $\mathbb{Q}$-martingale:
 
$$\frac{V_t}{B_t} = \mathbb{E}^{\mathbb{Q}}\!\left[\frac{V_T}{B_T} \,\bigg|\, \mathcal{F}_t\right]$$
 
We want the same asset, now expressed in units of $S^2_t$, to be a $\mathbb{Q}^2$-martingale. Using the change of measure formula $\mathbb{E}^{\mathbb{Q}^2}[X] = \mathbb{E}^{\mathbb{Q}}\!\left[X \cdot \frac{d\mathbb{Q}^2}{d\mathbb{Q}}\right]$, this requires:
 
$$\frac{V_t}{S^2_t} = \mathbb{E}^{\mathbb{Q}^2}\!\left[\frac{V_T}{S^2_T} \,\bigg|\, \mathcal{F}_t\right] = \mathbb{E}^{\mathbb{Q}}\!\left[\frac{V_T}{S^2_T} \cdot \frac{d\mathbb{Q}^2}{d\mathbb{Q}} \,\bigg|\, \mathcal{F}_t\right]$$
 
Comparing the two expressions, the unique process that makes this hold for every asset $V$ is:
 
$$\left.\frac{d\mathbb{Q}^2}{d\mathbb{Q}}\right|_{\mathcal{F}_T} = \frac{S^2_T / S^2_0}{B_T / B_0} = \frac{S^2_T}{S^2_0 \, e^{rT}}$$
 
The ratio of the two numéraires, normalised at time 0, is the only Radon-Nikodym derivative consistent with both martingale conditions simultaneously.
 
This is the discounted price of $S^2$, which is a $\mathbb{Q}$-martingale as required. Explicitly, since $S^2_T = S^2_0 \exp\!\left[(r - \frac{1}{2}\sigma_2^2)T + \sigma_2 \sqrt{T}\, \xi^2\right]$ under $\mathbb{Q}$:
 
$$\frac{d\mathbb{Q}^2}{d\mathbb{Q}} = \exp\!\left(-\frac{1}{2}\sigma_2^2 T + \sigma_2 \sqrt{T}\, \xi^2\right)$$
 
This is the Girsanov exponential with $\theta = -\sigma_2$, which shifts the drift of $W^2$ by $\sigma_2$. Under $\mathbb{Q}^2$, the Brownian motion becomes $\widetilde{W}^2_t = W^2_t - \sigma_2 t$. Substituting into the SDEs, the drift of each asset picks up an additional term from the volatility of the numéraire:
 
$$dS^1_t = (r + \rho\sigma_1\sigma_2) S^1_t \, dt + \sigma_1 S^1_t \, d\widetilde{W}^1_t$$
$$dS^2_t = (r + \sigma_2^2) S^2_t \, dt + \sigma_2 S^2_t \, d\widetilde{W}^2_t$$
 
For $S^1$, the extra drift comes through the correlation with $W^2$. For $S^2$, it is exactly $\sigma_2^2$, the quadratic variation of the numéraire itself. Applying Itô's lemma to the ratio $R_t = S^1_t/S^2_t$, these extra drift terms cancel exactly and we recover the driftless SDE from Section 2:
 
$$dR_t = \sigma R_t \, d\widetilde{W}_t$$
 
confirming that $R_t = S^1_t/S^2_t$ is a martingale under $\mathbb{Q}^2$.

The key point is that the measure change is driven entirely by the randomness of the new numéraire. The Radon-Nikodym derivative depends only on the Brownian shocks of $S^2$. Changing numéraire therefore reweights paths according to how $S^2$ fluctuates. Assets that are correlated with those fluctuations inherit a drift adjustment under the new measure. The shift from $\mathbb{Q}$ to $\mathbb{Q}^2$ is therefore sourced by the volatility and correlation structure of the numéraire, not by investor risk preferences.
 
| Measure | $dS^1_t$ drift | $dS^2_t$ drift |
|---|---|---|
| $\mathbb{Q}$ | $r$ | $r$ |
| $\mathbb{Q}^2$ | $r + \rho\sigma_1\sigma_2$ | $r + \sigma_2^2$ |
 
---

## The General Change of Numéraire Formula

The exchange option gave us one instance of a numéraire change. The pattern it revealed generalises. Every time we price a derivative, we are implicitly choosing a numéraire. Choosing it to match the structure of the payoff is what makes the difference between a one-line derivation and a page of algebra.
 
**Definition.** A *numéraire* is any strictly positive adapted process $N_t$. The associated measure $\mathbb{Q}^N$ is the probability measure under which every asset price $V_t$, expressed in units of $N_t$, is a martingale:
 
$$\frac{V_t}{N_t} = \mathbb{E}^{\mathbb{Q}^N}\!\left[\frac{V_T}{N_T} \,\bigg|\, \mathcal{F}_t\right]$$
 
**Theorem (Change of Numéraire).** Let $M$ and $N$ be two numéraires. Then the measures $\mathbb{Q}^M$ and $\mathbb{Q}^N$ are related by:
 
$$\frac{d\mathbb{Q}^N}{d\mathbb{Q}^M}\bigg|_{\mathcal{F}_T} = \frac{N_T / N_0}{M_T / M_0}$$
 
**Pricing formula.** The price of a claim with payoff $V_T$ at time $T$ is the same under any numéraire:
 
$$V_0 = N_0 \, \mathbb{E}^{\mathbb{Q}^N}\!\left[\frac{V_T}{N_T}\right] = M_0 \, \mathbb{E}^{\mathbb{Q}^M}\!\left[\frac{V_T}{M_T}\right]$$
 
The practical question is which numéraire makes $V_T / N_T$ easiest to work with. In the exchange option, choosing $N_t = S^2_t$ turned $V_T/N_T = \max(R_T - 1, 0)$ into a call on a driftless GBM. This pattern shows up across many problems: find the numéraire that makes the key ratio driftless, and the pricing reduces to a standard calculation. In practice, the right numéraire is usually the one already hiding inside the payoff.

---

## Applications
 
### The Forward Measure
 
In the equity world, the risk-free rate is constant and discounting is trivial. In interest rate models, the discount factor $e^{-\int_0^T r_s ds}$ is stochastic and correlated with the payoff. This is what makes rates derivatives hard under $\mathbb{Q}$: the expectation $\mathbb{E}^{\mathbb{Q}}[e^{-\int_0^T r_s ds} V_T]$ does not factorise.
 
The fix is to use the zero-coupon bond $P(t, T)$ as numéraire. Under the associated $T$-forward measure $\mathbb{Q}^T$:
 
$$V_0 = P(0, T) \, \mathbb{E}^{\mathbb{Q}^T}[V_T]$$
 
The stochastic discounting disappears entirely, absorbed into the measure change. The forward price $F(t,T) = \mathbb{E}^{\mathbb{Q}^T}[X_T \,|\, \mathcal{F}_t]$ is a martingale under $\mathbb{Q}^T$ by construction, which is exactly the driftless ratio we saw in the exchange option, just with $P(t,T)$ playing the role of $S^2$.
 
Take caplets as an example, the payoff on the compounded SOFR rate over $[T,T+\tau]$, paid at $T+\tau$, becomes a call on the forward rate under the $T+\tau$-forward measure. Assuming log-normality, Black’s formula applies immediately.

### The Annuity Measure
 
An interest rate swap has a more complex structure: its value depends on multiple payment dates $T_1, \ldots, T_n$. The natural numéraire is not a single bond but the annuity:
 
$$A_t = \sum_{i=1}^{n} \tau_i P(t, T_i)$$
 
the present value of receiving one unit on each payment date. Under the swap measure $\mathbb{Q}^A$, the swap rate $S_t$ (the fair fixed rate on the swap) is a martingale. This is the analogue of $R_t = S^1_t/S^2_t$ being driftless under $\mathbb{Q}^2$: the swap rate is driftless under its natural measure.
 
A payer swaption pays $A_T\max(S_T - K, 0)$ at expiry. Under $\mathbb{Q}^A$:
 
$$V_0 = A_0 \, \mathbb{E}^{\mathbb{Q}^A}[\max(S_T - K, 0)]$$

The stochastic annuity disappears into the measure change, leaving a plain option on the swap rate. Assuming $S_T$ is log-normal under $\mathbb{Q}^A$, this becomes Black's formula for interest rate swaptions.
 
### The FX Measure

In the equity world, there is one currency and one risk-free rate. In FX, there are two, and any valuation across currencies necessarily introduces an exchange rate. Working directly under $\mathbb{Q}^d$ therefore becomes cumbersome.

Let $X_t$ be the spot rate (domestic per foreign). Under the domestic risk-neutral measure $\mathbb{Q}^d$:
$$
dX_t = (r_d - r_f) X_t \, dt + \sigma X_t \, dW_t
$$

The drift reflects the funding differential between the two economies: holding foreign currency earns $r_f$ instead of $r_d$, so the exchange rate must compensate through its drift under domestic pricing.

A key question is what object should play the role of numéraire in FX. Although $X_t$ enters directly into pricing problems, it is not itself a traded wealth process. Holding foreign currency is not static: it accrues interest at the foreign short rate $r_f$. Any self-financing position in foreign currency therefore grows like a money market account. This means the natural tradable benchmark is not $X_t$ alone, but the foreign money market account expressed in domestic currency:
$$
N_t = X_t e^{r_f t}
$$

This distinction is important. Using $X_t$ alone ignores the carry from holding foreign cash and therefore does not correspond to any replicable trading strategy.

The practical benefit of adopting the FX numéraire is that pricing becomes structurally simpler when payoffs naturally align with this choice of numéraire. Instead of carrying the interest rate differential through every expectation under $\mathbb{Q}^d$, it is absorbed into the measure itself. What remains is a cleaner decomposition between FX risk and the underlying payoff structure. The real power of this perspective becomes more visible in multi-currency derivatives, where FX interacts directly with foreign underlyings. These products, in particular quanto and composite structures, will be the focus of the next article.

---
## Conclusion
Pricing is not tied to a single measure but is a choice of representation, and the optimal numéraire is the one that matches the structure of the payoff. The exchange option, the forward measure, and the annuity measure all illustrate the same principle from different angles: the measure change does not alter the answer, it reveals structure the risk-neutral measure keeps hidden. In the next article we apply this framework to currency converted derivatives, where the choice of numéraire becomes not just convenient but essential.



