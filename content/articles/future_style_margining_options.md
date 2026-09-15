---
title: "Futures-Style Margined Options: The Absence of Early Exercise Premium"
date: 2026-04-23
draft: false
math: true
tags: ["Options", "Futures", "Commodities", "Margining", "American Options"]
---

## Why This Matters

When I first studied options, most textbook examples were equity-style:
you pay a premium upfront, and you receive the payoff at expiry, or whenever
you choose to exercise, if the option is American. That framing was so ingrained that I took
it for the general case.

When I started working on commodity derivatives, I encountered a different
world. Many options are traded under futures-style margining. No premium
changes hands at inception; instead, the option is margined daily like a
futures contract. The convention tends to split along venue lines rather than by
underlying. US exchanges are predominantly equity-style: options on WTI crude
futures at the CME and options on corn and wheat futures at the CBOT all
require an upfront premium. European venues lean the other way. Options on ICE
Brent futures and options on EUA carbon futures at ICE Endex and EEX are both
margined futures-style.

A standard assumption in practice is that American futures-style options are valued
identically to their European counterparts. When I first went looking for an explanation
on why early exercise has no benefit, the most common answer I found was something like:

> *Daily marking-to-market removes the time-value-of-money advantage that usually
> justifies early exercise for American options.*

That statement makes some intuitive sense, but it never gave me the mathematical comfort
I was looking for. To really understand why American and European options coincide under futures-style
margining, I found it helpful to break the problem into smaller steps along two
separate dimensions:

- **Margining convention**: futures-style margining (FSM) vs. equity-style margining
  (ESM), both applied to European options on futures.
- **Exercise style**: American vs. European.

The first step is to get a clear understanding of the margining dimension: what is the
difference between a futures-style and an equity-style European option on a futures
contract, and how does the change in margining convention affect the PDE and its valuation?
This is less obvious than it first appears. The second step, showing that the American early
exercise feature has no value under futures-style margining, turns out to require no
additional machinery. It follows from the cash flow mechanics of exercise together with one
property of $V$ that the first step has already established.


## European Options: Equity-Style vs. Futures-Style Margining

I want to compare the valuation difference from the PDE perspective. We will derive the
PDE from scratch. Let $F_t$ denote the futures price at time $t$, assumed to follow
geometric Brownian motion under the risk-neutral measure $\mathbb{Q}$:

$$dF = \sigma FdW^{\mathbb{Q}}$$

There is no drift term. Under the risk-neutral measure, futures prices are martingales since entering a futures contract requires no capital (other than the initial margin required by the exchange). Throughout, the risk-free rate $r$ is taken to be deterministic, and margin balances are assumed to accrue at that same rate.

Consider a European option with value $V = V(F, t)$. We hedge it with a short position in
$\Delta$ futures contracts and follow the P&L of the combined position, which we write as
$d\Pi$. The futures leg contributes its daily settlement $-\Delta dF$ and requires no cash to
put on.

### The Equity-Style Case

In the equity-style world, $V$ is the cash premium paid upfront. Since the futures leg
requires no cash, the strategy's entire initial outlay is the option premium $V$.

Applying Itô's lemma to $V(F, t)$:

$$dV = \frac{\partial V}{\partial t}dt + \frac{\partial V}{\partial F}dF +
\frac{1}{2}\frac{\partial^2 V}{\partial F^2}(dF)^2$$

The P&L of the strategy is:

$$d\Pi = dV - \Delta dF = \left(\frac{\partial V}{\partial t} +
\frac{1}{2}\sigma^2 F^2 \frac{\partial^2 V}{\partial F^2}\right)dt +
\left(\frac{\partial V}{\partial F} - \Delta\right)\sigma F dW^{\mathbb{Q}}$$

Setting $\Delta = \frac{\partial V}{\partial F}$ eliminates the stochastic term. The
strategy is now instantaneously risk-free:

$$d\Pi = \left(\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 F^2
\frac{\partial^2 V}{\partial F^2}\right)dt$$

**No-arbitrage condition**: a riskless position must earn the risk-free rate $r$ on the
capital it ties up. That capital is the premium $V$, so we require $d\Pi = rV dt$. Setting the two
expressions equal and rearranging:

$$\boxed{\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 F^2
\frac{\partial^2 V}{\partial F^2} - rV = 0}$$

with terminal condition $V(F, T) = $ option payoff.

The $-rV$ term is the cost of carry on the cash investment $V$. It is present because
$V$ is the money the holder has paid out and it must earn the risk-free rate to break
even. The solution is the Black formula:

$$V^{\text{equity}}(F, t) = e^{-r(T-t)}\left[F N(d_1) - K N(d_2)\right]$$

where $d_1, d_2$ are the standard Black expressions.

### The Futures-Style Case

Under futures-style margining, no cash premium is paid at inception. Instead, the option
is margined daily: if the exchange's settlement price moves from $V_t$ to $V_{t+dt}$,
the holder receives (or pays) $dV = V_{t+dt} - V_t$ through their margin account. The
option position itself requires zero initial cash outlay. How does this change the PDE and
its valuation?

#### Derivation of the Futures-Style PDE

Follow the same delta-hedged strategy. Setting $\Delta = \frac{\partial V}{\partial F}$
eliminates the stochastic term as before, and the instantaneous risk-free P&L is:

$$d\Pi = \left(\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 F^2
\frac{\partial^2 V}{\partial F^2}\right)dt$$

Now apply the no-arbitrage condition. Neither leg requires cash at inception. The option is
margined futures-style and carries no premium, and the futures hedge requires no upfront
payment either. The strategy therefore ties up no capital at all, and since it is
instantaneously riskless, any non-zero deterministic drift would imply arbitrage. Its drift
must vanish:

$$d\Pi = 0$$

$$\boxed{\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 F^2
\frac{\partial^2 V}{\partial F^2} = 0}$$

with the same terminal condition. The $-rV$ term is gone because there is no capital tied up
to carry.

#### $V$ as a Martingale Under the Risk-Neutral Measure

The absence of the $-rV$ term has a direct probabilistic interpretation. Applying
Itô's lemma to $V(F, t)$ under the risk-neutral measure and substituting
$dF = \sigma F dW^{\mathbb{Q}}$:

$$dV = \left(\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 F^2
\frac{\partial^2 V}{\partial F^2}\right)dt + \frac{\partial V}{\partial F}\sigma F dW^{\mathbb{Q}}$$

The $dt$ term is exactly the left-hand side of the futures-style PDE, so it vanishes,
leaving:

$$dV = \frac{\partial V}{\partial F}\sigma F dW^{\mathbb{Q}}$$


Therefore $V$ is a local martingale under $\mathbb{Q}$, and under the assumed dynamics a true
martingale. This is the
direct counterpart to the futures price $F$ itself being a martingale under $\mathbb{Q}$:
just as $F$ requires no discounting because entering a futures contract requires no cash
outlay, $V$ requires no discounting because the futures-style option requires no upfront
premium. Being a martingale, $V$ satisfies:

$$V(F_t, t) = \mathbb{E}^{\mathbb{Q}}\left[V(F_T, T) \middle|\mathcal{F}_t\right]
= \mathbb{E}^{\mathbb{Q}}\left[\text{Payoff}(F_T) \middle|\mathcal{F}_t\right]$$

where

$$\text{Payoff}(F_T) = \begin{cases}
\max(F_T - K, 0) & \text{call} \\
\max(K - F_T, 0) & \text{put}
\end{cases}$$

That is, the futures-style MTM at any point in time is the risk-neutral expectation of
the terminal payoff with no discount factor applied.

#### Comparing the Two Worlds

We can now contrast the two margining conventions clearly.

| | Equity-Style | Futures-Style |
|---|---|---|
| Premium at inception | Paid upfront in cash | Zero, no cash changes hands |
| PDE | $V_t + \frac{1}{2}\sigma^2F^2V_{FF} - rV = 0$ | $V_t + \frac{1}{2}\sigma^2F^2V_{FF} = 0$ |
| What $V$ represents | Present value of the option | Exchange MTM settlement price |
| Probabilistic form | $e^{-r(T-t)}\mathbb{E}^{\mathbb{Q}}[\text{payoff}]$ | $\mathbb{E}^{\mathbb{Q}}[\text{payoff}]$ |
| Cash flow to holder | Premium $V$ paid at $t_0$, payoff received at $T$ | Daily margin flows $dV$, summing to payoff at $T$ |

In the **equity-style** world, $V(F, t)$ is the fair cash amount to exchange today for
the right to receive the option payoff at expiry. It is a present value in the
traditional sense. 

In the **futures-style** world, $V(F, t)$ is the exchange's mark-to-market settlement quote, used to compute each day's margin flow. It is not paid or received as a lump sum.

Since $V$ is the undiscounted expectation of the payoff, the closed-form solution is Black's
formula with the discount factor removed:

$$V^{\text{futures}}(F, t) = F N(d_1) - K N(d_2)$$

Comparing with the equity-style solution:

$$V^{\text{futures}} = e^{r(T-t)} V^{\text{equity}}$$

The futures-style MTM exceeds the equity-style present value by exactly $e^{r(T-t)}$. The futures-style holder collects the same economic cash flows as the equity-style holder but without paying anything upfront, so the quoted price is scaled up by the
cost of carry that the equity-style holder effectively prepays.

## American Options Under Futures-Style Margining

We now turn to the central question. In the equity-style world, American options can be worth more than European options. Early exercise can be optimal when the intrinsic value in hand, reinvested at $r$, exceeds the value of waiting (discussed in [Early Exercise of American Options: Call Equivalence and the Put Premium]({{< relref "american_vs_european_options.md" >}})). Does the same logic apply under futures-style margining?

### Setup and Notation

Consider an American put option on a futures contract, traded under futures-style
margining, with strike $K$, expiry at time $T$, and current time $t_0$. The exchange
publishes a daily MTM settlement price for the option, which we denote $V_i = V(F_i, t_i)$ on
day $i$. Recall that $V_i$ is not a present value but the exchange-quoted settlement price
used to compute each day's margin flow. The holder receives $V_i - V_{i-1}$ on day $i$
through their margin account.

At expiry on day $n$, the settlement price converges to intrinsic value:

$$V_n = \max(K - F_n, 0)$$

### Exercise Mechanics

When the holder of a futures-style American put exercises on day $m$, the following
happens in sequence:

1. The regular daily margin flow $V_m - V_{m-1}$ is settled as usual through the margin account. This happens regardless of exercise.
2. The option position is submitted for exercise. The exchange assigns the holder a short futures position at the strike price $K$. Since the current futures price is $F_m$, this newly assigned position is immediately marked to market, and the margin account is credited with $K - F_m$ (assuming the put is in the money). The holder may then close out the short futures position at $F_m$ at no further cost, or carry it forward.
3. The option is extinguished. No further option margin flows occur from day $m+1$ onward.

### Early Exercise on Day $m$

Suppose the holder reaches day $m$, before expiry on day $n$, and considers exercising.

| As of day $m$ | Hold | Exercise |
|---|---|---|
| Day $m$ margin flow | $V_m - V_{m-1}$ | $V_m - V_{m-1}$ |
| Exercise payoff | none | $\max(K - F_m, 0)$ |
| Days $m+1$ to $n$ | flows telescoping to $V_n - V_m$ | none, the option is terminated |
| Value as of day $m$ | $V_m$ | $\max(K - F_m, 0)$ |

The question is how the intrinsic value received at exercise compares to the option value
$V_m$ given up. The European section already gives us the relationship we need:

$$V(F_t, t) = \mathbb{E}^{\mathbb{Q}}\left[\text{Payoff}(F_T) \middle|\mathcal{F}_t\right]$$

What stands in the way is the expectation itself, and we can use Jensen's inequality to remove
it since the payoff is a convex function of $F_T$.

$$\mathbb{E}^{\mathbb{Q}}\left[\text{Payoff}(F_T) \middle| \mathcal{F}_t\right] \geq
\text{Payoff}\left(\mathbb{E}^{\mathbb{Q}}\left[F_T \middle| \mathcal{F}_t\right]\right) =
\text{Payoff}(F_t)$$

Two separate facts are at work in that line. Jensen moves the expectation inside the payoff
function, and the martingale property of $F$ then replaces the expected terminal futures price
with $F_t$ itself. As a result, the comparison lands on intrinsic value at today's price,
which is the quantity exercise actually delivers.

What Jensen gives is a weak inequality, that holding is worth at least as much as exercising.
For $V_m$ and intrinsic value to be worth the same, one of two conditions has to hold.

- $F_T$ is no longer random. With $\sigma = 0$ the futures price is already known, and there is
  no time value left to give up.
- The payoff is linear over every value $F_T$ can reach. Suppose $F_T$ were certain to land
  below $K$, so the put finishes in the money in every state. The payoff is then $K - F_T$
  throughout, and averaging it gives $K$ minus the expected futures price, which is $K - F_m$.
  That is intrinsic value today. How the probability is spread out below the strike makes no
  difference, because averaging a straight line and evaluating it at the average are the same
  operation.

However, neither condition holds in practice. Volatility is not zero, and under lognormal dynamics
$F_T$ can land anywhere in $(0, \infty)$, so the range of outcomes always straddles $K$ and the
payoff always bends somewhere inside it. The inequality is strict, and on day $m$:

$$V_m > \max(K - F_m, 0) \quad  \text{for all } m \lt n$$

So the American feature has no value. We worked with a put throughout, but the argument rests
only on convexity of the payoff, so the call case needs no separate treatment.

$$\boxed{V^{\text{American, futures-style}} = V^{\text{European, futures-style}}}$$


## Takeaway

An American futures-style option can be valued with the Black model with discounting removed.
The usual machinery for equity-style American options, a numerical PDE solve or a closed-form
American approximation, is unnecessary here. That matters most for calibration. Backing out an
implied volatility still needs a root finder, but against a numerical pricer each iteration
becomes a closed-form evaluation rather than a numerical solve, which makes building a vol
surface from these quotes considerably cheaper.