---
title: "Futures-Style Margined Options: The Absence of Early Exercise Premium"
date: 2026-04-23
draft: false
math: true
tags: ["Options", "Futures", "Derivatives", "Margining"]
---

## Why This Matters

When I first studied options, most textbook examples were equity-style:
you pay a premium upfront, and at expiry (or whenever you choose to exercise, if the
option is American), you receive the payoff. That framing stuck with me for a long time.

When I started working on commodity derivatives, I encountered a different world. Many
options in commodity markets are traded under futures-style margining. No premium
changes hands at inception, and instead the option is margined daily like a futures
contract. This is common across a wide range of exchange-traded products: options on WTI
crude oil futures at the CME, options on Henry Hub natural gas futures, options on corn
and wheat futures, and options on carbon emissions futures, to name a few.

A standard assumption in practice is that American futures-style options are valued
identically to their European counterparts. When I first went looking for an explanation
on why early exercise has no benefit, the most common answer I found was something like:

> *Daily marking-to-market removes the time-value-of-money advantage that usually
> justifies early exercise for American options.*

That statement makes some intuitive sense, but it never gave me the mathematical comfort
I was looking for.

To really understand why American and European options coincide under futures-style
margining, I found it helpful to break the problem into two smaller steps along two
separate dimensions:

- **Margining convention**: futures-style margining (FSM) vs. equity-style margining
  (ESM), both applied to options on futures.
- **Exercise style**: American vs. European.

The first step is to get a clear understanding of the margining dimension: what is the
difference between a futures-style and an equity-style European option on a futures
contract, and how does the change in margining convention affect the PDE and the meaning
of the quantities in it? This is less obvious than it first appears, and getting it right
is the key to everything that follows.

The second step, showing that the American early exercise feature has no value under
futures-style margining, turns out to require no additional mathematical heavy lifting. 
It follows almost immediately from a simple cash flow argument.

---

## European Options — Equity-Style vs. Futures-Style Margining

I want to compare the valuation difference from the PDE perspective. We will derive the
PDE from scratch. Let $F_t$ denote the futures price at time $t$, assumed to follow
geometric Brownian motion under the risk-neutral measure $\mathbb{Q}$:

$$dF = \sigma FdW^{\mathbb{Q}}$$

There is no drift term. Under the risk-neutral measure, futures prices are martingales since entering a futures contract requires no capital (other than the initial margin required by the exchange).

Consider a European option with value $V = V(F, t)$. We construct a delta-hedged
portfolio $\Pi$ consisting of a long position in the option and a position in
$\Delta$ futures contracts:

$$\Pi = V - \Delta F$$

where $\Delta F$ denotes the notional of the futures hedge, not its market value (futures
have zero value at inception).

### The Equity-Style Case

In the equity-style world, $V$ is the cash premium paid upfront. Since futures require
no upfront investment, the cost of the portfolio is simply the option value $V$.

Applying Itô's lemma to $V(F, t)$:

$$dV = \frac{\partial V}{\partial t}dt + \frac{\partial V}{\partial F}dF +
\frac{1}{2}\frac{\partial^2 V}{\partial F^2}(dF)^2$$

The change in the portfolio value is:

$$d\Pi = dV - \Delta dF = \left(\frac{\partial V}{\partial t} +
\frac{1}{2}\sigma^2 F^2 \frac{\partial^2 V}{\partial F^2}\right)dt +
\left(\frac{\partial V}{\partial F} - \Delta\right)\sigma F dW^{\mathbb{Q}}$$

Setting $\Delta = \frac{\partial V}{\partial F}$ eliminates the stochastic term. The
portfolio is now instantaneously risk-free:

$$d\Pi = \left(\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 F^2
\frac{\partial^2 V}{\partial F^2}\right)dt$$

**No-arbitrage condition**: a risk-free portfolio must earn the risk-free rate $r$.
Since the cash invested equals $V$, we require $d\Pi = rV dt$. Setting the two
expressions equal and rearranging:

$$\boxed{\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 F^2
\frac{\partial^2 V}{\partial F^2} - rV = 0}$$

with terminal condition $V(F, T) = $ option payoff.

The $-rV$ term is the cost of carry on the cash investment $V$. It is present because
$V$ is the money the holder has paid out and it must earn the risk-free rate to break
even. The solution is the Black (1976) formula with a discount factor:

$$V^{\text{equity}}(F, t) = e^{-r(T-t)}\left[F N(d_1) - K N(d_2)\right]$$

where $d_1, d_2$ are the standard Black expressions.

---
### The Futures-Style Case

#### What Changes

Under futures-style margining, no cash premium is paid at inception. Instead, the option
is margined daily: if the exchange's settlement price moves from $V_t$ to $V_{t+dt}$,
the holder receives (or pays) $dV = V_{t+dt} - V_t$ through their margin account. The
option position itself requires zero initial cash outlay.

This changes the no-arbitrage argument in a subtle but critical way.

#### Derivation of the Futures-Style PDE

Construct the same delta-hedged portfolio.
Setting $\Delta = \frac{\partial V}{\partial F}$ eliminates the stochastic term as
before. The instantaneous risk-free P&L of this portfolio is:

$$d\Pi = \left(\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 F^2
\frac{\partial^2 V}{\partial F^2}\right)dt$$

Now apply the no-arbitrage condition. The portfolio consists of:

- Option position (no cash outlay — futures-style, no premium paid)
- Futures position (no cash outlay — futures require no upfront payment)

The total cash invested in this portfolio is zero. Since the hedged portfolio requires
no initial capital and is instantaneously riskless, any non-zero deterministic drift
would imply arbitrage. Therefore its drift must vanish:

$$d\Pi = 0$$

$$\boxed{\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 F^2
\frac{\partial^2 V}{\partial F^2} = 0}$$

with the same terminal condition. The $-rV$ term is gone because there is no cash
investment to carry.

#### $V$ as a Martingale Under the Risk-Neutral Measure

The absence of the $-rV$ term has a direct probabilistic interpretation. Applying
Itô's lemma to $V(F, t)$ under the risk-neutral measure and substituting
$dF = \sigma F dW^{\mathbb{Q}}$:

$$dV = \left(\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 F^2
\frac{\partial^2 V}{\partial F^2}\right)dt + \frac{\partial V}{\partial F}\sigma F dW^{\mathbb{Q}}$$

$V$ satisfies the futures-style PDE, which requires
$\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 F^2
\frac{\partial^2 V}{\partial F^2} = 0$. The $dt$ term therefore vanishes,
leaving:

$$dV = \frac{\partial V}{\partial F}\sigma F dW^{\mathbb{Q}}$$


Therefore $V$ is a martingale under $\mathbb{Q}$. This is the
direct counterpart to the futures price $F$ itself being a martingale under $\mathbb{Q}$:
just as $F$ requires no discounting because entering a futures contract requires no cash
outlay, $V$ requires no discounting because the futures-style option requires no upfront
premium. Being a martingale, $V$ satisfies:

$$V(F_t, t) = \mathbb{E}^{\mathbb{Q}}\left[V(F_T, T) \middle|\mathcal{F}_t\right]
= \mathbb{E}^{\mathbb{Q}}\left[\text{Payoff}(F_T) \middle|\mathcal{F}_t\right]$$

where $\text{Payoff}(F_T) = \max(F_T - K, 0)$ for a call and $\text{Payoff}(F_T) = \max(K - F_T, 0)$ for a put.

That is, the futures-style MTM at any point in time is the risk-neutral expectation of
the terminal payoff with no discount factor applied. We will rely on one consequence of
this in the American options section:

$$V(F_t, t) > \text{Intrinsic Value}(F_t) \quad \text{for all } t < T$$

The strict inequality holds by a simple no-arbitrage argument for the lower bound, plus an intuitive observation for the strictness. First, $V$ can never fall below intrinsic value. If it did, one could buy the option and immediately exercise it for a riskless profit. Second, $V$ must be strictly greater than intrinsic value as long as time and volatility remain, because the option holder benefits from any further favourable move in $F_T$ before expiry, while being fully protected against unfavourable moves by the
payoff floor at zero. This asymmetry, unlimited upside participation and truncated downside, always commands a strictly positive premium above intrinsic value.

#### Comparing the Two Worlds

We can now contrast the two margining conventions clearly.

| | Equity-Style | Futures-Style |
|---|---|---|
| Premium at inception | Paid upfront in cash | Zero — no cash changes hands |
| **What $V$ represents** | **Present value of the option** | **Exchange MTM settlement price** |
| PDE | $V_t + \frac{1}{2}\sigma^2F^2V_{FF} - rV = 0$ | $V_t + \frac{1}{2}\sigma^2F^2V_{FF} = 0$ |
| Probabilistic form | $e^{-r(T-t)}\mathbb{E}^{\mathbb{Q}}[\text{payoff}]$ | $\mathbb{E}^{\mathbb{Q}}[\text{payoff}]$ |
| Cash flow to holder | Premium $V$ paid at $t_0$, payoff received at $T$ | Daily margin flows $dV$, summing to payoff at $T$ |

In the **equity-style** world, $V(F, t)$ is the fair cash amount to exchange today for
the right to receive the option payoff at expiry. It is a present value in the
traditional sense. 

In the **futures-style** world, $V(F, t)$ is the exchange's mark-to-market settlement quote, used to compute each day's margin flow. It is not paid or received as a lump sum. It always strictly exceeds intrinsic value before expiry. 

The closed-form solution to the futures-style PDE is:

$$V^{\text{futures}}(F, t) = F N(d_1) - K N(d_2)$$

Comparing with the equity-style solution:

$$V^{\text{futures}} = e^{r(T-t)} V^{\text{equity}}$$

The futures-style MTM exceeds the equity-style present value by exactly $e^{r(T-t)}$. The futures-style holder collects the same economic cash flows as the equity-style holder but without paying anything upfront, so the quoted price is scaled up by the
cost of carry that the equity-style holder effectively prepays.

## American Options Under Futures-Style Margining

We now turn to the central question. In the equity-style world, American options can be worth more than European options. Early exercise can be optimal when the intrinsic value in hand, reinvested at $r$, exceeds the value of waiting (discussed in [Early Exercise of American Options: Call Equivalence and the Put Premium]({{< relref "american_vs_european_options.md" >}})).
 Does the same logic apply under futures-style margining?

The answer is no, and we can see why by simply looking at what happens to the cash flows and payoffs when the holder exercises early.

### Setup and Notation

Consider an American put option on a futures contract, traded under futures-style
margining, with strike $K$, expiry at time $T$, and current time $t_0$. The exchange
publishes a daily MTM settlement price for the option. We denote this settlement price on
day $i$ as $V_i$, where:

$$V_i = V(F_i, t_i)$$

is the futures-style option MTM as defined above. Recall that $V_i$ is not a present
value but the exchange-quoted settlement price used to compute each day's margin
flow. The holder receives $V_i - V_{i-1}$ on day $i$ through their margin account.

At expiry on day $n$, the settlement price converges to intrinsic value:

$$V_n = \max(K - F_n, 0)$$

### Exercise Mechanics

When the holder of a futures-style American put exercises on day $m$, the following
happens in sequence:

1. The regular daily margin flow $V_m - V_{m-1}$ is settled as usual through the margin account. This happens regardless of exercise.
2. The option position is submitted for exercise. The exchange assigns the holder a short futures position at the strike price $K$. Since the current futures price is $F_m$, this newly assigned position is immediately marked to market, and the margin account is credited with $K - F_m$ (assuming the put is in the money). The holder may then close out the short futures position at $F_m$ at no further cost, or carry it forward.
3. The option is extinguished. No further option margin flows occur from day $m+1$ onward.

It is the combination of steps 2 and 3 that determines whether early exercise is
beneficial. Step 2 delivers the intrinsic value $\max(K - F_m, 0)$ through futures
assignment, while step 3 forfeits all future option margin flows. The question is
whether the amount received in step 2 compensates for what is given up in step 3.

### Early Exercise on Day $m$

Suppose the holder exercises the American put early on day $m$, where $m < n$.
The complete cash flows over the life of the position are:

| Day | Cash Flow | Sign |
|-----|-----------|------|
| 1 | $V_1 - V_0$ | positive or negative |
| 2 | $V_2 - V_1$ | positive or negative |
| $\vdots$ | $\vdots$ | $\vdots$ |
| $m$ (regular margin flow) | $V_m - V_{m-1}$ | positive or negative |
| $m$ (exercise payoff) | $\max(K - F_m, 0) - V_m$ | **always $<0$** |
| $m+1, \ldots, n$ | 0 (option is extinguished) | |

On the exercise day, after the regular margin flow is settled, the holder receives the
intrinsic value $\max(K - F_m, 0)$ through futures assignment.

From the no-arbitrage argument established in the previous section:

$$V_m > \max(K - F_m, 0) \quad \text{for all } m<n$$

The exercise payoff row in the table is therefore always $< 0$.There is never a reason to exercise early, so the American feature
has no value. Although we have used a put option to illustrate the mechanics, the same argument holds symmetrically for call options.

$$\boxed{V^{\text{American, futures-style}} = V^{\text{European, futures-style}}}$$

---

## Conclusion

The central insight of this article is that the quantity $V$ means something
fundamentally different under each margining convention. In the equity-style world, $V$ is a present value, which creates a tradeoff between immediate exercise and continued optionality. In the futures-style world, $V$ is a martingale under the risk-neutral measure, a forward-like quantity that always strictly exceeds intrinsic value before expiry. Once that is understood, the conclusion for American options follows directly from the cash flow mechanics: exercising early forfeits the option's remaining time value. The early exercise feature is contractually present but economically worthless, and American and European futures-style options are priced identically.

## References

- Black, F. (1976). The pricing of commodity contracts. *Journal of Financial Economics*, 3(1), 167–179.
