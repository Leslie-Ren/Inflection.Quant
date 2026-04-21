---
title: "Early Exercise of American Options: Call Equivalence and the Put Premium"
date: 2026-04-03
draft: false
math: true
tags: ["American Option", "European Option", "Put-Call Parity", "Black–Scholes PDE"]
---

## Why This Matters

While practitioners price American puts correctly in production systems, the deeper question of *why* early exercise is sometimes optimal, and the precise conditions under which it occurs, is less often articulated rigorously. This article works through the argument, starting with why early exercise is never optimal for calls, and then showing, using the Black–Scholes PDE, when and why it becomes mandatory for puts.

For those working with options pricing, hedging, or products with embedded American optionality, a rigorous understanding of the early exercise boundary can offer useful intuition beyond what standard pricing tools provide.


---

## Put–Call Parity

For European options on a non-dividend-paying stock, put–call parity states:

$$C_0 + Ke^{-rT} = P_0 + S_0$$

Rearranging, we get:

$$C_0 = S_0 + P_0 - Ke^{-rT}$$

This shows that a European call can be thought of as:

- Owning the stock
- Holding a European put (downside protection)
- Deferring payment of the strike (earning interest on $K$ until maturity)

The value of a put option must be non-negative. So from the rearranged parity equation, we can get a lower bound for the European call option:

$$C_0 \geq S_0 - Ke^{-rT}$$

---

## Why Early Exercise of an American Call Is Suboptimal

Consider an American call. Its holder can exercise early, but is it ever optimal?

- Exercising early gives a payoff of $S_t - K$ at time $t < T$.
- If the holder does not exercise, the option is worth at least $S_t - Ke^{-r(T-t)}$, which is more valuable than the exercised payoff.

Since the stock pays no dividends, there is no economic benefit to holding the stock earlier. Exercising the call early would **forfeit the time value of the option and the interest on $K$**, making it suboptimal. Therefore, the American call is never exercised early, and its price equals that of the European call.

---

### When This Breaks Down

The result changes if the stock pays dividends. Dividends reduce the stock price on the ex-dividend date, and option holders do not receive dividends. Exercising just before a dividend can therefore be advantageous sometimes — the early exercise premium becomes 
positive, and the American call is worth more than its European counterpart.

---

## Why American Puts Are Worth More Than European Puts

For puts, the logic reverses. Waiting is no longer free — the holder defers 
**receiving** $K$, rather than deferring paying it. Every period the put is held 
unexercised, the holder foregoes interest on $K$. When that cost exceeds the 
remaining benefit of holding, early exercise is optimal.

To see this concretely, consider a put deep in-the-money with $S_t \approx 0$:

- Exercising immediately yields $K - S_t \approx K$, which can be invested at rate $r$
- Waiting until maturity to receive $K - S_T$ means forgoing interest on $K$, 
while the stock can barely fall further

The cost of waiting is real and quantifiable. The benefit of waiting has nearly 
vanished. Early exercise dominates.

### The General Early Exercise Condition

More generally, think of holding the put as a trade-off between two things:

**The cost of waiting**: every period the option is held unexercised, the holder 
foregoes interest proportional to $r(K-S)$

**The benefit of waiting**: the stock could fall further, increasing the payoff. 
This is the option's remaining time value — the value of continued optionality.

Early exercise is optimal whenever the cost exceeds the benefit. When $S \approx 0$, the remaining optionality collapses to zero and this condition 
is trivially satisfied. But the same crossover occurs more broadly — whenever the 
put is sufficiently deep in-the-money, interest rates are high, or volatility is 
low enough that the time value of waiting no longer justifies the interest cost.

---

## The Black–Scholes PDE Perspective

The Black–Scholes framework gives a precise, formal account of why early exercise 
becomes optimal. Under the Black–Scholes assumptions, any derivative $P(S,t)$ on 
a non-dividend-paying stock satisfies the PDE:

<div>
\[
\underbrace{\frac{\partial P}{\partial t}}_{\text{Time Decay}}
+ \underbrace{\frac{1}{2}\sigma^2 S^2 \frac{\partial^2 P}{\partial S^2}}_{\text{Convexity Gain (Gamma)}}
+ \underbrace{rS\frac{\partial P}{\partial S}}_{\text{Drift of Underlying}}
- \underbrace{rP}_{\text{Carry Cost}}
= 0
\]
</div>


This PDE is derived by constructing a delta-hedged portfolio and requiring that its 
value grows at the risk-free rate under no-arbitrage. Each term has a financial meaning:

- **Time Decay** $\frac{\partial P}{\partial t}$: the rate at which the option loses 
value as expiry approaches, holding $S$ fixed
- **Convexity Gain** $\frac{1}{2}\sigma^2 S^2 \frac{\partial^2 P}{\partial S^2}$: 
the gain from being long gamma — because the put is convex in $S$, the holder 
benefits on average from large moves in either direction
- **Drift of Underlying** $rS\frac{\partial P}{\partial S}$: since 
$\frac{\partial P}{\partial S} < 0$ for a put, the risk-neutral upward drift of 
the stock works against the put holder
- **Carry Cost** $-rP$: the opportunity cost of holding the option rather than 
investing its value at the risk-free rate

The PDE is a **balance condition**: the convexity gain from being long gamma exactly 
offsets the combined drag from time decay, adverse drift, and carry cost.

### The Hold Region and the Exercise Region

For a European put, the PDE holds everywhere — the holder has no choice but 
to wait. For an American put, the holder has agency, and the picture splits into 
two regions.

In the **hold (continuation) region**, where $P(S,t) > K - S$, the option is worth 
more alive than dead. The convexity gain is sufficient to justify the carry cost and 
the drag from drift. The PDE holds with **equality**:

$$\frac{\partial P}{\partial t} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 P}{\partial S^2} + rS\frac{\partial P}{\partial S} - rP = 0$$

In the **exercise (stopping) region**, where $P(S,t) = K - S$, the balance breaks 
down. Substituting $P = K - S$:

- $\frac{\partial P}{\partial t} = 0$ — intrinsic value has no time decay
- $\frac{1}{2}\sigma^2 S^2 \frac{\partial^2 P}{\partial S^2} = 0$ — intrinsic value 
is linear in $S$, so gamma is zero
- $rS\frac{\partial P}{\partial S} = -rS$ — the delta of $K - S$ is $-1$
- $-rP = -r(K-S)$

The PDE evaluates to:

$$0 + 0 + (-rS) - r(K - S) = -rK < 0$$

The equality strictly fails. With gamma zero and drift negligible, the carry cost 
$rK$ is entirely uncompensated. The PDE inequality signals that continuation is 
dominated by immediate exercise.

### The $S \to 0$ Case: A Formal Illustration

To see this most clearly, assume $S \approx 0$ and suppose the holder continues to 
hold. As $S \to 0$, the gamma and drift terms vanish:

$$\frac{1}{2}\sigma^2 S^2 \frac{\partial^2 P}{\partial S^2} \to 0, \qquad rS\frac{\partial P}{\partial S} \to 0$$

The PDE reduces to:

$$\frac{\partial P}{\partial t} - rP = 0 \implies \frac{\partial P}{\partial t} \approx rK$$

This says the option's value must be **growing** at rate $rK$ per unit time. But 
that is impossible — when $S \approx 0$, the put is already worth approximately $K$, 
its maximum possible value. There is no room left to grow. The assumption of 
holding cannot be sustained, and the option must be exercised.

Note that $S \to 0$ is an extreme case that makes the argument unambiguous. The 
same logic applies more broadly: whenever the put is sufficiently deep in-the-money, 
interest rates are high, or volatility is low enough that the remaining optionality 
has eroded below the the opportunity cost of delaying the receipt of $K$, early exercise is optimal.

---

## Practical Implications

This result has direct consequences for practitioners:

**Options pricing**: American puts must be priced using methods that account for 
early exercise — binomial trees, finite difference methods, or approximation 
formulas. Using Black–Scholes directly will 
systematically underprice them, with the error growing as the put goes deeper 
in-the-money or as interest rates rise.

**Hedging**: The delta of an American put in the stopping region is $-1$ — the 
option moves one-for-one with the stock. A hedger treating it as a live option 
with a partial delta will be systematically underhedged.

**Structured products**: Any product with embedded American put optionality requires careful treatment of the early exercise boundary. Ignoring it introduces 
model risk that can be material in high rate environments.

---