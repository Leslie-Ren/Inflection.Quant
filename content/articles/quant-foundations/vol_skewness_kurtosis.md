---
title: "From Option Prices to the Shape of Returns: A Model-Free Construction of Volatility, Skewness and Kurtosis"
date: 2026-04-16
draft: false
math: true
tags: ["Implied Moments", "Volatility", "Skewness", "Kurtosis", "Static Replication", "Model-Free Methods"]
---
## Why This Matters

Most of my early intuition about options came from the Black-Scholes model, which is 
clean and widely used. But once I started working with real option data, it becomes 
clear that the Black-Scholes assumption of a lognormal distribution is too restrictive. 
For a given maturity, the implied volatility is not constant across strikes, and its 
shape suggests asymmetry and heavy tails in the risk-neutral distribution.

That leads to a more basic question. Instead of imposing a parametric model and 
calibrating its parameters, is there a way to extract volatility, skewness, and kurtosis 
directly from option prices in a model-free way? This is where the Bakshi, Kapadia, and 
Madan (2003) framework becomes useful. Their key idea is that smooth payoff functions 
can be represented as a continuum of vanilla options across strikes. In this view, 
volatility, skewness, and kurtosis are not model assumptions or calibration outputs. 
They are quantities that can be recovered from market prices through static option 
replication.

This article is my attempt to explain the math behind their ideas. But before getting 
into their formulas, I first revisit a simple form of Taylor's theorem, which turns out 
to be the key foundation behind their representation.


## Taylor's Theorem with Integral Remainder

Let $H$ be a twice continuously differentiable function and let $\bar{S}$ be a reference 
point. The first-order Taylor expansion with integral remainder states:

$$H(S) = H(\bar{S}) + H'(\bar{S})(S - \bar{S}) + \int_{\bar{S}}^{S} H''(K)(S - K)dK 
\tag{1}$$

**Derivation.** Start with the fundamental theorem of calculus:

$$H(S) = H(\bar{S}) + \int_{\bar{S}}^{S} H'(t)dt$$

Apply the fundamental theorem of calculus again to $H'(t)$:

$$H'(t) = H'(\bar{S}) + \int_{\bar{S}}^{t} H''(u)du$$

Substitute into the first equation:

$$H(S) = H(\bar{S}) + \int_{\bar{S}}^{S} \left[ H'(\bar{S}) + \int_{\bar{S}}^{t} H''(u)du 
\right] dt$$

Since $H'(\bar{S})$ is constant with respect to $t$, this separates into:

$$H(S) = H(\bar{S}) + H'(\bar{S})(S - \bar{S}) + \int_{\bar{S}}^{S}\int_{\bar{S}}^{t} 
H''(u)dudt$$

It remains to simplify the double integral. The integration runs over the triangular 
region:

$$\{(u, t) : \bar{S} \leq u \leq t \leq S\}$$

Switching the order of integration — integrating over $t$ first, then $u$ — this same 
region is described by $\bar{S} \leq u \leq S$ and $u \leq t \leq S$:

$$\int_{\bar{S}}^{S}\int_{\bar{S}}^{t} H''(u)\,du\,dt = \int_{\bar{S}}^{S}\int_{u}^{S} 
H''(u)dtdu$$

Since $H''(u)$ is constant with respect to $t$, the inner integral evaluates to 
$(S - u)$:

$$\int_{\bar{S}}^{S}\int_{u}^{S} H''(u)\,dt\,du = \int_{\bar{S}}^{S} H''(u)(S - u)du$$

Substituting back and renaming the dummy variable $u \to K$:

$$H(S) = H(\bar{S}) + H'(\bar{S})(S - \bar{S}) + \int_{\bar{S}}^{S} H''(K)(S - K)dK 
\qquad \blacksquare$$

---

## Pricing a Claim via Static Replication

Equation (1) holds for any $S$ and any reference point $\bar{S}$, but the remainder 
integral runs between $\bar{S}$ and $S$, which depends on the realized value of $S$. To 
express the remainder in terms of traded instruments, i.e. calls and puts at fixed strikes, we rewrite it by separating the call and put payoffs. Noting that for any $K$:

$$(S - K)^+ - (K - S)^+ = S - K$$

we can split the single integral over $[\bar{S}, S]$ into an integral of call payoffs over 
all strikes above $\bar{S}$ and an integral of put payoffs over all strikes below $\bar{S}$:

$$H(S) = H(\bar{S}) + H'(\bar{S})(S - \bar{S}) + \int_{\bar{S}}^{\infty} H''(K)(S - K)^+ dK + \int_0^{\bar{S}} H''(K)(K - S)^+dK \tag{2}$$

This works because when $S > \bar{S}$, the term $(S - K)^+$ is zero for all $K > S$, so 
the first integral effectively only collects contributions from $K \in [\bar{S}, S]$, 
which matches the original remainder in equation (1). The second integral vanishes 
entirely since $(K - S)^+ = 0$ for $K \leq \bar{S} < S$. The case $S \leq \bar{S}$ is 
symmetric: only the put integral contributes.

Equation (2) decomposes any smooth payoff into three components: a bond position of size 
$H(\bar{S}) - \bar{S} H'(\bar{S})$, a stock position of size $H'(\bar{S})$, and a 
continuum of calls and puts at every strike, each weighted by the second derivative 
$H''(K)$ of the payoff evaluated at that strike.

Now set $\bar{S} = S_0$, where $S_0$ is the known stock price at the pricing date and 
$S_t$ is the unknown stock price at maturity. This departs slightly from BKM's original 
notation, where they use $S_t$ for the known current price and $S$ for the unknown 
future price. I find the $S_0 / S_t$ convention cleaner for distinguishing constants 
from random variables. 

Applying risk-neutral valuation to both sides of equation (2) gives the arbitrage-free price of the claim:

$$e^{-rt}E^Q[H(S_t)]=[H(S_0) - S_0 H'(S_0)]e^{-rt} + H'(S_0)S_0 + \int_{S_0}^{\infty} H''(K)C(0,K)dK+ \int_0^{S_0} H''(K)P(0,K)dK \tag{3}$$


where $C(0,K)$ and $P(0,K)$ are the time-$0$ prices of European calls and puts with 
strike $K$ and maturity $t$. The bond and stock terms follow from the facts that 
$H(S_0)$ and $H'(S_0)$ are constants at time $0$, and that the no-arbitrage forward 
price satisfies $e^{-rt}\mathbb{E}^Q[S_t] = S_0$.[^1]
[^1]: This relation holds for non-dividend-paying assets. For dividend-paying equities or
indices such as the S&P 500, the no-arbitrage forward price is $S_0 e^{(r-q)t}$ where
$q$ is the continuous dividend yield. In that case, the stock and bond terms in equation
(3) require adjustment. Throughout this article we maintain the no-dividend assumption
for clarity; the extension is straightforward.

## Payoff Functions for the Volatility, Cubic, and Quartic Contracts

Define the log return over the horizon $t$ as:

$$R_t = \ln\frac{S_t}{S_0}$$

BKM propose three contracts whose payoffs are defined by powers of this return:

| Contract | Payoff $H(S_t)$ |
|---|---|
| Volatility | $\left(\ln \dfrac{S_t}{S_0}\right)^2$ |
| Cubic | $\left(\ln \dfrac{S_t}{S_0}\right)^3$ |
| Quartic | $\left(\ln \dfrac{S_t}{S_0}\right)^4$ |

The cubic contract captures asymmetry in the return distribution and the quartic 
contract captures tail heaviness. Their risk-neutral discounted values are:

$$V(0,t) = e^{-rt}\mathbb{E}^Q\left[R_t^2\right], \quad 
W(0,t) = e^{-rt}\mathbb{E}^Q\left[R_t^3\right], \quad 
X(0,t) = e^{-rt}\mathbb{E}^Q\left[R_t^4\right]$$

By equation (3), each can be replicated by a static portfolio of options once we compute $H'(S_0)$ and $H''(K)$ for each payoff and substitute into the replication formula. That is the task of the 
next section.

---
## Computing the Option Weights

Replicating each contract requires two ingredients: the first derivative 
$H'(S_0)$ evaluated at the current stock price, and the second derivative $H''(K)$ 
evaluated at each strike $K$. The first derivative determines the stock position and the 
second derivative determines the weight on each option. We compute these for all three 
contracts in turn and show a graph of the option weights at the end of the section. I plan to discuss the intuition of the option weights profile in a separate article.

### 1. The Volatility Contract

$$H(S_t) = \left(\ln\frac{S_t}{S_0}\right)^2$$

Taking derivatives with respect to $S_t$, evaluated at strike $K$:

$$H'(S_t) = \frac{2\ln\frac{S_t}{S_0}}{S_t}, \qquad H''(K) = \frac{2 - 2\ln\frac{K}{S_0}}{K^2} = \frac{2\left(1 - \ln\frac{K}{S_0}\right)}{K^2}$$

At $S_t = S_0$, the log term vanishes: $H'(S_0) = 0$. This means the stock position 
is zero and the bond position $H(S_0) - S_0 H'(S_0) = 0$ vanishes as well. The 
replication formula reduces entirely to the option integrals:

$$V(0,t) = \int_{S_0}^{\infty} \frac{2\left(1 - \ln\frac{K}{S_0}\right)}{K^2} C(0,K)dK + \int_0^{S_0} \frac{2\left(1 + \ln\frac{S_0}{K}\right)}{K^2} P(0,K)dK$$

where for the put integral we used the fact that $-\ln(K/S_0) = \ln(S_0/K)$ for
$K < S_0$. The weight on OTM calls, $2(1 - \ln(K/S_0))/K^2$, is positive for
near-the-money strikes but turns negative for strikes above $eS_0$, approximately 2.7
times the current stock price. In practice this sign change has little numerical
consequence: option prices at such extreme strikes are close to zero, and the $K^2$
denominator suppresses the contribution further. The put weight $2(1 + \ln(S_0/K))/K^2$
is always positive. Overall, the volatility contract is long calls and puts across
nearly all practically relevant strikes.

### 2. The Cubic Contract

$$H(S_t) = \left(\ln\frac{S_t}{S_0}\right)^3$$

Taking derivatives:

$$H'(S_t) = \frac{3\left(\ln\frac{S_t}{S_0}\right)^2}{S_t}, \qquad H''(K) = 
\frac{6\ln\frac{K}{S_0} - 3\left(\ln\frac{K}{S_0}\right)^2}{K^2}$$

At $S_t = S_0$, the log term again vanishes: $H'(S_0) = 0$, so the stock and bond 
positions are both zero. The replication formula is:

$$W(0,t) = \int_{S_0}^{\infty} \frac{6\ln\frac{K}{S_0} - 3\left(\ln\frac{K}{S_0}
\right)^2}{K^2} C(0,K)\,dK - \int_0^{S_0} \frac{6\ln\frac{S_0}{K} + 3\left(\ln\frac{S_0}
{K}\right)^2}{K^2} P(0,K)\,dK$$

The sign pattern here is economically meaningful. The cubic contract is long calls and 
short puts. When the risk-neutral distribution is left-skewed, OTM puts are expensive 
relative to OTM calls, so the cost of the short put position exceeds the long call 
position, driving the cubic contract value, $W$, to negative. A more negative $W$ corresponds to a more left-skewed distribution.

### 3. The Quartic Contract

$$H(S_t) = \left(\ln\frac{S_t}{S_0}\right)^4$$

Taking derivatives:

$$H'(S_t) = \frac{4\left(\ln\frac{S_t}{S_0}\right)^3}{S_t}, \qquad H''(K) = 
\frac{12\left(\ln\frac{K}{S_0}\right)^2 - 4\left(\ln\frac{K}{S_0}\right)^3}{K^2}$$

At $S_t = S_0$: $H'(S_0) = 0$, so again the stock and bond positions vanish. The 
replication formula is:

$$X(0,t) = \int_{S_0}^{\infty} \frac{12\left(\ln\frac{K}{S_0}\right)^2 - 4\left(
\ln\frac{K}{S_0}\right)^3}{K^2} C(0,K)\,dK + \int_0^{S_0} \frac{12\left(\ln\frac{S_0}
{K}\right)^2 + 4\left(\ln\frac{S_0}{K}\right)^3}{K^2} P(0,K)\,dK$$

The quartic contract is long both calls and puts. Unlike the volatility contract, however, the weights grow with distance from $S_0$: deep out-of-the-money options receive 
progressively larger weights. This makes the quartic contract especially sensitive to 
tail options, which is exactly what we want from a kurtosis measure.

### Summary of Option Weights

It is useful to collect the results in a single table. Define $m = \ln(K/S_0)$ for 
calls ($K > S_0$) and $m = \ln(S_0/K)$ for puts ($K < S_0$), so $m > 0$ in both cases.

| Contract | Weight on OTM calls | Weight on OTM puts | Direction |
|---|---|---|---|
| Volatility | $\dfrac{2(1-m)}{K^2}$ | $\dfrac{2(1+m)}{K^2}$ | Long calls and puts |
| Cubic | $\dfrac{6m - 3m^2}{K^2}$ | $-\dfrac{6m + 3m^2}{K^2}$ | Long calls, short puts |
| Quartic | $\dfrac{12m^2 - 4m^3}{K^2}$ | $\dfrac{12m^2 + 4m^3}{K^2}$ | Long calls and puts |


{{< bkm_weights>}}
## Relating Contract Prices to Risk-Neutral Moments

We now have the three contract prices $V(0,t)$, $W(0,t)$, and $X(0,t)$, each recoverable 
from observed option prices. The final step is to assemble these into the familiar 
statistical moments: variance, skewness, and kurtosis of the risk-neutral return 
distribution.

### The Risk-Neutral Mean

Before computing the centered moments, we need the risk-neutral mean log return:

$$\mu_t = \mathbb{E}^Q[R_t]$$

$\mu_t$ can be recovered from option prices using exactly the same static replication approach as the other three contracts — simply apply equation (3) to the payoff $$H(S_t) = \ln(S_t/S_0) = R_t$$. 
The result is:

$$\mu_t = e^{rt} - 1 - e^{rt}\int_{S_0}^{\infty} \frac{C(0,K)}{K^2}dK - 
e^{rt}\int_0^{S_0} \frac{P(0,K)}{K^2}dK$$

So $\mu_t$ is fully determined by observed option prices with no model assumptions, 
just like the volatility, cubic, and quartic contracts.

With $\mu_t$ in hand, the centered moments follow from standard moment algebra. Write 
$\hat{R}_t = R_t - \mu_t$ for the demeaned return.

**Variance.** The second centered moment is:

$$\mathbb{E}^Q\left[\hat{R}_t^2\right] = \mathbb{E}^Q\left[R_t^2\right] - \mu_t^2 
= e^{rt}V(0,t) - \mu_t^2$$

**Third centered moment.** Expanding $(R_t - \mu_t)^3$ and taking expectations under 
$\mathbb{Q}$:

$$\mathbb{E}^Q\left[\hat{R}_t^3\right] = \mathbb{E}^Q\left[R_t^3\right] - 
3\mu_t\mathbb{E}^Q\left[R_t^2\right] + 2\mu_t^3 = e^{rt}W(0,t) - 
3\mu_t e^{rt}V(0,t) + 2\mu_t^3$$

**Fourth centered moment.** Expanding $(R_t - \mu_t)^4$ and taking expectations under 
$\mathbb{Q}$:

$$\mathbb{E}^Q\left[\hat{R}_t^4\right] = \mathbb{E}^Q\left[R_t^4\right] - 
4\mu_t\mathbb{E}^Q\left[R_t^3\right] + 6\mu_t^2\mathbb{E}^Q\left[R_t^2\right] - 
3\mu_t^4 = e^{rt}X(0,t) - 4\mu_t e^{rt}W(0,t) + 6\mu_t^2 e^{rt}V(0,t) - 3\mu_t^4$$

### The Moment Formulas

Combining the above, the risk-neutral variance, skewness, and kurtosis are:

$$\text{Var}^Q_t = e^{rt}V(0,t) - \mu_t^2$$

$$\text{SKEW}^Q_t = \frac{e^{rt}W(0,t) - 3\mu_t e^{rt}V(0,t) + 2\mu_t^3}
{\left(e^{rt}V(0,t) - \mu_t^2\right)^{3/2}}$$

$$\text{KURT}^Q_t = \frac{e^{rt}X(0,t) - 4\mu_t e^{rt}W(0,t) + 
6\mu_t^2 e^{rt}V(0,t) - 3\mu_t^4}{\left(e^{rt}V(0,t) - \mu_t^2\right)^{2}}$$

Every quantity on the right hand side — $V$, $W$, $X$, and $\mu_t$ — is recoverable from the OTM option prices via static replication. No model has been assumed beyond the existence of a risk-neutral measure. In a typical implementation, the integrals are approximated numerically using a discrete set of observed option prices across available strikes. 

In practice, risk-neutral skewness extracted for an equity or equity index (e.g. SPX) is typically negative, reflecting the market’s tendency to assign higher probability and higher price to large downside moves than to symmetric upside moves. In intuitive terms, skewness captures crash asymmetry: downside moves are sharper, more abrupt, and more expensive to insure than upside moves.

Risk-neutral kurtosis is typically elevated relative to a normal distribution, capturing the market’s expectation that extreme moves, both positive and negative, occur more frequently than Gaussian assumptions would imply. In this sense, kurtosis measures the frequency and severity of extreme outcomes, independent of direction.

Together, these imply that the risk-neutral distribution is left-skewed and fat-tailed, with asymmetric crash risk and elevated probability of extreme events relative to the lognormal benchmark of Black–Scholes.
 
 ## Conclusion
 The BKM framework shows that variance, skewness, and kurtosis of the risk-neutral distribution can be recovered directly from vanilla option prices without calibrating a parametric model. The key is Taylor's theorem with integral remainder, which decomposes any smooth payoff into a static portfolio of calls and puts weighted by the payoff's second derivative. In practice, the continuum of strikes is replaced by a discrete set of observed option prices, requiring numerical integration across available strikes. 

 ## References

- Bakshi, G., Kapadia, N., & Madan, D. (2003). Stock return characteristics, skew laws, and the differential pricing of individual equity options. *The Review of Financial Studies*, 16(1), 101–143.
- Carr, P., & Madan, D. (2001). Optimal positioning in derivative securities. *Quantitative Finance*, 1(1), 19–37.
