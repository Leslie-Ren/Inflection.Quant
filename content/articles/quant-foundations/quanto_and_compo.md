---
title: "Quanto and Compo Commodity Options: FX's Hidden Role in Pricing and Risk"
date: 2026-05-19
draft: false
math: true
tags: ["options", "commodities", "FX", "structured-products", "risk-management", "WTI", "Black-Scholes"]
---
## Why This Matters
 
Many of the world's most actively traded commodities are priced in USD, yet end investors and corporates often operate in other currencies. A Canadian oil producer hedging output, a European airline managing jet fuel costs, or an Asian sovereign wealth fund allocating to commodity exposure all face the same underlying issue: commodity risk does not exist in isolation from FX risk. The standard approach is to hedge the commodity leg with USD-denominated futures or swaps and manage FX separately through forwards or options. This works, but it treats the two risks as independent. Quanto and compo options take a different approach by packaging both risks into a single instrument, but the way each handles FX risk creates some pricing and hedging subtleties that I find are easy to miss.

This article works through both structures using WTI/CAD as an example and focuses on a few practical questions:
 
> - Why does the quanto adjustment exist and what drives its magnitude?
> - How do FX volatility and correlation enter the pricing of each structure, and why do they enter differently?
> - When does each structure suit a given participant?
> - How do dealers hedge each structure, and which risks are hardest to manage?

---
 
## Payoff Structures
 
We work with two assets throughout this article. Let F denote the USD price of WTI crude futures, and let X denote the spot USDCAD exchange rate, quoted as Canadian dollars per one US dollar. A call option struck at K has payoff at expiry T as follows.
 
**Quanto Call**
 
$$\text{Quanto Payoff} = \bar{X} \cdot \max(F_T - K, 0)$$
 
where $\bar{X}$ is a fixed exchange rate agreed at inception. The buyer receives the intrinsic value of the oil option, converted at a predetermined rate, regardless of where USDCAD trades at expiry. The FX rate is contractually frozen.
 
**Compo Call**
 
$$\text{Compo Payoff} = \max(X_T F_T - K, 0)$$
 
where K is denominated in CAD. Here the oil price is first converted to CAD at the prevailing spot rate $X_T$, and the option is exercised based on whether that CAD-denominated oil price exceeds the CAD strike. The option is in the money only if the full currency-converted price clears the hurdle. Both oil moves and FX moves determine whether and by how much the option pays.
 
The distinction is obvious. In the quanto, the buyer has pure oil exposure with zero residual FX risk: FX only affects the fixed conversion notional. In the compo, the strike itself is in CAD, so a weakening USD can push the option out of the money even if oil rises in USD terms, and a strengthening USD can push the option into the money even if oil is flat. The exercise decision and the payoff magnitude are both affected by FX. The compo is inherently a two-dimensional product.
 
---
 
## Dynamics of F and X Under the CAD Risk-Neutral Measure
 
To price these instruments consistently we must work under a single measure. Since payoffs are denominated in CAD, we use the CAD risk-neutral measure $\mathbb{Q}^{CAD}$.
 
Let $r_d$ be the CAD interest rate and $r_f$ be the USD interest rate. Under $\mathbb{Q}^{CAD}$:
 
$$\frac{dX}{X} = (r_d - r_f)\,dt + \sigma_X\,dW_X^{CAD}$$
 
This is standard: the drift of USDCAD under the domestic (CAD) measure equals the interest rate differential, and $\sigma_X$ is the FX volatility.
 
For WTI futures, the natural starting point is the USD risk-neutral measure $\mathbb{Q}^{USD}$. Under $\mathbb{Q}^{USD}$, futures are martingales by no-arbitrage, so F is driftless:
 
$$\frac{dF}{F} = \sigma_F\,dW_F^{USD}$$
 
The two Brownian motions $W_F^{USD}$ and $W_X^{USD}$ have instantaneous correlation $\rho$:

$$dW_F^{USD}\,dW_X^{USD} = \rho\,dt$$

Since X is quoted as CAD per USD, a rally in oil that strengthens CAD causes USDCAD to fall. This is not incidental: Canada is one of the world's largest oil exporters, and CAD is widely regarded as a petrocurrency whose value is closely tied to energy prices. Oil returns and USDCAD returns therefore move in opposite directions systematically, giving $\rho < 0$ for this pair.

Since measure change only adds a drift and leaves quadratic covariation unchanged, $\rho$ is invariant under the measure change. The same correlation holds under $\mathbb{Q}^{CAD}$:

$$dW_F^{CAD}\,dW_X^{CAD} = \rho\,dt$$
 
Since our payoffs are denominated in CAD, we need to express $F$ under $\mathbb{Q}^{CAD}$ rather than $\mathbb{Q}^{USD}$. The key result is that under $\mathbb{Q}^{CAD}$, $F$ is no longer driftless but acquires a drift of $-\rho\,\sigma_F\,\sigma_X$, known as the quanto adjustment. This drift arises entirely from the correlation between $F$ and $X$ and would vanish if the two were independent.

The derivation below shows how this drift emerges from the Radon-Nikodym derivative linking the two measures. Readers comfortable with the result can skip ahead to the next section.

<details>
<summary>Derivation: measure change from $\mathbb{Q}^{USD}$ to $\mathbb{Q}^{CAD}$</summary>

To move from $\mathbb{Q}^{USD}$ to $\mathbb{Q}^{CAD}$ we need the Radon-Nikodym derivative linking the two measures (See [this article]({{< relref "change_of_numeraire.md" >}}) for the general change-of-numeraire framework). The USD and CAD risk-neutral measures are both obtained by discounting with their respective money market accounts, and the exchange rate X connects them. The Radon-Nikodym derivative is proportional to the ratio of the USD and CAD numeraires expressed in a common currency:

$$\frac{d\mathbb{Q}^{USD}}{d\mathbb{Q}^{CAD}}\bigg|_T = \frac{X_T / X_0}{e^{(r_d - r_f)T}}$$

This is the value at time T of one unit of USDCAD forward, normalized to start at 1. Since X follows geometric Brownian motion under $\mathbb{Q}^{CAD}$, we can write this Radon-Nikodym derivative explicitly as:

$$\left.\frac{d\mathbb{Q}^{USD}}{d\mathbb{Q}^{CAD}}\right|_T = \exp\!\left(-\frac{1}{2}\sigma_X^2 T + \sigma_X W_X^{CAD}(T)\right)$$

This is the standard Girsanov density for a constant shift $\sigma_X$. The Radon-Nikodym derivative is driven entirely by $W_X^{CAD}$, the Brownian motion of the exchange rate process. By Girsanov's theorem, $W_X^{CAD}$ acquires a drift when we move to $\mathbb{Q}^{USD}$:

$$dW_X^{USD} = dW_X^{CAD} - \sigma_X\, dt$$

**Why F Acquires a Drift**

The Brownian motions $W_F^{USD}$ and $W_X^{USD}$ have instantaneous correlation $\rho$, meaning we can always decompose:

$$dW_F^{USD} = \rho\, dW_X^{USD} + \sqrt{1-\rho^2}\, dW_\perp$$

where $W_\perp$ is a Brownian motion independent of $W_X^{USD}$. Applying the measure change $dW_X^{USD} = dW_X^{CAD} - \sigma_X\,dt$:

$$dW_F^{USD} = \rho(dW_X^{CAD} - \sigma_X\,dt) + \sqrt{1-\rho^2}\,dW_\perp = -\rho\,\sigma_X\,dt + \rho\,dW_X^{CAD} + \sqrt{1-\rho^2}\,dW_\perp$$

The diffusion part $\rho\,dW_X^{CAD} + \sqrt{1-\rho^2}\,dW_\perp$ has instantaneous correlation $\rho$ with $W_X^{CAD}$, which is exactly the definition of $dW_F^{CAD}$, so we write:

$$dW_F^{USD} = -\rho\,\sigma_X\,dt + dW_F^{CAD}$$

Substituting into the SDE for $F$:

$$\frac{dF}{F} = \sigma_F\,dW_F^{USD} = \sigma_F\!\left(dW_F^{CAD} - \rho\,\sigma_X\,dt\right)$$

Under $\mathbb{Q}^{CAD}$, $F$ therefore follows:

$$\frac{dF}{F} = -\rho\,\sigma_F\,\sigma_X\,dt + \sigma_F\,dW_F^{CAD}$$

The drift $-\rho\,\sigma_F\,\sigma_X$ is the quanto adjustment. It is a direct consequence of changing the measure of a correlated asset.

</details>
 
---
 
## The Quanto Adjustment: Why It Exists

Although the quanto adjustment $-\rho\,\sigma_F\,\sigma_X$ emerges naturally from the measure change derivation, it is worth building an intuition for why it exists and why its magnitude takes the form it does. The hedging argument provides a clean economic explanation.

Consider a bank that has sold a quanto forward to a client: at maturity the bank pays $F_T \cdot \bar{X}$ in CAD, where $\bar{X}$ is the fixed exchange rate agreed at inception. To hedge, the bank goes long a regular WTI futures contract and converts the USD proceeds at the market rate $X_T$ at maturity, receiving $F_T \cdot X_T$ in CAD.

The bank's hedging P&L is:

$$F_T \cdot X_T - F_T \cdot \bar{X} = F_T(X_T - \bar{X})$$

Setting $\bar{X} = \mathbb{E}[X_T]$ to simplify the illustration[^1], the expected hedging cost becomes:

$$\mathbb{E}[F_T(X_T - \bar{X})] = \mathbb{E}[F_T X_T] - \mathbb{E}[F_T]\,\mathbb{E}[X_T] = \text{Cov}(F_T, X_T)$$

Since $\rho < 0$ for this pair, the covariance is negative: the simple hedge bleeds in expectation. To break even, the bank must charge the client a forward price above $F_0$. The required markup is determined by the covariance $\text{Cov}(F_T, X_T) = \rho \sigma_F \sigma_X T$

In continuous time, the bleeding accumulates at rate $\rho\,\sigma_F\,\sigma_X$ per unit time proportionally to the current level of $F$, giving the SDE under $\mathbb{Q}^{CAD}$:

$$\frac{dF}{F} = -\rho\,\sigma_F\,\sigma_X\,dt + \sigma_F\,dW_F^{CAD}$$

Compounding this proportional drift over $T$ gives the quanto-adjusted forward, which we denote $F_0^*$:

$$F_0^* \equiv \mathbb{E}^{\mathbb{Q}^{CAD}}[F_T] = F_0 \cdot e^{-\rho\,\sigma_F\,\sigma_X\,T}$$

Since $\rho < 0$, the exponent is positive and $F_0^* > F_0$. The more negative $\rho$ is, and the larger $\sigma_F$ and $\sigma_X$ are, the more the simple hedge bleeds and the greater the adjustment required.

---
[^1]: Without this simplification, the expected hedging cost contains an additional term $\mathbb{E}[F_T](\mathbb{E}[X_T] - \bar{X})$. This term is purely deterministic since both $\mathbb{E}[F_T]$ and $\mathbb{E}[X_T]$ are known at inception and $\bar{X}$ is a fixed contractual rate, and it vanishes when $\bar{X} = \mathbb{E}[X_T]$. The covariance term is the only irreducible source of hedging cost and is the true quanto effect.

## Pricing the Quanto Option
 
The quanto call eliminates FX risk by fixing the conversion rate at $\bar{X}$. The payoff is linear in F alone. Under $\mathbb{Q}^{CAD}$, we need to price:
 
$$V_{quanto} = \bar{X} \cdot e^{-r_d T} \cdot \mathbb{E}^{CAD}\left[\max(F_T - K, 0)\right]$$
 
Since F under $\mathbb{Q}^{CAD}$ is lognormal with forward $F_0^*$ as derived above, applying the Black formula directly:
 
$$V_{quanto} = \bar{X} \cdot e^{-r_d T} \cdot \left[F_0^*\,N(d_1) - K\,N(d_2)\right]$$
 
where:
 
$$d_1 = \frac{\ln(F_0^*/K) + \frac{1}{2}\sigma_F^2\,T}{\sigma_F\sqrt{T}}, \qquad d_2 = d_1 - \sigma_F\sqrt{T}$$
 
This is simply a Black formula with the oil future price replaced by $F_0^*$. The vol input is $\sigma_F$ alone: FX volatility enters only through the correlation term absorbed into $F_0^*$, and a larger negative covariance results in a higher adjusted forward and a higher call value.
 
---
 
## Pricing the Compo Option

The compo payoff is $\max(X_T F_T - K, 0)$ where K is in CAD. The key observation is that $X_T F_T$ is itself a lognormal under $\mathbb{Q}^{CAD}$, since it is the product of two correlated lognormals. We define the CAD-denominated oil price:

$$S_T = X_T F_T$$

The SDE for $S_T$ follows from Itô's lemma applied to the product of $X_T$ and $F_T$. Under $\mathbb{Q}^{CAD}$:

$$\frac{dS}{S} = (r_d - r_f)\,dt + \sigma_X\,dW_X^{CAD} + \sigma_F\,dW_F^{CAD}$$

The cross term $\rho\sigma_X\sigma_Fdt$ from Itô's lemma exactly cancels the quanto drift. As a result, $S_T$ drifts at $(r_d - r_f)$. This reflects that $S_t$ is a USD-denominated commodity price expressed in CAD units, whose drift is governed by the relative pricing of USD and CAD under the CAD measure. The diffusion is driven jointly by $W_X^{CAD}$ and $W_F^{CAD}$.

The compo call is then simply a call on $S_T$ struck at K, all in CAD:

$$V_{compo} = e^{-r_d T} \cdot \mathbb{E}^{\mathbb{Q}^{CAD}}\left[\max(S_T - K, 0)\right]$$

To apply Black's formula we need the forward and the volatility of $S_T$. 

$$S_0^{fwd} = F_0 \cdot X_0\,e^{(r_d - r_f)T}$$

$$\sigma_{compo} = \sqrt{\sigma_F^2 + 2\rho\,\sigma_F\,\sigma_X + \sigma_X^2}$$

Applying Black's formula directly to $S_T$:

$$V_{compo} = e^{-r_d T}\left[S_0^{fwd}\,N(d_1^c) - K\,N(d_2^c)\right]$$

where:

$$d_1^c = \frac{\ln(S_0^{fwd}/K) + \frac{1}{2}\sigma_{compo}^2\, T}{\sigma_{compo}\sqrt{T}}, \qquad d_2^c = d_1^c - \sigma_{compo}\sqrt{T}$$

The compo is a standard Black call on the CAD-denominated oil forward, with a composite vol that blends oil vol, FX vol, and their covariance. FX volatility enters quadratically, unlike through the drift as in the quanto. 

While it may seem that a higher FX vol always increases the compo vol, this is not always the case. Taking the derivative of $\sigma_{compo}$ with respect to $\sigma_X$:

$$\frac{\partial\,\sigma_{compo}}{\partial\,\sigma_X} = \frac{\sigma_X + \rho\,\sigma_F}{\sigma_{compo}}$$

This is positive only when $\sigma_X > -\rho\,\sigma_F$. When $\rho$ is negative, increasing $\sigma_X$ can decrease $\sigma_{compo}$ because the negative cross term $2\rho\,\sigma_F\,\sigma_X$ grows in magnitude faster than the $\sigma_X^2$ term. Intuitively, when $\rho < 0$, oil and FX move in opposite directions, and large FX moves increasingly offset the oil moves in CAD terms. In an extreme case of large $\sigma_X$ and very negative $\rho$, the FX leg almost perfectly hedges the oil leg and $S_T$ barely moves at all. The compo option can therefore be cheaper than a plain oil option when correlation is sufficiently negative, reflecting the natural hedge between oil and CAD that was discussed in the quanto adjustment section.
 
---
 
# Which Structure Suits Which Participant?

### Currency of Exposure

The most fundamental question is what currency the participant's exposure actually lives in. A Canadian producer or refiner whose budget constraint is a CAD breakeven price is asking "is oil above C$\$130$?" rather than "is oil above \$100 USD?". For that participant the compo is the natural fit: the strike is set directly in their decision currency and the exercise decision aligns with their actual P&L. The quanto can hedge the same oil exposure but the exercise is made in USD terms, introducing a mismatch against a CAD budget that the participant must then manage separately. Conversely, a fund reporting in CAD whose mandate is pure commodity exposure benefits from the quanto's fixed conversion rate, which removes USD/CAD as a P&L variable entirely and keeps the oil attribution clean.

### Relative Cost

Once the currency question is settled, cost becomes the next consideration, and neither structure is universally cheaper. At $\rho = 0$ the compo tends to be more expensive because it embeds FX risk directly into the payoff, raising the total volatility $\sigma_{compo}$. But as $\rho$ becomes more negative, the cross term in $\sigma_{compo}$ works in the buyer's favour, and for WTI/CAD where $\rho$ is meaningfully negative and $\sigma_F$ is substantially larger than $\sigma_X$, the compo can be cheaper than the quanto. The crossover point where the compo premium falls below the quanto premium depends on the interplay between the compo vol reduction and the quanto forward adjustment $F_0^*$, and participants who are indifferent between the two payoff structures should price both under current market inputs before deciding.
 

### Operational Complexity
Operational simplicity favours the quanto, particularly for corporate treasuries and smaller counterparties. Both structures require $\sigma_F$, $\sigma_X$, and $\rho$, none of which are directly observable. But in the quanto, $\rho$ enters only through the drift adjustment in $F_0^*$, and once that adjusted forward is computed the valuation reduces to a standard single-underlying Black formula. In the compo, $\rho$ enters $\sigma_{compo}$ directly and the sensitivity of the option value to correlation is more immediate and material, making independent marking harder for a treasury without a dedicated quant function. Beyond valuation, the cross-gamma and correlation risks discussed in the next section mean dealers may charge a wider bid-offer spread on the compo, which partially offsets any premium saving from the lower composite vol and should be factored into the all-in cost comparison.

The pricer below allows direct comparison of both structures under user-specified inputs. The hedging section that follows explains how dealers manage each once the trade is on.


{{<compo_quanto_pricer>}}

---
 
## Hedging Each Structure
 
Understanding which structure fits a given exposure is only half the picture. Once a dealer has sold either instrument, the more operationally demanding question is how to manage the risk through the life of the trade. The two structures present meaningfully different hedging problems.
 
**Hedging the Quanto**
 
The quanto call has only one source of price risk: the level of WTI. Because the exchange
rate is fixed contractually, USD/CAD is not a risk factor and the dealer runs no FX delta.
The delta hedge is therefore a straightforward position in WTI futures, sized to the
Black delta evaluated at the quanto-adjusted forward.
As WTI moves, the futures position is rebalanced in the usual way.
 
Although there is no FX delta, a USD/CAD forward is still required for currency
translation. The WTI futures position generates P&L in USD while the liability to the
option holder is in CAD. The dealer enters a USD/CAD forward sized to the expected dollar
value of the futures position to convert those proceeds into CAD at a known rate. This
forward is updated as the delta is rebalanced. It is a funding hedge rather than a
risk-factor hedge. It does not arise from any sensitivity of the option value to the
exchange rate, but from the operational mismatch between the currency of the hedging
instrument and the currency of the liability.
 
The more subtle risk in the quanto is correlation between WTI returns and USD/CAD returns.
Correlation enters the pricing formula through the drift of the oil forward under the CAD
measure, and the dealer who sold the option carries residual exposure to shifts in this
parameter. This is difficult to hedge because correlation is not directly traded. Pure
correlation products such as covariance swaps exist but are often illiquid in the commoidty/FX market. In practice dealers manage correlation exposure within
book limits, accepting that residual exposure will sit on the book as a managed risk.
The sensitivity is relatively contained, however, because correlation enters only through
the drift adjustment rather than through the volatility of the underlying itself.

**Hedging the Compo**
 
The compo is more complex because both WTI and USD/CAD are live risk factors. The payoff
depends on the product $F \cdot X$, so the dealer must run two delta hedges simultaneously:
a WTI futures position and a USD/CAD forward position. The two hedges are coupled: the
WTI delta is proportional to the current spot FX rate, and the FX delta is proportional
to the current WTI forward. Every time one underlying moves, the hedge notional in the
other leg must be updated. This cross-gamma cannot be fully eliminated with vanilla
instruments, and is typically treated as a cost of carry priced into the bid-offer spread
at inception.
 
On the volatility side, both WTI vol and FX vol enter $\sigma_{compo}$, so the dealer
carries independent vega in each. WTI vega is hedged with WTI options, and FX vega with
USD/CAD options. As discussed in the pricing section, the sign of FX vega depends on
$\sigma_X + \rho\sigma_F$. For WTI/CAD where $\rho < 0$, this quantity can be
negative, meaning the dealer who sold the option is long rather than short FX vega. Dealers
must verify this sign before structuring the USD/CAD options overlay, as the
positive-correlation intuition would produce a hedge in the wrong direction.
 
Correlation risk in the compo is more material than in the quanto because $\rho$ enters
$\sigma_{compo}$ directly rather than just the drift. The instinctive view is that a dealer
who sold a call is short correlation, since higher $\rho$ increases $\sigma_{compo}$ and
raises the option value. But this reasoning imports a positive-correlation assumption
silently. For WTI/CAD where $\rho$ is negative, a move toward more positive correlation
does hurt the dealer, while a further decline in correlation reduces $\sigma_{compo}$ and
benefits them. The direction of the exposure is not fixed: it depends on where $\rho$
currently sits and which way it moves. As with the quanto,
residual correlation risk is carried on the book and priced into the spread at inception.


**Comparative Summary**
 
| | Quanto | Compo |
|---|---|---|
| WTI delta | WTI futures at $F_0^*$ | WTI futures; notional scales with spot FX |
| FX hedge | USD/CAD forward for P&L translation only; not a risk-factor hedge | USD/CAD forward as risk-factor delta hedge; notional scales with WTI forward |
| Cross-gamma | None | Cannot be fully hedged with vanilla instruments; treated as cost of carry and priced into the spread |
| WTI vega | WTI options | WTI options |
| FX vega | Enters only via $F_0^*$; less material than WTI vega for short-dated trades, but grows with tenor and correlation  | USD/CAD options; sign of exposure determined by $\sigma_X + \rho\sigma_F$ |
| Correlation risk | Via drift; relatively contained | Via $\sigma_{compo}$; more material; sign depends on level of $\rho$ |
| Correlation hedge | Residual book risk; priced into spread | Residual book risk; priced into spread |
| Overall complexity | Moderate | Higher; two coupled dynamic hedges |
 
---


## A Note on Calibrating $\rho$

Throughout this article, $\rho$ appears in every formula and drives many of the key risk management decisions. In practice, calibrating it is less straightforward than calibrating $\sigma_F$ or $\sigma_X$, both of which can be implied from liquid option markets. Correlation has no directly quoted instrument. Dealers typically estimate $\rho$ from historical return series, often using rolling windows of varying length and weighting schemes that emphasise recent observations. The choice of window, frequency, and whether to use spot or futures returns can produce meaningfully different estimates. Some desks supplement historical estimates with implied correlation backed out from traded cross-asset products where available, though for WTI/CAD such instruments are sparse. The resulting uncertainty in $\rho$ is itself a source of model risk, and given how directly it enters both the quanto drift adjustment and the compo composite volatility, even modest miscalibration can shift prices and hedge ratios materially. This is one of the reasons dealers price correlation risk into the spread rather than attempting to hedge it precisely.