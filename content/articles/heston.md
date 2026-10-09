---
title: "Heston Model: What the Square Root Bought and What It Cost"
date: 2026-09-30
draft: false
math: true
tags: ["stochastic volatility", "Heston", "CIR", "Feller condition", "characteristic function", "Fourier pricing", "forward start options", "barrier options", "Monte Carlo", "local volatility"]
ShowToc: true
TocOpen: false
---

## Why This Matters

The local vol model we discussed in the [previous article]({{< ref "local_vol.md#linking-densities-to-option-quotes" >}}) assumes the vol is a deterministic function of spot and time, so the vol has no randomness of its own. But the market does not treat vol that way. The VIX is the market's implied vol for SPX, and options on the VIX are themselves quoted with an implied vol, a vol of vol. The market is pricing vol itself as a random quantity.

One consequence of ignoring the randomness in the vol is that the forward vol skew comes out wrong: it flattens as the start date moves out, so a forward starting option is mispriced, as we explained in the local vol article. A common fix is the stochastic vol model, which treats vol as a second stochastic process. Heston is among the most commonly used models of this kind. When I first learned about Heston, it was not obvious to me why we landed on that particular form, so before stating its assumptions and working through its pricing, I want to first trace the earlier models that led to Heston.

## The Road to Heston

### The easiest way to make the variance random

Black-Scholes assumes the spot is lognormal with a constant vol, and we ask what the simplest change is to make the variance random. The answer is to recycle the same geometric Brownian motion we already have for the spot and apply it to the variance $v_t = \sigma_t^2$, which is the Hull and White (1987) model:

$$dS_t = \mu S_t\,dt + \sqrt{v_t}\,S_t\,dW_t^S$$
$$dv_t = \mu_v v_t\,dt + \xi v_t\,dW_t^v$$
$$d\langle W^S, W^v\rangle_t = \rho\,dt$$

The variance stays positive, and it moves in proportion to a shock because of the $\xi v_t\,dW_t^v$ term. Since the variance is lognormal, the vol is also lognormal and carries the same proportional scaling. A 100% shock doubles the vol whether it starts at 10% or at 30%.

Setting $\rho = 0$ gives a useful conditioning result shared by stochastic vol models:

$$C = \mathbb{E}\left[C_{BS}(\bar v)\right], \qquad \bar v = \frac{1}{T}\int_0^T v_t\,dt.$$

Conditional on the path of $v$, the vol is no longer random and the spot is lognormal, so Black-Scholes applies. The path enters only through its average variance $\bar v$, giving each path a Black-Scholes price at its own $\bar v$, and the option price averages those over the distribution of $\bar v$. Since $\bar v$ is a time average of a lognormal process, its density function has no closed form. To compute $C$, Hull and White expand $C_{BS}(\bar v)$ in a Taylor series around $\mathbb{E}[\bar v]$. The price becomes a one dimensional expectation and can be handled easily by numerical methods.

Two things about this model need improvement. The variance does not mean revert, so nothing holds it near a long run level. That is not how a vol market behaves, where a high vol regime decays and a quiet one eventually picks up.

The second is that the conditioning argument above needs $\rho = 0$, and that assumption removes the vol skew. With independent drivers, a high variance path is just as likely when the spot rises as when it falls, so the smile it produces is symmetric: neither tail is favored. In equity markets $\rho$ is typically negative. A fall in the spot comes with a rise in variance, so once the spot is down it becomes more volatile and can end far below its starting level. A rise in the spot comes with a smaller variance, so once it is up it becomes calmer and tends to stay close. This puts more mass far below the forward and less mass far above it, giving the fat left tail and the negative skew the market quotes. To keep the vanilla price tractable, the Hull-White model assumes zero correlation, so it cannot produce the negative skew observed in equity markets.

### Making the vol mean revert

Write the Hull-White vol in log form via Ito's lemma,

$$d\log\sigma_t = m\,dt + \tfrac{\xi}{2}\,dW_t^v, \qquad m = \tfrac{\mu_v}{2} - \tfrac{\xi^2}{4},$$

so the log vol is a Brownian motion with a constant drift, which is where the problem sits: the drift does not depend on where the log vol currently is. To make it mean revert, we write

$$d\log\sigma_t = \kappa(\theta - \log\sigma_t)\,dt + \tfrac{\xi}{2}\,dW_t^v$$

This is the Ornstein-Uhlenbeck process, a canonical mean reverting process, and putting it on the log vol is a specification considered by Scott (1987) and Wiggins (1987). $\theta$ is the long run mean and $\kappa$ is the speed of mean reversion. The gap between the current log vol and $\theta$ is expected to decay exponentially at rate $\kappa$:

$$\mathbb{E}\left[\log\sigma_t - \theta\right] = \left(\log\sigma_0 - \theta\right)e^{-\kappa t}.$$

Take $\theta = \log 0.20$ and $\kappa = 2$. The vol is pulled back toward 20% and the gap halves every four months, since $e^{-2 \times 4/12} \approx 0.51$.

The disadvantage of this model is that even with $\rho = 0$ the vanilla price has no tractable form. It still reduces to an average of Black-Scholes prices over $\bar v$, but that average cannot be evaluated in closed form. Prices can still be computed by Monte Carlo or by solving the PDE numerically, but a calibration runs hundreds of those solves inside an optimizer, which is too slow to be practical.

### Getting a tractable price

By now we have the mean reversion but are stuck with an intractable model. Stein and Stein (1991) worked with the Ornstein-Uhlenbeck process on the vol itself, which gives a model that is priceable at speed:

$$d\sigma_t = \kappa(\theta - \sigma_t)\,dt + \xi\,dW_t^v.$$

Now $\sigma_t$ is Gaussian, so the integrated variance $\int_0^T \sigma_t^2\,dt$ is a quadratic functional of a Gaussian process. This model does not have a closed form solution for the vanilla option, but it is quick to solve via Fourier transform. Stein and Stein derived it for the $\rho = 0$ case. Schöbel and Zhu (1999) later derived the pricing through the characteristic function of the log spot, which admits $\rho \neq 0$. We will not go over the details here but will do a deep dive of the same Fourier mechanism later in the Heston model.

What this model gives up is positivity: the vol can cross zero and go negative. Since only $\sigma_t^2$ enters the log spot's variance and drift, this is not a problem for the spot dynamics, but it is still a problem for the vol. A vol is a standard deviation, so a negative value has no meaning by itself. A subtler issue is what happens when the vol is close to zero. The drift is negligible against the diffusion over a short interval $\Delta t$: the drift moves the vol by an amount of order $\Delta t$, while the diffusion moves it by order $\sqrt{\Delta t}$, which is much larger. So the vol crosses zero back and forth many times in quick succession. These zero variance episodes come in dense bursts, so the model manufactures artificial calm periods.

### Keeping the vol positive

Now we have a dilemma. Both models mean revert, but the log OU guarantees positivity and loses a tractable option price, while the OU is solvable but allows the vol to go negative. That is not a problem unique to vol. When OU was used to model the short rate, its ability to go negative was seen at the time as a defect,[^negrates] and the fix was to scale the shock by the square root of the state,

$$dr_t = \kappa(\theta - r_t)\,dt + \xi\sqrt{r_t}\,dW_t.$$

This is the Cox, Ingersoll and Ross (1985) process. It solves the dilemma cleanly: the mean reversion survives, the process never goes negative, and the bond price comes out in closed form. We will explain where the price tractability comes from in the Heston section. Whether the rate can reach zero or stays strictly positive is set by the Feller condition (Feller, 1951), which gives strict positivity when

$$2\kappa\theta \ge \xi^2.$$

The idea is that the process touches zero only if it takes a downward shock larger than its current distance from the boundary. Near zero, the OU and CIR processes have the same upward pull, a drift close to $\kappa\theta > 0$. The OU keeps the shock size fixed, so even with a strong pull, a large enough shock eventually carries the process through zero. CIR shrinks the shock as the distance narrows, so near zero the pull faces an ever smaller shock. If the pull is strong enough, the process never reaches zero. The Feller condition is derived below, which readers already familiar with it can skip.

<details style="margin-bottom: 1.5em;">
<summary>Where the Feller condition comes from</summary>

We are essentially trying to find the probability that the rate ever reaches zero, given the current rate $r_0 = r$,

$$u(r) = P(\text{the rate ever reaches } 0 \mid r_0 = r).$$

A probability may appear difficult to solve for. But we can convert it into an expectation of a payoff and view it as an option pricing problem, which we are quite familiar with: $u(r)$ is what we would pay for a claim that returns a dollar if the rate ever reaches zero,

$$u(r) = \mathbb{E}\big[\mathbf{1}\{\text{the rate reaches } 0\} \mid r_0 = r\big].$$

The value of a claim on a diffusion satisfies the [backward Kolmogorov equation]({{< ref "forward_backward_pde.md" >}}). For CIR that reads

$$\partial_t u + \kappa(\theta - r)\,\partial_r u + \tfrac{1}{2}\xi^2 r\,\partial_{rr} u = 0.$$

The payoff is a pure probability with no deadline attached, so the value does not depend on the calendar: starting from a given rate, the chance of eventually reaching zero is the same today as a year from now. That makes $\partial_t u = 0$, and the equation drops to an ODE in $r$:

$$\tfrac{1}{2}\xi^2 r\,u''(r) + \kappa(\theta - r)\,u'(r) = 0.$$

There is no $u$ term, so this is a separable first order equation in $u'$, and solving it gives

$$u'(r) = C\,r^{-q}\,e^{2\kappa r/\xi^2}, \qquad q = \frac{2\kappa\theta}{\xi^2}.$$

Integrating once more recovers $u$ itself,

$$u(r) = A + C\int^r s^{-q}e^{2\kappa s/\xi^2}\,ds,$$

with a second constant $A$. We need two conditions to fix both constants. One is readily available: $u(0) = 1$, since a rate already at zero has reached it. For the other we impose an artificial condition $u(b) = 0$ at a high level $b$, which measures the probability of reaching zero before rising to $b$, and then let $b \to \infty$ to recover the probability of ever reaching zero. With $u(0) = 1$ and $u(b) = 0$,

$$u(r) = \frac{\int_r^b s^{-q}e^{2\kappa s/\xi^2}\,ds}{\int_0^b s^{-q}e^{2\kappa s/\xi^2}\,ds}.$$

When $q < 1$, both numerator and denominator are finite, and $u(r)$ is a probability strictly between zero and one: there is a positive chance the rate reaches zero before $b$. As $b \to \infty$ the two integrals grow without bound and their ratio tends to one, so the rate is certain to reach zero eventually. When $q \ge 1$, that is $2\kappa\theta \ge \xi^2$, the numerator is still finite but the denominator diverges, since near zero the integrand behaves like $s^{-q}$, so $u(r) = 0$ and zero is not reached before $b$. This holds for every $b$, so it holds in the limit as well, and the rate never reaches zero at any time. The threshold $2\kappa\theta \ge \xi^2$ is the Feller condition.

</details>

CIR also has an interesting property: two independent CIR processes with the same $\kappa$ and $\xi$ add up to another CIR process. Take

$$dr^{(i)}_t = \kappa(\theta_i - r^{(i)}_t)\,dt + \xi\sqrt{r^{(i)}_t}\,dW^{(i)}_t, \qquad i = 1, 2,$$

with $W^{(1)}$ and $W^{(2)}$ independent, and let $r_t = r^{(1)}_t + r^{(2)}_t$. The drifts add directly,

$$dr_t = \kappa(\theta_1 + \theta_2 - r_t)\,dt + \xi\left(\sqrt{r^{(1)}_t}\,dW^{(1)}_t + \sqrt{r^{(2)}_t}\,dW^{(2)}_t\right).$$

The noise term does not look like a single square root yet, but its quadratic variation is $\xi^2(r^{(1)}_t + r^{(2)}_t)\,dt = \xi^2 r_t\,dt$. So

$$dW_t = \frac{\sqrt{r^{(1)}_t}\,dW^{(1)}_t + \sqrt{r^{(2)}_t}\,dW^{(2)}_t}{\sqrt{r_t}}$$

is a continuous martingale with $d\langle W\rangle_t = dt$. Therefore

$$dr_t = \kappa(\theta_1 + \theta_2 - r_t)\,dt + \xi\sqrt{r_t}\,dW_t.$$

The sum is again a CIR process, and its long run level is simply $\theta_1 + \theta_2$. This reminded me of the chi-square distribution, which is also closed under addition: adding two independent chi-squares gives another chi-square, and their degrees of freedom add up. And it turns out that CIR indeed follows a scaled noncentral chi-square distribution with $\delta$ degrees of freedom and noncentrality $\lambda$,

$$r_t = c\,Y, \qquad Y \sim \chi'^2(\delta, \lambda),$$

where

$$c = \frac{\xi^2(1 - e^{-\kappa t})}{4\kappa}, \qquad \delta = \frac{4\kappa\theta}{\xi^2}, \qquad \lambda = \frac{r_0 e^{-\kappa t}}{c}.$$

The degrees of freedom $\delta$ are proportional to $\theta$, so adding levels adds degrees of freedom. Since the noncentral chi-square has a known density, $r_t$ also has a closed form density by the chain rule, $f_{r_t}(r) = \frac{1}{c}\,f_Y(r/c)$.

<details style="margin-bottom: 1.5em;">
<summary>Building CIR from squared Gaussian processes</summary>

A chi-square is a sum of squared normals, so we look for independent processes $X_1, \dots, X_n$ that are normal at every $t$ and whose squares sum to a CIR process. The standard Gaussian diffusion is the linear SDE

$$dX_i = (b_i + a_i X_i)\,dt + \sigma_i\,dW_i.$$

Set $Y = \sum_i X_i^2$. By Ito,

$$dY = \sum_i \left(2b_i X_i + 2a_i X_i^2 + \sigma_i^2\right)dt + \sum_i 2\sigma_i X_i\,dW_i,$$

and we want this to match $dY = \kappa(\theta - Y)\,dt + \xi\sqrt{Y}\,dW$. Four conditions follow.

1. The drift must depend on $Y$ alone. The term $2b_i X_i$ depends on the sign of $X_i$, not only on $X_i^2$, so $b_i = 0$ for all $i$.
2. The remaining drift must give the mean reversion, $\sum_i 2a_i X_i^2 = -\kappa\sum_i X_i^2$ for every value of the $X_i$. That forces $a_i = -\kappa/2$ for all $i$.
3. The noise has quadratic variation $\sum_i 4\sigma_i^2 X_i^2\,dt$, which must equal $\xi^2 Y\,dt$. That forces $\sigma_i = \xi/2$ for all $i$, and then, by the same argument as for the sum above, the noise can be written as $\xi\sqrt{Y}\,dW$.
4. What is left of the drift is the constant $\sum_i \sigma_i^2 = n\xi^2/4$, which comes from squaring the noise, and it plays the role of $\kappa\theta$. That fixes the number of components, $n = 4\kappa\theta/\xi^2 = \delta$.

So the components are identical OU processes, $dX_i = -\tfrac{\kappa}{2}X_i\,dt + \tfrac{\xi}{2}\,dW_i$, and there are $\delta$ of them. This construction directly explains the distribution when $\delta$ is an integer; the noncentral chi-square law extends the result to arbitrary positive $\delta$.

</details>

### Getting all three at once

Heston (1993) takes the CIR process from short rates, applies it to the variance, and lets the two Brownian motions correlate:

$$dS_t = \mu S_t\,dt + \sqrt{v_t}\,S_t\,dW_t^S$$
$$dv_t = \kappa(\theta - v_t)\,dt + \xi\sqrt{v_t}\,dW_t^v$$
$$d\langle W^S, W^v\rangle_t = \rho\,dt$$

The result is a model with all three properties at once. The drift $\kappa(\theta - v_t)$ mean reverts, the $\sqrt{v_t}$ diffusion keeps the variance positive, strictly so under Feller, and the vanilla option price stays tractable even with $\rho \neq 0$. Each earlier model was missing one of the three. Hull-White prices quickly and stays positive, but does not mean revert, and its pricing argument needs the two Brownian motions to be independent, which is what leaves its smile symmetric. The log OU mean reverts and stays positive, but has no fast price, so its correlation, though unrestricted, cannot be calibrated in practice. Stein-Stein mean reverts and prices quickly, and Schöbel and Zhu later recovered $\rho \neq 0$ for it, but its vol can cross zero and go negative, which no market vol does. Heston was the first widely influential model to combine all three.

| Model | Variance process | Mean reverts | Stays positive | Tractable vanilla option price | Skew, which requires $\rho \neq 0$ |
|---|---|---|---|---|---|
| Hull-White (1987) | GBM on $v$ | ✗ | ✓ | ✓ (Taylor expansion, but only if $\rho = 0$) | ✗ symmetric smile, since the tractable price needs $\rho = 0$ |
| Scott, Wiggins (1987) | OU on $\log\sigma$ | ✓ | ✓ | ✗ (Monte Carlo or PDE) | ✓ $\rho \neq 0$ allowed, but without a tractable price |
| Stein-Stein (1991) | OU on $\sigma$ | ✓ | ✗ (vol crosses zero) | ✓ (Fourier) | ✓ $\rho \neq 0$ via Schöbel-Zhu (1999) |
| **Heston (1993)** | CIR on $v$ | ✓ | ✓ (strictly, under Feller) | ✓ (Fourier) | ✓ $\rho \neq 0$ allowed, with the same tractable price |

## Where the Tractability Comes From

A model whose vanilla price needs a PDE solve or a Monte Carlo run is time consuming to work with. Calibrating its parameters means computing prices many times over, and the cost adds up. This is what held back the log OU model. Heston has an advantage: its vanilla option price is tractable and faster to evaluate than either PDE or Monte Carlo, and this is what turned it into a practical reference model. We spend this section on where that tractability comes from.

The price of a call is the discounted payoff integrated against the density of the terminal log spot $X_T = \log S_T$,

$$C = e^{-rT}\int_{-\infty}^{\infty} (e^{x} - K)^{+} p(x)\,dx.$$

Under Black-Scholes $p$ is the Gaussian density and the integral evaluates to the Black-Scholes formula. Under Heston, conditional on the variance path, $X_T$ is still Gaussian, with variance $(1 - \rho^2)\int_0^T v_t\,dt$ and a mean that depends on the path through $\rho\int_0^T \sqrt{v_t}\,dW_t^v$. Removing the conditioning means averaging the conditional Gaussian over the joint distribution of those two path integrals, and the mixture it produces has no closed form. The density route to the option price is blocked.

### The Fourier Route to the Price

Notice the pricing integral has only two components, the payoff and the density. Although the density is out of reach, we can work on the payoff instead. We turn to a powerful tool discussed in an earlier article, the [Fourier transform]({{< ref "fourier_transform.md" >}}). The idea is to write the payoff as an integral of exponentials. Pricing the payoff then comes down to taking the expectation of an exponential of the log spot, which is its characteristic function. For many models the characteristic function has a closed form even though the density does not, and we will see that Heston is one of them.

First we want to reconstruct the payoff from exponential frequencies, the standard inverse Fourier transform. But the transform is only defined for functions that are integrable, meaning $\int_{-\infty}^{\infty}|f(x)|\,dx$ is finite. The call payoff grows like $e^{x}$, so it is not integrable. The way around this is to multiply by a factor that damps the payoff before transforming it, an idea that goes back to Carr and Madan (1999):

$$g(x) = e^{-\alpha x}\left(e^{x} - K\right)^{+}, \qquad \alpha > 1.$$

The damped payoff decays like $e^{-(\alpha - 1)x}$, so it is integrable and its transform exists. Projecting it onto a specific exponential frequency $u$ gives the weight that frequency carries,

$$\hat{g}(u) = \int_{-\infty}^{\infty} e^{-iux}g(x)\,dx = \frac{K^{1-s}}{s(s-1)}, \qquad s = \alpha + iu.$$

Note that $\hat{g}$ is only a function of $\alpha$, the strike and the frequency, so it can be reused for every parameter set during pricing and calibration.

Reconstructing $g$ from its frequencies and undoing the damping recovers the payoff itself,

$$\left(e^{x} - K\right)^{+} = \frac{e^{\alpha x}}{2\pi}\int_{-\infty}^{\infty}\hat{g}(u)\,e^{iux}\,du.$$

Substituting this into the pricing integral leaves a double integral in $x$ and $u$, and exchanging the order of integration puts the $x$ integral on the inside,

$$C = \frac{e^{-rT}}{2\pi}\int_{-\infty}^{\infty}\hat{g}(u)\left[\int_{-\infty}^{\infty} e^{(\alpha + iu)x}p(x)\,dx\right]du.$$

The inner integral is the density paired with a single exponential, and that is the characteristic function

$$\phi(z) = \mathbb{E}\left[e^{izX_T}\right] = \int_{-\infty}^{\infty} e^{izx}p(x)\,dx.$$

Matching the exponents, $iz = \alpha + iu$, so it is evaluated at $z = u - i\alpha$. The price becomes a single integral of the contract's transform against the model's characteristic function,

$$C = \frac{e^{-rT}}{2\pi}\int_{-\infty}^{\infty}\hat{g}(u)\,\phi(u - i\alpha)\,du.$$

If Heston has a closed form characteristic function, the call is valued by numerical integration of a single integral, and we never need the density at all. Let us work out that characteristic function next.

### The Backbone Without the Vol of Vol

The difficulty in computing the characteristic function comes from the randomness of the variance, which makes the log spot non-Gaussian, so we can simplify the problem by switching that randomness off first. With $\xi = 0$ the variance follows its mean reverting drift deterministically,

$$v_t = \theta + (v_0 - \theta)e^{-\kappa t}.$$

The spot still diffuses, but with a known variance at every instant, so the log spot is Gaussian. Under the risk neutral measure used for pricing, the spot drifts at $r$ in place of $\mu$:

$$X_T \sim \mathcal{N}\left(X_0 + rT - \tfrac{1}{2}V_T,\; V_T\right), \qquad V_T = \int_0^T v_t\,dt = \theta\left(T - \frac{1 - e^{-\kappa T}}{\kappa}\right) + v_0\,\frac{1 - e^{-\kappa T}}{\kappa}.$$

Note $V_T$ is a constant plus $v_0$ times a factor that depends only on $T$. A function of this form, linear plus a constant, is called affine, so $V_T$ is affine in $v_0$.

The characteristic function of a Gaussian is known,

$$\phi(u) = \exp\left(iu(X_0 + rT) - \tfrac{1}{2}(u^2 + iu)V_T\right).$$

Substituting $V_T$, the characteristic function also reduces to an exponential affine in $v_0$,

$$\phi(u) = \exp\left(A(T) + B(T)v_0 + iuX_0\right),$$

$$A(T) = iurT - \frac{u^2 + iu}{2}\,\theta\left(T - \frac{1 - e^{-\kappa T}}{\kappa}\right), \qquad B(T) = -\frac{u^2 + iu}{2}\cdot\frac{1 - e^{-\kappa T}}{\kappa}.$$

Assuming the vol is deterministic, we end up with a closed form characteristic function in exponential affine form. This is only a special case of Heston, but it gives us a reasonable guess for what the characteristic function could look like in the general case, when we switch the vol's randomness back on.

### Putting the Vol of Vol Back

The characteristic function is an expectation of a function of the terminal log spot. Seen from an earlier time $t$, its value depends on the current state, the log spot $x$ and the variance $v$, and on the time left to maturity, $\tau = T - t$,

$$\phi(\tau, x, v) = \mathbb{E}\left[e^{iuX_T} \,\middle|\, X_t = x,\ v_t = v\right].$$

An expectation of this kind must satisfy the [backward Kolmogorov PDE]({{< ref "forward_backward_pde.md" >}}) in the two state variables,

$$\partial_\tau \phi = \left(r - \frac{v}{2}\right)\partial_x \phi + \frac{v}{2}\,\partial_{xx}\phi + \kappa(\theta - v)\,\partial_v \phi + \frac{\xi^2 v}{2}\,\partial_{vv}\phi + \rho\xi v\,\partial_{xv}\phi.$$

At maturity, $\tau = 0$, the log spot $x$ is known, so $\phi$ is just the payoff $e^{iux}$. Writing the equation in $\tau$ turns this payoff into an initial condition. The characteristic function seen from today is $\phi$ at $\tau = T$, $x = X_0$ and $v = v_0$, and the backbone version from the previous section is the case $\xi = 0$. We now check whether the PDE admits a solution of the same exponential affine form, with $A$ and $B$ as functions of $\tau$:

$$\phi = \exp\left(A(\tau) + B(\tau)v + iux\right), \qquad A(0) = B(0) = 0.$$

Every derivative of $\phi$ returns $\phi$ times a factor:

$$\partial_\tau \phi = (A' + B'v)\,\phi, \qquad \partial_x \phi = iu\,\phi, \qquad \partial_{xx} \phi = -u^2\,\phi,$$

$$\partial_v \phi = B\,\phi, \qquad \partial_{vv} \phi = B^2\,\phi, \qquad \partial_{xv} \phi = iuB\,\phi.$$

Substituting these into the PDE and dividing through by $\phi$ leaves

$$A'(\tau) + B'(\tau)\,v = \big(iur + \kappa\theta B\big) + \left(\frac{\xi^2}{2}B^2 - (\kappa - \rho\xi iu)B - \frac{u^2 + iu}{2}\right)v.$$

Both sides are affine in $v$. Since this has to hold at every level of the variance, matching the coefficients gives us two ODEs in $\tau$:

$$A' = iur + \kappa\theta B, \qquad B' = \frac{\xi^2}{2}B^2 - (\kappa - \rho\xi iu)B - \frac{u^2 + iu}{2}.$$

So the exponential affine guess solves the Heston PDE exactly when $A$ and $B$ solve these two ODEs from $A(0) = B(0) = 0$. The second equation is a first order ODE that is quadratic in $B$, known as a Riccati equation, and it has a closed form solution, derived in the [appendix](#appendix-solving-the-riccati-equation). $A$ then follows by integrating $B$. The solutions are

$$A(\tau) = iur\tau + \frac{\kappa\theta}{\xi^2}\left[(b - d)\tau - 2\log\frac{1 - g e^{-d\tau}}{1 - g}\right], \qquad B(\tau) = \frac{b - d}{\xi^2}\cdot\frac{1 - e^{-d\tau}}{1 - g e^{-d\tau}},$$

where

$$b = \kappa - \rho\xi iu, \qquad d = \sqrt{b^2 + \xi^2(u^2 + iu)}, \qquad g = \frac{b - d}{b + d}.$$

### Pricing a Vanilla Option

We now put together everything we have derived so far to price a call option. From the Fourier route, the price is

$$C = \frac{e^{-rT}}{2\pi}\int_{-\infty}^{\infty}\frac{K^{1-s}}{s(s-1)}\,\phi(u - i\alpha)\,du, \qquad s = \alpha + iu, \qquad \alpha > 1,$$

and the characteristic function seen from today is the solution above at $\tau = T$, $x = X_0 = \log S_0$ and $v = v_0$,

$$\phi(u) = \exp\left(A(T) + B(T)v_0 + iuX_0\right).$$

Evaluating it at $u - i\alpha$ means replacing $u$ by $u - i\alpha$ inside $b$, $d$, $g$, $A$ and $B$.

This is a single integral of a smooth function of $u$, which is quick to evaluate by numerical integration. The characteristic function does not depend on the strike, so one set of its values prices every strike at a given maturity, and only the payoff transform changes from strike to strike. The call price does not depend on $\alpha$ in theory. The damping is applied to the payoff before the transform and undone when the payoff is reconstructed. In practice, $\alpha$ affects how smooth the integrand is, and the smoother the curve, the fewer sample points are needed to evaluate the integral accurately. The selection of $\alpha$ will be discussed in a follow up article on Heston calibration.

The numerical computation can be made twice as fast. Replacing $u$ by $-u$ leaves the real part of the integrand unchanged and flips the sign of its imaginary part.[^halfrange] Over the whole range the imaginary part therefore cancels and the real part counts twice, so we only need the real part over positive $u$,

$$C = \frac{e^{-rT}}{\pi}\int_{0}^{\infty}\operatorname{Re}\left[\frac{K^{1-s}}{s(s-1)}\,\phi(u - i\alpha)\right]du.$$

Below is a minimal Heston pricer in Python.

```python
import numpy as np
from scipy.integrate import quad

def heston_phi(z, T, x0, v0, r, kappa, theta, xi, rho):
    b = kappa - rho * xi * 1j * z
    d = np.sqrt(b**2 + xi**2 * (z**2 + 1j * z))
    g = (b - d) / (b + d)
    e = np.exp(-d * T)
    A = 1j * z * r * T + kappa * theta / xi**2 * ((b - d) * T - 2 * np.log((1 - g * e) / (1 - g)))
    B = (b - d) / xi**2 * (1 - e) / (1 - g * e)
    return np.exp(A + B * v0 + 1j * z * x0)

def heston_call(S0, K, T, r, v0, kappa, theta, xi, rho, alpha=1.5):
    x0 = np.log(S0)
    def integrand(u):
        s = alpha + 1j * u
        phi = heston_phi(u - 1j * alpha, T, x0, v0, r, kappa, theta, xi, rho)
        return (K**(1 - s) / (s * (s - 1)) * phi).real
    return np.exp(-r * T) / np.pi * quad(integrand, 0, np.inf, limit=200)[0]

# Heston with kappa = 2, theta = v0 = 0.04, xi = 0.3, rho = -0.7
print(heston_call(100, 100, 1.0, 0.03, 0.04, 2.0, 0.04, 0.3, -0.7))   # 9.2425

# vol of vol close to zero recovers Black-Scholes at 20% vol
print(heston_call(100, 100, 1.0, 0.03, 0.04, 2.0, 0.04, 1e-4, -0.7))  # 9.4134
```

As a check, setting $\xi$ close to zero with $v_0 = \theta = 0.04$ should recover Black-Scholes at 20% vol. With $S_0 = K = 100$, $T = 1$, $r = 3\%$ and $\alpha = 1.5$, the formula gives 9.4134, matching the Black-Scholes price to four decimal places.

### What the Square Root Bought

Heston's characteristic function has a closed form, the exponential affine, and that form is what makes the vanilla option quick to price: the pricing PDE collapses into two ODEs, and the option price becomes a single numerical integral. The **<span style="color:green">square root</span>** bought us the exponential affine, and hence the fast valuation.

$$dv_t = \kappa(\theta - v_t)\,dt + \xi\,{\color{green}\boldsymbol{\sqrt{v_t}}}\,dW_t^v.$$

Under the exponential affine guess, every derivative on the right side of the PDE returns $\phi$ times a factor ($iu$, $-u^2$, $B$, $B^2$, $iuB$). None of these factors contains $v$, because the exponent is linear in $v$. So $v$ enters the right side of the equation only through the five PDE coefficients:

- drift of $x$: $r - v/2$, affine, because the spot's variance is $v$
- variance of $x$: $v$, affine
- drift of $v$: $\kappa(\theta - v)$, affine by design, that is just mean reversion
- variance of $v$: $\xi^2 v$, and this is where the square root shows up, as the diffusion enters the PDE only through its square
- covariance of $x$ and $v$: $\rho\xi v$, the square root again, multiplied against the spot's own $\sqrt{v}$ diffusion so the two half powers combine into a whole one

All five are affine in $v$, and matching against $A' + B'v$ on the left side of the PDE gives the two ODEs. If a vol diffusion is linear in the variance, $\xi v_t\,dW_t^v$, it breaks the last two coefficients: its square is $\xi^2 v^2$ and its covariance with the log return is $\rho\xi v^{3/2}$, and the left side has no $v^2$ or $v^{3/2}$ slot to match.

So the square root is not one of several convenient diffusions: it is the one that lands both the variance of $v$ and the covariance at degree one. In fact, the characteristic function stays exponentially affine exactly when the drift and covariance stay affine in the state (Duffie, Filipović and Schachermayer, 2003). CIR, the process that inspired Heston, belongs to the same family, which is why its bond price comes out in the same exponential affine form. Although the exponential affine delivers the tractability, the reverse is not true: a model does not have to be affine to be tractable. For example, the 3/2 stochastic vol model, with variance diffusion $\xi v_t^{3/2}\,dW_t^v$, has coefficients that are not affine, yet its characteristic function also has a closed form. That comes through a separate route relying on its own special structure, perhaps a topic worth covering in a future article.

## Pricing a Forward Starting Option

With vanilla pricing under Heston in hand, we turn to the product local vol struggles with, the forward starting option. A forward starting call fixes its strike at a future date $T_1$ as a fraction $k$ of the spot on that date, and pays at $T_2$,

$$\left(S_{T_2} - kS_{T_1}\right)^{+}.$$

Today's vanilla prices give the distribution of the spot at each maturity, but this payoff depends on the return from $T_1$ to $T_2$, a window that opens in the future. Its value is set by the smile the market will quote at $T_1$, and that depends on how the spot and the vol evolve jointly between now and then.

### Reducing It to a Vanilla Option

Notice the strike is quoted as a fraction of the spot, which makes the spot the natural unit to price in, the same move we used for the exchange option in the [change of numéraire article]({{< ref "change_of_numeraire.md" >}}). Standing at $T_1$, the strike is known and what remains is a vanilla call with maturity $\tau = T_2 - T_1$, worth $S_{T_1}\,C(k, \tau; v_{T_1})$, where $C(k, \tau; v)$ is the Heston price of a call on a unit spot with strike $k$ and maturity $\tau$, starting from variance $v$. Taking the spot as numéraire simplifies the problem: we no longer have to deal with the correlation between $S_{T_2}$ and $S_{T_1}$, nor with the correlation between $v_{T_1}$ and $S_{T_1}$,

$$P = \mathbb{E}\left[e^{-rT_1}S_{T_1}\,C(k, \tau; v_{T_1})\right] = S_0\,\mathbb{E}^{S}\left[C(k, \tau; v_{T_1})\right],$$

so the forward start option is a mixture of Heston vanillas over the variance at $T_1$, and all that remains is the distribution of $v_{T_1}$ under the share measure. Switching from the risk neutral measure to the stock measure shifts the spot's driver by its own vol, $dW_t^{S,\mathbb{Q}} = dW_t^{S,\mathbb{Q}^S} + \sqrt{v_t}\,dt$, as derived in the [change of numéraire article]({{< ref "change_of_numeraire.md" >}}). Through the correlation between the spot and the variance, this carries into the variance driver and adds $\rho\xi v_t$ to its drift, so under the share measure

$$dv_t = \kappa^*(\theta^* - v_t)\,dt + \xi\sqrt{v_t}\,dW_t^{v,\mathbb{Q}^S}, \qquad \kappa^* = \kappa - \rho\xi, \qquad \theta^* = \frac{\kappa\theta}{\kappa - \rho\xi}.$$

The variance process is still CIR, so its value at $T_1$ is a scaled noncentral chi-square. The measure change moves the long run mean to $\theta^*$ and the reversion speed to $\kappa^*$, but the degrees of freedom $4\kappa^*\theta^*/\xi^2 = 4\kappa\theta/\xi^2$ stay the same. Since this density is known in closed form, the option price becomes a one dimensional integral of the vanilla pricer against the density of $v_{T_1}$. Below is the code to price the forward start call.

```python
from scipy.stats import ncx2

def forward_start_call(S0, k, T1, T2, r, v0, kappa, theta, xi, rho):
    tau = T2 - T1
    ks = kappa - rho * xi                                   # reversion speed under the share measure
    c = xi**2 * (1 - np.exp(-ks * T1)) / (4 * ks)          # scale of the noncentral chi-square
    df = 4 * kappa * theta / xi**2                          # degrees of freedom
    nc = v0 * np.exp(-ks * T1) / c                          # noncentrality
    lo, hi = c * ncx2.ppf([1e-10, 1 - 1e-10], df, nc)
    def integrand(v):
        return heston_call(1.0, k, tau, r, v, kappa, theta, xi, rho) * ncx2.pdf(v / c, df, nc) / c
    return S0 * quad(integrand, lo, hi, limit=200)[0]

# one year forward start, one year option, same parameters as before
print(forward_start_call(100, 1.0, 1.0, 2.0, 0.03, 0.04, 2.0, 0.04, 0.3, -0.7))   # 9.0328
```

As checks, a start date close to zero recovers the vanilla price of 9.2425, a vol of vol close to zero recovers the Black-Scholes price of 9.4134, and a Monte Carlo run with 200,000 paths gives 9.0350 for the ATM case, agreeing to within 0.01.

### What Sets the Forward Smile

From forward start option prices we can back out the implied forward vol and see how the forward smile behaves. The inversion works the same way as for an ordinary implied vol: we set the Heston price equal to the Black-Scholes price, both written with $S_0$ as the numéraire. Unlike Heston, where the variance at $T_1$ has a distribution and the price is an expectation over it, the vol under Black-Scholes is a constant, so no expectation is needed: a forward start is simply worth $S_0$ times a call on a unit spot with strike $k$ and maturity $\tau$.

$$S_0\,\mathbb{E}^{S}\left[C(k, \tau; v_{T_1})\right] = S_0\,C_{BS}(1, k, \tau, r, \sigma_{\text{fwd}}).$$

Let us vary the parameters one at a time and see what each does to the forward smile. The base case is $v_0 = \theta = 0.04$, $\kappa = 2$, $\xi = 0.3$ and $\rho = -0.7$, with $r = 3\%$, and each panel moves one parameter while the rest stay there. The option starts at $T_1 = 1$ and runs for $\tau = 1$ year. The panels vary the risk neutral parameters. The share measure used in the pricing only rewrites the variance drift in terms of $\kappa^*$ and $\theta^*$, and does not change the forward smile itself.

{{< forward_smile_params >}}

**Level, $v_0$ and $\theta$.** The first two panels shift the whole curve up or down. $v_0$ and $\theta$ are the variance now and the variance in the long run, and the forward start sees a blend of the two, weighted by how much of the gap the mean reversion has closed by $T_1$. At a one year start date with $\kappa = 2$ the gap is nearly gone, so $\theta$ moves the whole curve up or down while $v_0$ barely registers. The shorter the start date, the more $v_0$ matters.

**Skew, $\rho\xi$.** With $\rho = 0$ there is no skew: the forward smile is symmetric around the forward, both wings lifted equally. With $\rho < 0$ a falling spot comes with rising variance, which fattens the left tail and tilts the smile, the same mechanism we explained in the [Hull-White section](#the-easiest-way-to-make-the-variance-random). The skew enters through the covariance of the spot and the variance, $d\langle S, v\rangle_t = \rho\xi v_t S_t\,dt$, so it is the product $\rho\xi$ that matters: the correlation sets the sign, and both parameters set the size.

**Curvature, $\xi$ and $\kappa$.** $\xi$, the vol of vol, is what lets the variance wander and $\kappa$, the mean reversion speed, is what pulls it back, so the two work against each other. A larger $\xi$ or a slower $\kappa$ leaves the variance at $T_1$ spread over a wider range, which is why their two panels look so alike. To see how a wider dispersion of the variance raises the curvature, we set the correlation to zero to isolate the impact of the vol of vol. The conditioning result from the [Hull-White section](#the-easiest-way-to-make-the-variance-random) then applies: the option is worth the average of Black-Scholes prices over the paths of the variance. Writing $\bar\sigma = \sqrt{\bar v}$ for a path's average vol,

$$C = \mathbb{E}\left[C_{BS}(\bar\sigma)\right].$$

At the money, the Black-Scholes price is nearly linear in vol, the $0.4\,\sigma F$ rule from the [Brownian motion article]({{< ref "understanding_brownian_motion.md#why-the-atm-black-scholes-formula-reduces-to-a-function-of-sigmasqrtt" >}}), so

$$C \approx \mathbb{E}\left[0.4\,\bar\sigma F\right] = 0.4\,\mathbb{E}[\bar\sigma]\,F.$$

Only the average vol moves the price, and the vol of vol barely matters. Far out of the money, the price is convex in vol,[^wing] and expanding it around the average vol gives

$$C \approx C_{BS}\left(\mathbb{E}[\bar\sigma]\right) + \tfrac{1}{2}\,C_{BS}''\left(\mathbb{E}[\bar\sigma]\right)\operatorname{Var}(\bar\sigma).$$

Since $C_{BS}''$ is positive in the wings, a higher variance of $\bar\sigma$, which is what a higher vol of vol brings, raises the wing option price. A higher implied vol is needed to match it. This lift in the wing implied vols, with the at the money vol almost unchanged, is the curvature.

### Comparison with Local Vol Forward Smiles

The flattening forward skew is one of the main defects of local vol, and it is what sent us to Heston in the first place, so it is natural to compare the forward start implied vols under the two models. To have a fair comparison, the two models must at least agree on the vanilla prices. We can already price those under Heston, and applying the [Dupire formula]({{< ref "local_vol.md#linking-densities-to-option-quotes" >}}) to that surface backs out the local vols. Dupire needs the derivatives of the call price with respect to $T$ and $K$, and both are tractable under Heston.

<details style="margin-bottom: 1.5em;">
<summary>The Dupire derivatives in closed form</summary>

Dupire needs

$$\sigma_{\text{loc}}^2(K, T) = \frac{\partial_T C + rK\,\partial_K C}{\tfrac{1}{2}K^2\,\partial_{KK} C}.$$

In the vanilla option price, the strike sits in one place only, the factor $K^{1-s}$ with $s = \alpha + iu$, so differentiating in $K$ acts on that factor alone.

$$\partial_K C = -\frac{e^{-rT}}{\pi}\int_0^\infty \operatorname{Re}\left[\frac{K^{-s}}{s}\,\phi(u - i\alpha)\right]du, \qquad \partial_{KK} C = \frac{e^{-rT}}{\pi}\int_0^\infty \operatorname{Re}\left[K^{-1-s}\,\phi(u - i\alpha)\right]du.$$

Maturity enters through the discount factor and through $\phi = \exp(A(T) + B(T)v_0 + iuX_0)$, so $\partial_T\phi = (A'(T) + B'(T)v_0)\,\phi$, and $A'$ and $B'$ are given by the two ODEs we already solved,

$$A' = iur + \kappa\theta B, \qquad B' = \frac{\xi^2}{2}B^2 - (\kappa - \rho\xi iu)B - \frac{u^2 + iu}{2},$$

evaluated at $u - i\alpha$. Therefore

$$\partial_T C = -rC + \frac{e^{-rT}}{\pi}\int_0^\infty \operatorname{Re}\left[\frac{K^{1-s}}{s(s-1)}\,\big(A' + B'v_0\big)\,\phi(u - i\alpha)\right]du.$$

All three integrals can be solved the same way as the price, by numerical integration.

</details>

{{< forward_smile_twin >}}

Under local vol the forward smile flattens as the start date moves out, while Heston keeps its skew at every start date. This is because Heston is time invariant. Recall its two SDEs: the drift $\kappa(\theta - v_t)$, the diffusion $\xi\sqrt{v_t}$ and the correlation $\rho$ are all built from constants and the current variance, and the calendar date appears nowhere. So a variance of $v$ behaves the same way whether it is reached today or in two years, and the only thing an option starting at $T_1$ inherits from the intervening period is the level it starts from, which is why the skew survives at every start date. As the start date moves further out, the forward smile settles into a fixed shape instead of flattening.

Local vol is the opposite, since its vol is an explicit function of time. As the [local vol article]({{< ref "local_vol.md#local-vol-matches-the-marginals-not-the-path" >}}) showed, its slices flatten quickly with maturity, so the forward smile is set mainly by the local vols between the start date and maturity, where the skew has largely faded. Under Heston the skew is not inherited from today's surface but rebuilt over each window: the covariance of the spot and the variance, $\rho\xi v_t S_t\,dt$, is driven by the same $\rho\xi$ over $[T_1, T_2]$ as over $[0, T_2 - T_1]$, so the same $\rho\xi$ that tilts today's smile tilts the forward one.

## Pricing a Barrier Option

A barrier option depends not only on the terminal value but also on whether the spot ever touched the barrier along the way. A down and out call pays $(S_T - K)^+$ only if the spot never traded below a barrier $B$ over the option's life. We take the barrier as monitored continuously here, though discrete monitoring is also traded. The Fourier route we used for vanillas is built on the terminal distribution alone, so it no longer works here. We need to turn to numerical methods: PDE and Monte Carlo.

The PDE handles the barrier naturally: at every time step, the value is set to zero once the spot falls below $B$. The greeks come off the same grid at no extra cost. For a single underlying the PDE is decently fast, since its grid only spans two state variables, spot and variance, and marching that grid back through time is fast. But it suffers from the curse of dimensionality. Each additional source of randomness adds a dimension to the grid, and the number of grid points grows exponentially with it. Once there are several underlyings, each carrying its own spot and variance, Monte Carlo takes over as the more practical approach, since its cost grows with the number of paths rather than exponentially with the dimension.

For our example we go with Monte Carlo, as it is more general and easier to demonstrate. The Heston PDE deserves its own article, and I plan to write about it in the future. Throughout this section, we keep the base parameters from the forward start section, $v_0 = \theta = 0.04$, $\kappa = 2$, $\xi = 0.3$, $\rho = -0.7$ and $r = 3\%$.

### When the Simulated Variance Goes Negative

We need to simulate two correlated random variables, the spot and the variance, driven by standard normals $Z_1$ and $Z_2 = \rho Z_1 + \sqrt{1 - \rho^2}\,Z_\perp$. For the spot, work with the log and step it forward over $\Delta t$ using the variance at the start of the step,

$$\log S_{t+\Delta t} = \log S_t + \left(r - \tfrac{1}{2}v_t\right)\Delta t + \sqrt{v_t}\sqrt{\Delta t}\,Z_1.$$

The variance follows the same recipe,

$$v_{t+\Delta t} = v_t + \kappa(\theta - v_t)\Delta t + \xi\sqrt{v_t}\sqrt{\Delta t}\,Z_2.$$

Both updates are Euler steps: each one holds the drift and the shock at their values from the start of the step. Simulating the variance needs some special care because of the $\sqrt{v_t}$ term. When $v_t$ is already close to zero, nothing stops a large negative $Z_2$ from pushing $v_{t+\Delta t}$ below zero, and the next step then asks for the square root of a negative number.

You may wonder whether the Feller condition from the [CIR section](#keeping-the-vol-positive) rules this out, and our base parameters do satisfy it, $2\kappa\theta = 0.16$ against $\xi^2 = 0.09$. It turns out the negative variance is a flaw of the numerical scheme rather than of the Heston model itself. In continuous time, as $v$ heads toward zero the noise $\xi\sqrt{v}$ shrinks along the way, so every small move down is made with a smaller shock, while the drift keeps pushing up at close to $\kappa\theta$. Under the Feller condition the upward push wins, and the variance never reaches zero. The simulation instead fixes the shock for the whole step at $\xi\sqrt{v_t}\sqrt{\Delta t}\,Z_2$, using the variance at the start. If the path is heading down, the true diffusion would be shrinking during the step, but the scheme keeps the larger starting value until the step ends, and since $Z_2$ is unbounded a single draw can carry it straight past zero. This means the scheme goes negative regardless of the parameters, and Feller only changes how often it happens. In a simulation at the base parameters, 1.7% of paths still go negative at 250 steps a year. Calibrated to an equity surface, the Feller condition often fails, and $\xi$ can be near 1 or higher. Rerunning the same simulation at $\xi = 1$, where $\xi^2$ far exceeds $2\kappa\theta = 0.16$, around 95% of paths go negative.

A simple repair is full truncation, with $v_t^+ = \max(v_t, 0)$,

$$v_{t+\Delta t} = v_t + \kappa(\theta - v_t^+)\Delta t + \xi\sqrt{v_t^+}\sqrt{\Delta t}\,Z_2.$$

This works, but it converges slowly, meaning we need small time steps before the Monte Carlo price settles on the true price, and the problem gets worse once Feller fails. Flooring at zero discards the negative draws while keeping the large positive ones, so the spot sees a higher variance on average than it should, and options come out overpriced. At $\xi = 1$, a one year ATM vanilla call whose exact Fourier price is 8.016 comes out at 8.570 with monthly steps and 8.118 with weekly steps, and is still about 0.02 too high at daily steps.

A widely used alternative for simulating Heston is the quadratic exponential scheme of Andersen (2008). Rather than taking an Euler step, it draws the next variance directly from a distribution that resembles the true distribution. Over a step, the exact transition of the variance is the scaled noncentral chi-square, and its mean and variance are known in closed form. The scheme picks a simple positive distribution and matches those two moments. When the dispersion of the next variance is low relative to its level, a squared shifted normal is used to approximate the noncentral chi-square, and when it is high, as happens near zero when the Feller condition fails, an exponential with a point mass at zero is used instead. In both cases the two parameters are solved from the two moments, and the spot is then stepped using the two variance endpoints, which keeps its correlation with the variance intact. Each step of the scheme costs more to compute than an Euler step, but it needs far fewer steps to converge, so the total run time is lower, which is why it is preferred over full truncation in practice.

**Bonus question.** The exact variance transition, the noncentral chi-square distribution, can be sampled directly, so why approximate it? The answer is not speed. Hint: think about how a desk computes greeks.

<details style="margin-bottom: 1.5em;">
<summary>Answer</summary>

A greek is usually computed by bumping a parameter and repricing with the same random numbers, so that the noise in the two prices cancels in their difference. That requires each random number to play the same role in both runs. The usual way to sample a noncentral chi-square relies on rejection sampling, which consumes a varying number of random numbers per draw. A small bump can change whether a candidate is accepted, after which the base and bumped paths run on different draws, their noise no longer cancels, and dividing by a small bump size makes the greek very noisy. The noncentral chi-square can also be sampled by inverting its distribution function, which keeps the mapping, but that inverse has no closed form and needs a numerical root find for every draw.

</details>

### The Monitoring Bias

The barrier is monitored continuously, but a simulation only checks the spot at the grid points. A path that dips below the barrier between two steps and comes back above is still counted as alive, so the simulation keeps paths that should have been knocked out, and the price comes out too high. The bias shrinks only in proportion to $\sqrt{\Delta t}$, because that is how far a path typically moves within one time step. This problem is not unique to Heston. It comes from the discrete grid, and it is present even with a constant vol. We met this problem in the [Monte Carlo variance reduction article]({{< ref "mc_variance_reduction.md#conditional-monte-carlo" >}}), where the fix was the Brownian bridge. Conditional on the log spot at both ends of an interval lying above the barrier, the probability that the path crossed the barrier in between has a closed form under a constant vol assumption,

$$p_{\text{cross}} = \exp\left(-\frac{2\,(x_0 - \log B)(x_T - \log B)}{\sigma^2 T}\right).$$

The same approach can be adapted for Heston by applying it at each time step, with the variance assumed constant over that small interval,

$$p_{\text{cross}} = \exp\left(-\frac{2\,(x_t - \log B)(x_{t+\Delta t} - \log B)}{v_t\,\Delta t}\right).$$

Using the Brownian bridge to correct the monitoring bias, the down and out call prices at 7.454 with 25 steps and 7.451 with 250, against 8.163 and 7.718 without it. Even 1,000 steps without the bridge leave the price at 7.580, still about 0.13 too high.

Below is the Monte Carlo pricer for the down and out call, with full truncation for the variance and a flag to switch the Brownian bridge correction on and off.

```python
import numpy as np

def down_and_out_call(S0, K, B, T, r, v0, kappa, theta, xi, rho,
                      n_steps, n_paths, bridge=True, seed=0):
    rng = np.random.default_rng(seed)
    dt = T / n_steps
    x = np.full(n_paths, np.log(S0))
    v = np.full(n_paths, v0)
    alive = np.ones(n_paths, dtype=bool)
    logB = np.log(B)
    for _ in range(n_steps):
        z1 = rng.standard_normal(n_paths)
        z2 = rho * z1 + np.sqrt(1 - rho**2) * rng.standard_normal(n_paths)
        vp = np.maximum(v, 0.0)                                   # full truncation
        x_new = x + (r - 0.5 * vp) * dt + np.sqrt(vp * dt) * z1
        v = v + kappa * (theta - vp) * dt + xi * np.sqrt(vp * dt) * z2
        alive &= x_new > logB                                     # knocked out at a grid point
        if bridge:                                                # crossed between grid points
            gap = np.maximum(x - logB, 0) * np.maximum(x_new - logB, 0)
            p_cross = np.exp(-2 * gap / np.maximum(vp * dt, 1e-300))
            alive &= rng.random(n_paths) > p_cross
        x = x_new
    payoff = np.where(alive, np.maximum(np.exp(x) - K, 0.0), 0.0)
    return np.exp(-r * T) * payoff.mean(), np.exp(-r * T) * payoff.std() / np.sqrt(n_paths)

# base parameters, S0 = K = 100, barrier at 90, one million paths
args = dict(S0=100, K=100, B=90, T=1.0, r=0.03, v0=0.04, kappa=2.0, theta=0.04, xi=0.3, rho=-0.7)
print(down_and_out_call(**args, n_steps=25, n_paths=1_000_000, bridge=False))   # 8.163, se 0.012
print(down_and_out_call(**args, n_steps=25, n_paths=1_000_000, bridge=True))    # 7.454, se 0.012
print(down_and_out_call(**args, n_steps=250, n_paths=1_000_000, bridge=True))   # 7.451, se 0.012
```

## What the Square Root Cost

We have seen what Heston offers: a variance that mean reverts and stays positive, and a vanilla price that is fast to compute. Beyond vanillas, it prices a forward starting option with a forward smile that keeps its skew regardless of the start date, which local vol could not do. And for path dependent products like the barrier option, where no semi closed form exists, it can still be priced by PDE or Monte Carlo. With those strengths in view, we turn to the model's drawbacks. I would like to examine it from two angles. The first is whether the model behaves the way the market does, both in how vol moves over time and in the shape of the surface it produces. The second is calibration. So far we have focused on pricing, but before any price can be produced, the model has to be fitted to the market, and it is worth checking whether the model is as easy to fit as it is to price.

### Departures from Market Behavior

**Vol of vol.** The vol of vol $\xi$ is a critical parameter that stochastic vol introduces, since it is what makes the vol itself random, so we examine it first. The simplest place to look is the VIX and its own vol index, the VVIX, which Cboe computes from VIX options in the same way the VIX is computed from SPX options. It measures how much the market expects the VIX to move over the next 30 days, as a percentage of its own level. Both histories can be downloaded free from Cboe. Here are five years of daily closes:

{{< vix_vvix >}}

The two move almost in sync, and every stress episode, such as August 2024 and April 2025, lifts both together. Plotting one against the other shows the pattern clearly: on calm days with the VIX around 13 the VVIX sits near 80, and when the VIX is at 30 or above it is typically between 110 and 150, a correlation of 0.61 over the period. When vol is high, the market expects the VIX to move by a larger percentage, not just by more points in absolute terms.

Now compare this to what Heston says. The VIX is quoted in vol while Heston models the variance, so we apply Ito's lemma to $\sigma_t = \sqrt{v_t}$ to get the dynamics of the vol,

$$d\sigma_t = \left(\frac{\kappa(\theta - v_t)}{2\sigma_t} - \frac{\xi^2}{8\sigma_t}\right)dt + \frac{\xi}{2}\,dW_t^v.$$

The shock to the vol is $\tfrac{\xi}{2}\,dW_t^v$, a fixed number that does not depend on the vol level. For example, with $\xi = 0.3$, the diffusion term of $d\sigma_t$ is $0.15\,dW_t^v$. Taking $dt$ as one year, where $dW_t^v$ has a standard deviation of 1, the shock has a standard deviation of 0.15, which is 15 vol points: a vol at 20% would typically be pushed to about 5% or 35% by the shock alone over a year, and a vol at 40% to about 25% or 55%, the same 15 points either way. If the current vol is 20%, Heston says the vol of vol should be $15/20 = 75\%$, and when the current vol is 40%, it should be only $15/40 = 37.5\%$. So in Heston the vol of vol falls as the vol rises, the opposite of what the market shows. Hull-White, the first stochastic vol model we discussed, assumes the variance follows a geometric Brownian motion, which makes the vol a geometric Brownian motion too, so its percentage vol of vol stays constant at any level. That is closer to the data, but still not right, since the market's percentage vol of vol rises with the vol level. Going from Heston to Hull-White, the diffusion coefficient of the variance changes from $\sqrt{v_t}$ to $v_t$, and the percentage vol of vol changes from falling to flat. This motivates raising the power of $v_t$ further, which leads to the 3/2 model, with variance diffusion $\xi v_t^{3/2}\,dW_t^v$.

Knowing that Heston's vol of vol dynamics deviate from the market's should caution us against applying the model to products whose value depends on the vol of vol, such as VIX options. A model whose vol of vol falls as vol rises will tend to underprice VIX calls in exactly the stressed scenarios their buyers pay for. In practice, one common approach avoids modeling how vol evolves altogether: each VIX option is priced with the Black model on the VIX future of the same expiry, using a smile fitted to that expiry's quotes, in the same way equity options are marked off their own implied vol surface. That VIX options exist at all told us vol is random, which led us to look into stochastic vol models. Heston models the vol as a random process, but the way it makes vol random fails to represent the market.

**Skew across maturities.** Another key parameter that stochastic vol introduces is the correlation $\rho$ between the spot and the vol, which produces the skew. We therefore check whether the skew Heston produces matches the market's across maturities.

Before any comparison, we first fit Heston to the market. We take the SPX implied vols on September 15, 2026 at five strikes from 95% to 105% moneyness, across eleven expiries from one week to just over four years, and fit all five parameters to all of them at once. The fit gives $v_0 = 0.023$, $\kappa = 4.56$, $\theta = 0.056$, $\xi = 1.65$ and $\rho = -0.70$. The fit is not perfect, since we are using five parameters to fit 55 quotes, but it is close: across the 55 quotes, the root mean square error is 0.47 vol points, and the largest single miss is 1.4 points. How to calibrate Heston will be explained in the next article. We then compare the ATM skew, the slope of implied vol between 97.5% and 102.5% moneyness, as defined below,

$$\text{ATM skew} = \frac{\sigma(102.5\%) - \sigma(97.5\%)}{\ln 1.025 - \ln 0.975},$$

and plot it against maturity, with maturity on a log scale.

{{< spx_skew_term >}}

From the graph, the market skew flattens from $-1.57$ at one week to $-0.15$ at four years. The fitted Heston curve has a different shape: it flattens at the short end, and at the long end it moves toward zero faster than the market, leaving the long dated smiles too flat. The market's short dated smiles have both more skew and more curvature than the longer dated ones. To match both shapes, the calibration pushes $\xi$ up to 1.65, since, as we saw in [What Sets the Forward Smile](#what-sets-the-forward-smile), the vol of vol affects both: it bends the smile, and it steepens it through the product $\rho\xi$. The downside is that a vol of vol large enough for the one week smile is too large for the rest: it overshoots the skew between one and six months, and the curve is still flatter than the market beyond two years.

Heston's skew levels off at the short end, a known limitation of the model (Gatheral, 2006). At the money and for short maturities, its skew is approximately

$$\text{ATM skew} \approx \frac{\rho\xi}{4\sigma_0}, \qquad \sigma_0 = \sqrt{v_0},$$

which does not depend on the maturity. The derivation is below, which readers already familiar with it can skip. A common fix is to add jumps to the spot, as in the Bates model, which we will explain in a future article.

<details style="margin-bottom: 1.5em;">
<summary>Where the short maturity skew comes from</summary>

The skew is the slope of the implied vol in log strike, $d\sigma/d\log K = K\,d\sigma/dK$. The model price of a call is the Black-Scholes formula at that strike's own implied vol, $C(K) = C_{BS}(K, \sigma(K))$, so its derivative in the strike has two terms,

$$\frac{dC}{dK} = \frac{\partial C_{BS}}{\partial K} + \text{Vega}\times\frac{d\sigma}{dK}.$$

Both strike derivatives are probabilities. Raising the strike by a small amount lowers the payoff $(S_T - K)^+$ by that amount on every path that finishes above $K$, so $-dC/dK = P(S_T > K)$, the model's probability of finishing above the strike. Under Black-Scholes with the vol held fixed, the same probability is $N(d_2)$. Solving for the slope,

$$K\,\frac{d\sigma}{dK} = K\,\frac{N(d_2) - P(S_T > K)}{\text{Vega}}.$$

Assume zero rates, so that the at the money strike is $K = S_0$. At the money, the Black-Scholes price is approximately $0.4\,\sigma S_0\sqrt{T}$, so the vega is $0.4\,S_0\sqrt{T}$, and the skew becomes

$$\text{ATM skew} \approx \frac{N(d_2) - P(S_T > S_0)}{0.4\sqrt{T}}.$$

The skew is now the difference between two probabilities, and both sit close to one half. For small $x$, $N(x) \approx \tfrac{1}{2} + 0.4\,x$ by a first order Taylor expansion at $x = 0$. We start with the lognormal probability. For a short maturity the ATM implied vol is close to $\sigma_0$, so $d_2 = -\tfrac{1}{2}\sigma_0\sqrt{T}$ and

$$N(d_2) \approx \frac{1}{2} - 0.4\,\frac{\sigma_0\sqrt{T}}{2}.$$

Next is Heston's probability of finishing above $S_0$. Over a short horizon the drift of the vol is negligible against its shock, so from the dynamics of $\sigma_t$ derived by Ito earlier, $\sigma_t \approx \sigma_0 + \tfrac{\xi}{2}W_t^v$. Substituting into the log return,

$$\log\frac{S_T}{S_0} \approx -\tfrac{1}{2}\sigma_0^2 T + \sigma_0 W_T^S + \frac{\xi}{2}\int_0^T W_t^v\,dW_t^S.$$

Split the vol driver into the part shared with the spot and an independent part, $W^v = \rho W^S + \sqrt{1 - \rho^2}\,W^\perp$. The independent part is as likely to raise the return as to lower it, so it does not move the probability. The shared part gives $\int_0^T W_t^S\,dW_t^S = \tfrac{1}{2}\big((W_T^S)^2 - T\big)$ by Ito's lemma. Writing $W_T^S = \sqrt{T}\,Z$ with $Z$ a standard normal,

$$\log\frac{S_T}{S_0} \approx -\tfrac{1}{2}\sigma_0^2 T + \sigma_0\sqrt{T}\,Z + \frac{\rho\xi}{4}\,T\,(Z^2 - 1).$$

The spot finishes above $S_0$ when this is positive. When $Z$ is away from zero, the $\sqrt{T}$ term is much larger than the two terms of order $T$ and decides the sign alone, so the $Z^2$ does not matter. When $Z$ is near zero, $Z^2$ is negligible against 1. In both cases $(Z^2 - 1)$ can be replaced by $-1$ without changing the sign, and the condition becomes

$$Z > \frac{\sigma_0\sqrt{T}}{2} + \frac{\rho\xi\sqrt{T}}{4\sigma_0}.$$

By the same Taylor approximation, the probability of finishing in the money under Heston is

$$P(S_T > S_0) \approx \frac{1}{2} - 0.4\left(\frac{\sigma_0\sqrt{T}}{2} + \frac{\rho\xi\sqrt{T}}{4\sigma_0}\right).$$

This gives

$$\text{ATM skew} \approx \frac{N(d_2) - P(S_T > S_0)}{0.4\sqrt{T}} \approx \frac{\rho\xi}{4\sigma_0}.$$

This is the slope exactly at the money in the limit $T \to 0$, and at the fitted parameters it is $-1.90$. The graph measures the slope across the wide range from 97.5% to 102.5%, with one week as its shortest maturity, hence the difference from the value of about $-1.5$ shown there.

</details>

Given that Heston's skew deviates from the market's, products whose value depends on the skew, such as digital options, will be mispriced by the model, even at expiries where its vols are close to the market's. A digital call pays 1 if the spot finishes above the strike $K$, so with zero rates its price is the probability $P(S_T > K)$, which the derivation above writes as $N(d_2) - \text{Vega}\times d\sigma/dK$. The first term is the flat vol price and the second is driven by skew. For a European digital, the better approach is to price it directly off the implied vol surface as a tight call spread, since the surface reproduces the vanilla quotes around the strike, which the fitted Heston model only approximates.

### Challenges in Calibration

Calibration is the prerequisite for any pricing. In the forward start and barrier examples, we assumed the parameters were known, but in practice they come from first fitting the model to vanilla option quotes. One issue is immediately apparent: Heston has only five free parameters, so a fit to a whole surface of quotes cannot be perfect. We already saw this in the discussion of skew across maturities, where the fitted parameters strike a compromise between the skew at the short, middle and long ends. As a result, vanillas priced under the calibrated model do not exactly match their market quotes. This matters for Heston's intended use: an exotic is priced under the model and hedged with vanillas, but the model values those vanillas slightly differently from the prices at which the desk actually trades them, so the hedges are not quite consistent with the market. Ideally, the model used for an exotic would reproduce the vanilla quotes as closely as possible, so that the model values its hedges at the same prices the desk trades them.

In addition to the imperfect fit, the calibrated parameters usually fail the Feller condition. To reproduce the steep equity skew and curvature, the optimizer needs a large vol of vol, and our fit to the SPX surface above shows the condition failing by a wide margin: $\xi = 1.65$ gives $\xi^2 = 2.72$ against $2\kappa\theta = 0.51$. As a result, the variance spends real time near zero, which is not how vol markets behave. The parameters that best fit the market are the ones that let the model's variance reach zero.

Another issue is less obvious, although we already caught a glimpse of it in [What Sets the Forward Smile](#what-sets-the-forward-smile), where raising $\xi$ and lowering $\kappa$ produced almost the same panels. Both act on how widely the variance spreads out, $\xi$ pushing it away from $\theta$ and $\kappa$ pulling it back, and option prices mostly see the balance between the two rather than each on its own. As a result, the parameters are only loosely pinned down. To see how loosely, we return to the fit to the 55 option quotes from the skew discussion, where $\kappa$ came out at 4.56 and $\xi$ at 1.65. Fixing $\kappa$ anywhere from 2.5 to 8, more than a threefold range, and refitting the rest moves $\xi$ from 1.33 to 2.09, while the fit error barely changes, from 0.47 to at most 0.55 vol points. Across that whole range, the best fit differs by less than a tenth of a vol point, so the quotes barely distinguish $\kappa = 2.5$ from $\kappa = 8$, and modest changes in the quotes move the parameters noticeably: raising the vols of the three longest expiries by 0.3 points shifts $\kappa$ from 4.56 to 3.94. For a desk, parameters that drift with small moves in the market make the model's risk and hedges unstable.

### Where This Leads

This article mainly focuses on the Heston model's pricing and properties, and we deferred calibration each time it came up. I will cover the calibration in the next article. Although we can study how to obtain the best fit Heston can offer, the fit will still be imperfect, since five parameters cannot reproduce a whole surface. Of Heston's calibration issues, this is the one that most directly affects the price of an exotic. We have already met a model with the opposite strength. Local vol is built from the vanilla surface through Dupire and reproduces the vanilla quotes closely by construction, but its vol dynamics are oversimplified and its forward smile flattens. The natural question is whether one model can combine the two, keeping Heston's stochastic variance while matching the surface as closely as local vol does. That model is local stochastic volatility, which multiplies the Heston vol by a leverage function built from the local vol surface through Dupire. The LSV article will follow after we dive into Heston calibration, since a calibrated Heston model is its starting point.

## Reflection

Looking back over the article, one tool did a lot of the heavy lifting, and it was not the Fourier transform. It was Feynman-Kac, the theorem that turns an expectation into a PDE, and we invoked it twice under the name of the backward Kolmogorov equation. In pricing, the characteristic function $\mathbb{E}[e^{iuX_T}]$ is the expectation of a payoff, and Feynman-Kac turned it into the two variable PDE that the exponential affine guess reduced to two ODEs. In the Feller derivation, the probability that the rate ever reaches zero looked hard to attack directly, so we recast it as the value of a claim paying a dollar at zero, and Feynman-Kac again turned it into an equation we could solve analytically.

What makes the theorem so useful is that it lets a question in probability be answered with the tools of analysis. A characteristic function and a hitting probability look like different objects, but each is an expectation over the path, so each becomes a PDE. I find it striking that a single theorem sits under both the Heston characteristic function and the Feller condition, and working through them gave me a new appreciation of why Feynman-Kac is one of the cornerstones of quantitative finance.

## Appendix: Solving the Riccati Equation

Recall from the main text that $b = \kappa - \rho\xi iu$ and $d = \sqrt{b^2 + \xi^2(u^2 + iu)}$. Writing $\omega = \frac{u^2 + iu}{2}$, the equation for $B$ is

$$B' = \frac{\xi^2}{2}B^2 - bB - \omega.$$

Its coefficients do not depend on $\tau$, so it can be solved by separating variables. Setting the right side to zero, the quadratic formula gives the two roots

$$\beta_\pm = \frac{b \pm \sqrt{b^2 + 2\xi^2 \omega}}{\xi^2} = \frac{b \pm d}{\xi^2},$$

so the equation factors as

$$B' = \frac{\xi^2}{2}(B - \beta_+)(B - \beta_-).$$

To separate the variables, we divide both sides by the product $(B - \beta_+)(B - \beta_-)$, so that everything in $B$ is on the left and $d\tau$ is on the right. One over a product of two factors cannot be integrated directly, but partial fractions split it into a difference of two simple fractions, each of which integrates to a logarithm. Using $\beta_+ - \beta_- = 2d/\xi^2$,

$$\frac{\xi^2}{2d}\left(\frac{1}{B - \beta_+} - \frac{1}{B - \beta_-}\right)dB = \frac{\xi^2}{2}\,d\tau.$$

Integrating gives

$$\log\frac{B - \beta_+}{B - \beta_-} = d\tau + \text{const}.$$

Exponentiating,

$$\frac{B - \beta_+}{B - \beta_-} = M e^{d\tau}.$$

At $\tau = 0$ we have $B = 0$, so $M = \beta_+/\beta_-$. Solving for $B$,

$$B(\tau) = \beta_+\,\frac{1 - e^{d\tau}}{1 - \frac{\beta_+}{\beta_-}\,e^{d\tau}}.$$

This is the exact solution, but it is not the form we want to compute with. It is written in terms of $e^{d\tau}$, and since the real part of $d$ is positive, that exponential grows very quickly with maturity. On a computer the large numbers eventually overflow, so we rewrite the same solution in terms of the decaying exponential $e^{-d\tau}$. Multiplying the numerator and the denominator by $\frac{\beta_-}{\beta_+}\,e^{-d\tau}$, and writing $g = \beta_-/\beta_+ = \frac{b - d}{b + d}$,

$$B(\tau) = \beta_-\,\frac{1 - e^{-d\tau}}{1 - g e^{-d\tau}} = \frac{b - d}{\xi^2}\cdot\frac{1 - e^{-d\tau}}{1 - g e^{-d\tau}}.$$

This is the form stated in the main text and used in the code.

For $A$, rewrite the fraction as $1 - \frac{(1 - g)e^{-d\tau}}{1 - g e^{-d\tau}}$. The second piece integrates to a logarithm, since the numerator is proportional to the derivative of the denominator, and with $\frac{1 - g}{g d} = \frac{2}{b - d}$,

$$\int_0^\tau B(s)\,ds = \frac{b - d}{\xi^2}\,\tau - \frac{2}{\xi^2}\log\frac{1 - g e^{-d\tau}}{1 - g}.$$

Integrating $A' = iur + \kappa\theta B$ from $A(0) = 0$ then gives the $A(\tau)$ stated in the main text.

As a check, let $\xi \to 0$. Then $b \to \kappa$, $d \approx b + \xi^2(u^2 + iu)/(2b)$, so $\beta_- \to -\frac{u^2 + iu}{2\kappa}$ and $g \to 0$, and

$$B(\tau) \to -\frac{u^2 + iu}{2}\cdot\frac{1 - e^{-\kappa\tau}}{\kappa},$$

which is the backbone's $B(T)$ at $\tau = T$.

## References

Andersen, L. (2008). Simple and efficient simulation of the Heston stochastic volatility model. *Journal of Computational Finance*, 11(3), 1-42.

Bates, D. S. (1996). Jumps and stochastic volatility: Exchange rate processes implicit in Deutsche Mark options. *Review of Financial Studies*, 9(1), 69-107.

Carr, P., and Madan, D. (1999). Option valuation using the fast Fourier transform. *Journal of Computational Finance*, 2(4), 61-73.

Cox, J. C., Ingersoll, J. E., and Ross, S. A. (1985). A theory of the term structure of interest rates. *Econometrica*, 53(2), 385-407.

Duffie, D., Filipović, D., and Schachermayer, W. (2003). Affine processes and applications in finance. *Annals of Applied Probability*, 13(3), 984-1053.

Feller, W. (1951). Two singular diffusion problems. *Annals of Mathematics*, 54(1), 173-182.

Gatheral, J. (2006). *The Volatility Surface: A Practitioner's Guide*. Wiley.

Heston, S. L. (1993). A closed-form solution for options with stochastic volatility with applications to bond and currency options. *Review of Financial Studies*, 6(2), 327-343.

Hull, J., and White, A. (1987). The pricing of options on assets with stochastic volatilities. *Journal of Finance*, 42(2), 281-300.

Schöbel, R., and Zhu, J. (1999). Stochastic volatility with an Ornstein-Uhlenbeck process: An extension. *European Finance Review*, 3(1), 23-46.

Scott, L. O. (1987). Option pricing when the variance changes randomly: Theory, estimation, and an application. *Journal of Financial and Quantitative Analysis*, 22(4), 419-438.

Stein, E. M., and Stein, J. C. (1991). Stock price distributions with stochastic volatility: An analytic approach. *Review of Financial Studies*, 4(4), 727-752.

Wiggins, J. B. (1987). Option values under stochastic volatility: Theory and empirical estimates. *Journal of Financial Economics*, 19(2), 351-372.

Data: VIX and VVIX daily closes from Cboe's historical index data; SPX implied volatilities as of September 15, 2026.

[^negrates]: That concern has aged poorly, with negative policy rates common across Europe and Japan through the 2010s.

[^halfrange]: The integrand is built from waves $e^{iux} = \cos(ux) + i\sin(ux)$. Replacing $u$ by $-u$ leaves the cosine unchanged and flips the sign of the sine, so the real part is even in $u$ and the imaginary part is odd.

[^wing]: A wing option pays only if the spot makes a move of many standard deviations, and the chance of that falls off like $e^{-k^2/(2\sigma^2 T)}$, where $k$ is the log distance to the strike. Each extra point of vol makes the rare move much more likely than the last one did.