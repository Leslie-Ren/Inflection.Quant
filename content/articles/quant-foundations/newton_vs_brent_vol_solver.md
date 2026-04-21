---
title: "Solving for Implied Volatility: Newton's Method vs Brent's Method"
date: 2026-04-09
draft: false
math: true
tags: ["Root Finding", "Newton-Raphson","Brent Method", "Secant Method", "Bisection Method", "Inverse Quadratic Interpolation", "Convergence Rate"]
---

## Why This Matters

I started writing about calibrating a full volatility surface and realised it first requires a clear understanding of a simpler problem: solving for implied vol from a single option price. At its core, this is a root-finding problem: given a market price, we need to find the volatility that makes the model match that price. Once framed this way, the question becomes how to solve this nonlinear problem efficiently and reliably.

I learned Newton's method in class, saw it in textbooks, and implemented it in the traders' tools at work. Most of the time it just works. I knew convergence wasn't guaranteed, but that stayed abstract until I actually hit an edge case. 

So the natural question becomes: why does Newton's method break, and what do you use when it does? That's where Brent's method comes in. It's robust, guaranteed to converge, and commonly used in production. This article is my attempt to build an intuition for both: how they converge, why they fail, and how to combine them into something you can actually trust in production.

---

## 1. Newton's Method (Newton–Raphson)

### The Idea

Newton's method is a gradient-based iterative scheme. The basic idea is to linearize the function around a current guess and step to the root of that linear approximation. For implied volatility, we define:

$$f(\sigma) = C_\text{BS}(\sigma) - C_\text{mkt}$$

Newton's update formula:

$$\sigma_{n+1} = \sigma_n - \frac{f(\sigma_n)}{f'(\sigma_n)}$$

comes directly from a first-order Taylor expansion:

$$f(\sigma) \approx f(\sigma_n) + f'(\sigma_n)(\sigma - \sigma_n)$$

Setting this approximation to zero and solving for $\sigma$ gives the iteration formula above. Intuitively, it moves the current guess toward where the tangent line crosses zero.

---

### Quadratic Convergence

Newton's method converges **quadratically** near the root. Let $\hat{\sigma}$ be the true root of $f(\hat{\sigma})=0$ and define the error $e_n = \sigma_n - \hat{\sigma}$. Using a Taylor expansion around the root:

$$f(\sigma_n) = f(\hat{\sigma} + e_n) = f'(\hat{\sigma})e_n + \frac{1}{2} f''(\xi_n)e_n^2$$

for some $\xi_n$ between $\sigma_n$ and $\hat{\sigma}$. The next step gives an updated error:

$$e_{n+1} = \sigma_{n+1} - \hat{\sigma} \approx -\frac{f''(\xi_n)}{2f'(\hat{\sigma})}e_n^2$$

Thus the error satisfies $e_{n+1} \propto e_n^2$, proving quadratic convergence: the number of correct digits roughly doubles each iteration once near the root.

---

### When Newton's Method Fails

Newton's method is not guaranteed to converge. In practice, I've seen it converge in a few iterations for a liquid ATM option but oscillate or diverge entirely on a 5-delta wing. Why?

This method is derived from a first-order Taylor expansion of the Black–Scholes price around the current volatility estimate:

$$C_{\text{BS}}(\sigma + \Delta\sigma) \approx C_{\text{BS}}(\sigma) + \mathcal{V}\,\Delta\sigma$$
Here vega $\mathcal{V} = \frac{\partial C_{\text{BS}}}{\partial \sigma}$ corresponds to $f'$ in the earlier notation.

This gives the Newton step:

$$\Delta\sigma = \frac{C_{\text{market}} - C_{\text{BS}}(\sigma)}{\mathcal{V}}$$

This update relies on two implicit assumptions: vega is well-defined and nonzero, and the price function is well-approximated by a linear function in volatility over the step size. When either breaks down, Newton's method becomes unstable or inaccurate.

---

#### 1. Vega Near Zero — Step Size Instability

Vega $\mathcal{V}$ appears in the denominator of the Newton step. When $\mathcal{V} \approx 0$, even a small pricing error produces a large volatility update, causing the iterate to move far from the solution. Vega becomes small when the option is **deep in-the-money or deep out-of-the-money**, or when **time to maturity is very short**. In these regimes, the option price is close to intrinsic value and has limited sensitivity to volatility.

In this regime, the update becomes ill-conditioned:

$$\Delta\sigma = \frac{C_{\text{market}} - C_{\text{BS}}(\sigma)}{\mathcal{V}} \quad \text{becomes unstable as } \mathcal{V} \to 0$$

This is not because the root does not exist, but because the inverse problem becomes poorly conditioned: small errors in price translate into large errors in volatility.

---

#### 2. Poor Initial Guess — Curvature (Vomma) Breakdown

Newton's method is a local approximation scheme. It relies on truncating the Taylor expansion:

$$C_{\text{BS}}(\sigma + \Delta\sigma) \approx C_{\text{BS}}(\sigma) + \mathcal{V}\Delta\sigma + \frac{1}{2}\frac{\partial^2 C_{\text{BS}}}{\partial \sigma^2}(\Delta\sigma)^2 + \cdots$$

The second-order term involves vomma, $\frac{\partial^2 C_{\text{BS}}}{\partial \sigma^2}$, which measures the curvature of the option price with respect to volatility. When the initial guess is far from the true implied volatility, the Newton step $\Delta\sigma$ becomes large. In this regime, the omitted quadratic term is no longer negligible and the linear approximation breaks down.

As a result, the tangent line no longer accurately predicts where the price curve intersects the market price level. The Newton update can overshoot the root, placing the next iterate on the opposite side of the solution. Repeated overshoots can lead to oscillation, and in extreme cases, failure to converge.

The interactive chart below lets you drag both the true implied vol and the initial guess to observe how the method behaves. It uses a European call with $S = 100$, $K = 130$, $r = 5\%$, $T = 0.5$ years, and no dividends.

{{< newton_method_initial_guess>}}

---

### Improving Newton's Method in Practice

- **Smart initial guess:** Use the Brenner-Subrahmanyam approximation for near-the-money options ([discussed in this article]({{< relref "understanding_brownian_motion.md" >}})) and the Corrado-Miller approximation for options away from the money. A start close to the true implied vol suppresses curvature error on the first Newton step and reduces iteration count significantly.
- **Bounded volatility updates:** Clamp the iterate to an admissible interval (e.g. $[10^{-4}, 5]$) after every update, and cap the step size to prevent large single-step overshoots.
- **Robust fallback:** When vega falls below a threshold, switch to a bracketing method such as Brent's method to guarantee convergence. This hybrid approach preserves Newton's fast convergence in well-behaved regimes while remaining robust at the boundaries.

---

## 2. Brent's Method

### The Idea

Brent's method is best understood not as a single algorithm, but as an **adaptive system that combines three root-finding methods** — bisection, the secant method, and inverse quadratic interpolation. At each iteration, it selects the most aggressive step that satisfies its safety conditions, and falls back to a more conservative method when those conditions are not met. The result is a solver that is simultaneously guaranteed to converge and capable of superlinear acceleration whenever the function is well-behaved.

The foundation is a **bracket** $[\sigma_a, \sigma_b]$ satisfying:

$$f(\sigma_a) \cdot f(\sigma_b) < 0$$

meaning the pricing error $f(\sigma) = C_\text{BS}(\sigma) - C_\text{mkt}$ changes sign across the interval, guaranteeing a root lies within. This bracket is maintained throughout every iteration — it is the safety guarantee that Newton's method lacks.

---

### The Three Building Blocks

Brent's method can be viewed as a safeguarded interpolation scheme: bisection guarantees global convergence, while interpolation provides local acceleration whenever the function behaves well.

---

#### Bisection — Slow but Unconditionally Converges

At each iteration, bisection evaluates $f$ at the midpoint of the current bracket:

$$\sigma_\text{mid} = \frac{\sigma_a + \sigma_b}{2}$$

and replaces whichever endpoint shares the same sign as $f(\sigma_\text{mid})$, halving the interval. The bracket shrinks by exactly half each step, regardless of the shape of $f$.

**Convergence intuition:** Bisection uses only the *sign* of the function — not its magnitude or slope. It discards almost all quantitative information at each step. This is why it is slow: halving the interval each time gives linear convergence at rate $\frac{1}{2}$, requiring $\lceil \log_2(W/\varepsilon) \rceil$ iterations to reduce a bracket of width $W$ to tolerance $\varepsilon$. For a bracket $[10^{-4}, 5]$ and tolerance $10^{-8}$, that is 29 iterations with no possibility of acceleration.

---

#### Secant Method — Fast but No Guarantee

The secant method fits a straight line through the two most recent iterates $(\sigma_a, f_a)$ and $(\sigma_b, f_b)$, and steps to where that line crosses zero:

$$\sigma_\text{new} = \sigma_b - f_b  \frac{\sigma_b - \sigma_a}{f_b - f_a}$$

Unlike bisection, it uses the actual values of $f$ at both endpoints — not just their signs — to estimate where the root is.

Although this looks similar to Newton's method — replacing the analytic derivative with a finite difference — the secant method is fundamentally a different approach. Newton's method always evaluates the derivative at the *current* point, anchoring each step to local curvature. The secant method instead draws a chord through the **two most recent iterates**, making it entirely derivative-free and driven purely by recent function history. 

**Convergence intuition:** By using function values, the secant method can take a much more informed step than bisection. Rather than depending on a single previous error term, the next error is influenced by both of the two most recent errors. In fact, a local asymptotic analysis shows that near the root,
$$
e_{n+1} \approx Ce_n e_{n-1}
$$
for some constant $C$ determined by the local behavior of $f$ . This coupling between successive errors leads to superlinear convergence. In particular, the asymptotic convergence rate is approximately $\varphi \approx 1.618$, the golden ratio. Intuitively, each iteration amplifies the effect of the two previous error reductions, producing faster decay than linear methods but slower than quadratic methods like Newton's method. The cost is that without a bracket constraint, the secant step can overshoot the root if $f$ is highly curved between the two points, and there is no convergence guarantee in general.

---

#### Inverse Quadratic Interpolation — Fastest but Most Fragile

Inverse quadratic interpolation (IQI) fits a **quadratic polynomial through the three most recent iterates** $(\sigma_a, f_a)$, $(\sigma_b, f_b)$, $(\sigma_c, f_c)$. Critically, it fits $\sigma$ as a function of $f$ — not $f$ as a function of $\sigma$ — so that evaluating at $f = 0$ directly yields the next iterate. Using Lagrange interpolation:

$$\sigma_\text{new} = \sigma_a \frac{f_b f_c}{(f_a - f_b)(f_a - f_c)} + \sigma_b \frac{f_a f_c}{(f_b - f_a)(f_b - f_c)} + \sigma_c \frac{f_a f_b}{(f_c - f_a)(f_c - f_b)}$$

**Convergence intuition:** By incorporating a third point and fitting a quadratic, IQI captures the local curvature of $f$ — the information that the secant method ignores. Near the root, a local asymptotic analysis shows that the error satisfies a higher-order nonlinear recurrence involving the three most recent iterates. This leads to an asymptotic convergence order of approximately $q \approx 1.839$. Intuitively, each additional interpolation point increases the amount of local structure captured by the model, allowing the method to reduce the error more aggressively than two-point methods. The trade-off is that fitting a quadratic through three points can be numerically unstable when the points are poorly spaced. If the corresponding function values are close together or do not adequately span the root, the interpolation becomes ill-conditioned: small differences in the data lead to large changes in the fitted curve. Because IQI effectively extrapolates to $f = 0$, this instability can produce a step that lies far outside the current bracket. This is why IQI is not used in isolation and instead requires a bracketing safeguard.

---

### How Brent Combines Them

At each iteration, Brent's method first checks whether the bracket width has reduced below a tolerance $\delta$ — if so, the method terminates. Otherwise, it proposes an IQI step (or secant step if only two distinct function values are available) and accepts it only if two safety conditions hold:

1. **The step lands inside the current bracket** $[\sigma_a, \sigma_b]$
2. **The step represents sufficient progress toward the root** — A common implementation heuristic is to reject the step if it is larger than roughly half the previous bracket width, treating such moves as insufficiently controlled.

If either condition fails, the method falls back to bisection. This gives Brent's method its core character: it runs as fast as IQI or the secant method when the interpolation is well-behaved, but is guaranteed to make at least the progress of bisection at every step.

$$\sigma_\text{new} = \begin{cases} \sigma_\text{IQI or secant} & \text{if both safety conditions are satisfied} \\\\ \sigma_\text{mid} & \text{otherwise (bisection fallback)} \end{cases}$$

---

#### Evaluating and Updating the Bracket

Once $\sigma_\text{new}$ is determined — regardless of which method proposed it — the following three steps are always executed:

**Step 1: Evaluate the pricing error at the new point**

$$f(\sigma_\text{new}) = C_\text{BS}(\sigma_\text{new}) - C_\text{mkt}$$

The sign of $f(\sigma_\text{new})$ determines which half of the bracket contains the root.

**Step 2: Narrow the bracket using the sign of $f(\sigma_\text{new})$**

Replace whichever endpoint shares the same sign as $f(\sigma_\text{new})$. The bracket always shrinks after every iteration — this is the convergence guarantee that no interpolation step can violate.

**Step 3: Promote $\sigma_\text{new}$ as the best current estimate**

$\sigma_\text{new}$ becomes the new endpoint — Brent always keeps the best estimate as one of the bracket endpoints. The previous endpoint becomes the third point available for the next IQI step.

---

### When Brent's Method Struggles

Brent's method is significantly more robust than Newton's, but it is not without limitations.

#### 1. Slower Convergence Near the Root

Once Newton's method enters its quadratic convergence regime it is faster than Brent. For applications requiring very high precision — such as calibrating a vol surface across many strikes simultaneously — the difference in iteration count can matter. Brent's superlinear rate means it typically requires more iterations than Newton to achieve the same terminal accuracy, assuming Newton does not encounter the failure modes described above.

#### 2. Requires a Valid Initial Bracket

Brent's method requires two initial values $\sigma_a$ and $\sigma_b$ such that $f(\sigma_a) f(\sigma_b) < 0$. In practice, for implied volatility, constructing such a bracket is rarely difficult because economically reasonable volatility ranges are well known (e.g. $[10^{-4}, 5]$). As a result, bracket initialization is usually a one-time evaluation rather than a significant computational overhead. However, it remains a structural requirement of the method, in contrast to Newton's method, which only requires a single starting point.

---

### Brent vs Newton: When to Use Which

Neither method dominates in all regimes. The choice depends on the structure of the problem.

| Scenario | Preferred Method |
|---|---|
| Near-the-money, liquid option | Newton — fast quadratic convergence |
| Deep ITM / OTM, short expiry | Brent — robust when vega is near zero |
| Poor or unknown initial guess | Brent — bracketing guarantees convergence |
| High-precision vol surface calibration | Newton with smart initial guess, provided vega is well-conditioned |
| Production solver requiring robustness | Hybrid: Newton with Brent fallback |

In practice, a well-engineered implied vol solver uses Newton's method as the primary engine and falls back to Brent when vega is too small or the Newton iterate leaves the admissible domain. This hybrid approach inherits the speed of Newton in normal regimes and the reliability of Brent at the boundaries.

A natural question is why not use bisection alone as the fallback rather than Brent. The answer is that bisection is reliable but slow. Brent already contains bisection as its internal worst-case fallback, but accelerates with inverse quadratic interpolation whenever safe to do so. In practice, using Brent as a fallback provides a more structured and generally faster alternative to pure bisection, while maintaining the same global convergence guarantees.


*For a full derivation of the convergence rate for secant method or IQI, feel free to  [contact me](mailto:Infelction.quant@gmail.com).*
