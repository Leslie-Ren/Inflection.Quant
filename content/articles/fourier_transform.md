---
title: "Fourier Transform: The Leap into Frequency"
date: 2026-06-16
draft: false
math: true
tags: ["fourier-transform", "pde", "black-scholes", "heat-equation", "characteristic-function", "heston"]
---

## Why This Matters

The previous article on [finite difference methods]({{< relref "fdm.md" >}}) solved the heat equation by brute force: lay down a grid, step it through time, let the solution emerge step by step. It works, and for many pricing problems it is the practical choice. But can we instead solve the equation analytically?

For a class of PDEs, the heat equation among them, we can. The idea is to stop viewing a function as a shape over space and instead see it as a combination of frequencies. That change of view is the Fourier transform. In this article, we explore it on the heat equation, where the diffusion structure shows through without the variable coefficients of Black-Scholes to clutter it, and then turn to its application in option pricing, where it works even when the distribution of prices has no closed form.

---

## The Heat Equation Again

Take the same setup as before: a metal rod of length $L$, insulated along its sides so heat only flows along its length. Its temperature $u(x, t)$ depends on position $x \in [0, L]$ and time $t$, and evolves under

$$\frac{\partial u}{\partial t} = \kappa \frac{\partial^2 u}{\partial x^2},$$

where $\kappa$ is the thermal diffusivity. We know the temperature profile at $t = 0$, call it $u(x, 0) = f(x)$, and we want the profile at any later time.

The numerical method took $f(x)$ on a grid and advanced it through time. To solve it analytically we have nothing yet, so the honest first move is just to look for solutions, any solutions, and see what the equation will allow.

---

## Looking for Solutions

What makes this equation hard to solve directly is that it ties space and time together: the rate of change in time at a point depends on the curvature in space around that point, so the whole profile evolves as one coupled object. Rather than attack that coupling head-on, watch what the rod actually does as it cools, and look for solutions that behave the same way.

Heat tends to fade in place. A hump of warmth in the middle of the rod stays a hump in the middle; it does not drift to one end or sprout a second peak. The shape of the profile holds while its height drains away. If we look for solutions that do exactly this, hold a fixed shape and simply shrink, we are looking for a spatial profile $X(x)$ multiplied by a time-dependent amplitude $T(t)$ that scales it up or down:

$$u(x, t) = X(x)\, T(t).$$

The figure below shows one such shape cooling: it loses height while holding its form, every point cooling in the same proportion. This is the behavior we are looking for.

{{< heat_decay_single >}}

With the form in hand, substitute it into the heat equation. The left side differentiates only $T$, the right side only $X$:

$$X(x)\, T'(t) = \kappa\, X''(x)\, T(t).$$

Move everything in $T$ to one side and everything in $X$ to the other:

$$\frac{T'(t)}{\kappa\, T(t)} = \frac{X''(x)}{X(x)}.$$

A function of $t$ alone can equal a function of $x$ alone, for all $x$ and $t$, only if both are the same constant. Call it $-\lambda$ (the sign is a convention chosen because we expect decay):

$$\frac{X''(x)}{X(x)} = -\lambda, \qquad \frac{T'(t)}{\kappa\, T(t)} = -\lambda.$$

Two ordinary differential equations, each in one variable, which we can now solve independently.

---

## The Spatial Part and the Role of Boundary Conditions

The spatial equation is $X''(x) = -\lambda X(x)$. For $\lambda > 0$ its solutions are oscillations,

$$X(x) = A \cos(\omega x) + B \sin(\omega x), \qquad \omega = \sqrt{\lambda}.$$

Both sine and cosine solve the equation. Which combination we keep is not a mathematical preference; it is forced by what we hold fixed at the ends of the rod.

Consider two different boundary conditions.

*Ends held at zero temperature.* If the rod is clamped to ice baths at both ends, then $u(0, t) = u(L, t) = 0$ for all time, so $X(0) = X(L) = 0$. The condition $X(0) = 0$ kills the cosine term, since $\cos(0) = 1$ and the cosine cannot vanish at the origin without disappearing entirely; only the sine, with $\sin(0) = 0$, survives. The condition $X(L) = 0$ then forces $\sin(\omega L) = 0$, so $\omega L$ is a multiple of $\pi$. The allowed shapes are

$$X_n(x) = \sin\left(\frac{n \pi x}{L}\right), \qquad n = 1, 2, 3, \dots$$

*Ends insulated.* If instead the ends are insulated so no heat flows across them, the condition is on the slope, not the value: $u_x(0, t) = u_x(L, t) = 0$. Now it is the sine that dies, because its derivative $\cos$ fails to vanish at the origin, and the cosine survives:

$$X_n(x) = \cos\left(\frac{n \pi x}{L}\right), \qquad n = 0, 1, 2, \dots$$

Same equation, same family of oscillations, opposite selection. Fixed-temperature ends give a sine basis; insulated ends give a cosine basis. A mix of conditions (one end fixed, one insulated) gives a basis of shifted oscillations that are part sine and part cosine. The general lesson: the boundary conditions choose which oscillations are admissible. We will carry the fixed-temperature (sine) case forward for simplicity, but everything that follows has a cosine twin.

In every case the admissible frequencies are discrete: a countable list $\omega_n = n\pi/L$, spaced $\pi/L$ apart.

---

## The Time Part

The time equation $T'(t) = -\lambda \kappa\, T(t)$ says the rate of change is proportional to the current value, which is solved by an exponential. The spatial side already fixed $\lambda = \omega_n^2$, so

$$T_n(t) = T_n(0)\, e^{-\kappa \omega_n^2 t}.$$

Each spatial shape decays exponentially, and the decay rate is $\kappa \omega_n^2$. Higher frequency means faster decay, and the dependence is on the *square* of the frequency. This single expression states the diffusion behavior outright: sharp features, which are built from high-frequency shapes, decay fast; smooth features, built from low frequencies, persist.

A single separated solution is therefore the chosen spatial shape times its decay, with the amplitude set aside for now,

$$u_n(x, t) = \sin\left(\frac{n \pi x}{L}\right) e^{-\kappa \omega_n^2 t}.$$

In the insulated-end problem the sine here is a cosine and nothing else changes: the time factor depends only on $\omega_n$, which is the same list of frequencies either way. Whatever the boundary conditions select for the spatial shape, that shape decays at rate $\kappa\omega_n^2$. This solution satisfies the heat equation and the boundary conditions. What it does not yet satisfy is the initial condition: at $t = 0$ it is a single sine, and our actual starting profile $f(x)$ is some arbitrary shape.

---

## Fourier's Bold Assumption

A single sine cannot match an arbitrary profile, so Fourier proposed adding many: that with enough sines and the right amplitudes, the sum can match essentially any starting profile we would meet in practice. For the ice-bath rod we are carrying forward, that means

$$f(x) = \sum_{n=1}^{\infty} b_n \sin\left(\frac{n \pi x}{L}\right),$$

and the insulated case is the same statement with cosines, $f(x) = \tfrac{1}{2}a_0 + \sum_{n=1}^{\infty} a_n \cos(n\pi x/L)$. This expansion of a function into its sine and cosine components is the **Fourier series**.

At the time this was a genuinely audacious claim: that a function with corners, jumps, and flat stretches, the kind no single formula describes, could be reconstructed entirely from smooth oscillations.

Part of the claim is straightforward. Adding solutions is allowed: if $u_1$ and $u_2$ each solve the heat equation, then because each derivative acts on the two pieces separately,

$$\frac{\partial (u_1 + u_2)}{\partial t} = \frac{\partial u_1}{\partial t} + \frac{\partial u_2}{\partial t} = \kappa\frac{\partial^2 u_1}{\partial x^2} + \kappa\frac{\partial^2 u_2}{\partial x^2} = \kappa\frac{\partial^2 (u_1 + u_2)}{\partial x^2},$$

so the sum solves it too, as does any amplitude-weighted sum of sines. The hard part is the reach of the claim: that the sum can be made to equal an arbitrary $f(x)$ in the first place. Fourier asserted this but did not prove it; the first rigorous proof came decades later from Dirichlet, who gave sufficient conditions on the function for the series to converge.

Those conditions are that the function has to be bounded, have at most finitely many jumps and corners, and a finite number of maxima and minima over the interval. A profile like $1/x$ near the origin, which runs off to infinity, cannot be represented. Physical temperature profiles satisfy all of this, so the conditions are no real restriction here.

Granting the assumption, the problem reduces to one question: given $f(x)$, what are the coefficients $b_n$?

---

## Extracting a Coefficient Is a Projection

Fourier's contribution was the assumption that the decomposition exists; that the sines are orthogonal was already known before him. Orthogonality is what turns the assumption into a recipe, making each coefficient a clean projection.

We will borrow the analogy from vectors. Two vectors are independent when their dot product is zero. For $\mathbf{u} = (u_1, u_2, \dots, u_k)$ and $\mathbf{v} = (v_1, v_2, \dots, v_k)$ the dot product multiplies them component by component and adds, $\mathbf{u}\cdot\mathbf{v} = \sum_i u_i v_i$.

A function is like a vector with a component at every point $x$ rather than finitely many, so the sum over components becomes an integral:

$$\langle g, h \rangle = \int_0^L g(x)\, h(x)\, dx.$$

Apply it to two different sines and the integral is zero:

$$\int_0^L \sin\left(\frac{n \pi x}{L}\right) \sin\left(\frac{m \pi x}{L}\right) dx = \begin{cases} 0 & n \neq m \\ L/2 & n = m. \end{cases}$$

The sines are orthogonal in the same sense as perpendicular vectors, and a sine dotted with itself gives a nonzero length, here $L/2$.

Beyond the integral, the orthogonality of two sines is visible in a picture.

{{< orthogonality >}}

Take $\sin(\pi x/L)$ and $\sin(2\pi x/L)$. Their product is positive where both have the same sign and negative where they differ, and over the rod those regions cancel exactly: for every stretch where the product is positive there is a matching stretch where it is equally negative. The integral adds up the product, so the cancellation drives it to zero.

This orthogonality is what lets us recover a single coefficient. To get $b_m$, dot both sides of the decomposition with the sine we want: multiply by $\sin(m\pi x/L)$ and integrate over the rod.

$$\int_0^L f(x) \sin\left(\frac{m \pi x}{L}\right) dx = \sum_{n=1}^{\infty} b_n \int_0^L \sin\left(\frac{n \pi x}{L}\right)\sin\left(\frac{m \pi x}{L}\right) dx.$$

Orthogonality collapses the entire sum on the right to its single $n = m$ term, which carries the factor $L/2$. Solving,

$$b_m = \frac{2}{L} \int_0^L f(x) \sin\left(\frac{m \pi x}{L}\right) dx.$$

The coefficient is the projection of $f$ onto the $m$-th sine, divided by that sine's self-overlap $L/2$. The cosine coefficients come out the same way, projecting onto cosines instead.

---

## Solving the Heat Equation: Bounded Case

We now have all three parts: the sine shapes, their decay rates, and the coefficients that match the initial profile. Combining them, the solution is the decomposition with each term carrying its time factor:

$$u(x, t) = \sum_{n=1}^{\infty} b_n \sin\left(\frac{n \pi x}{L}\right) e^{-\kappa (n\pi/L)^2 t}.$$

This is the ice-bath rod. The insulated rod has the same form with cosines, and a rod with mixed conditions carries both sines and cosines together; the boundary conditions decide which shapes appear, not whether the two kinds can share a solution.

Read it left to right and the physics is in plain view. At $t = 0$ the exponentials are all $1$ and we recover $f(x)$. As time advances, every sine shrinks, and the high-$n$ sines shrink fastest because the rate goes as $n^2$. The profile loses its sharp features first and relaxes toward smoothness. The grid used by the finite difference method would approach this decay one time step at a time; the analytical form gives the entire time evolution in closed form, with each term's decay rate $\kappa(n\pi/L)^2$ appearing explicitly in its exponent.

A concrete profile makes this tangible. Take a rod sitting at a uniform temperature $u_0$, then clamped to ice at both ends at $t = 0$, so $f(x) = u_0$. The coefficients are a single elementary integral,

$$b_n = \frac{2}{L}\int_0^L u_0 \sin\frac{n\pi x}{L}\,dx = \frac{2u_0}{n\pi}\left(1 - \cos n\pi\right),$$

which is $4u_0/(n\pi)$ for odd $n$ and zero for even $n$. The solution is then

$$u(x,t) = \frac{4u_0}{\pi}\sum_{n \text{ odd}} \frac{1}{n}\sin\frac{n\pi x}{L}\,e^{-\kappa(n\pi/L)^2 t}.$$

{{< uniform_decay >}}

The amplitudes fall off as $1/n$, so the $n=1$ sine is the largest from the start, and the $e^{-\kappa(n\pi/L)^2 t}$ factors then shrink the higher-$n$ sines far faster than it. Within a short time only the $n=1$ term carries any weight, and the rod relaxes toward a single smooth arch bowing down to the cold ends.

---

## What Happens When the Rod Grows Infinitely Long?

Everything so far rests on the rod being finite. The ends at $0$ and $L$ are what forced the frequencies into a discrete list $\omega_n = n\pi/L$. But what happens if the rod grows infinitely long, with no boundary conditions at all? Lengthen it and the spacing $\pi/L$ shrinks; the admissible sines crowd closer together. In the limit $L \to \infty$ the spacing goes to zero and the allowed frequencies fill in to a continuum: every $\omega$, not just the integer multiples of $\pi/L$, is now in play.

Two things change in that limit. The first is mechanical. With the frequencies packed infinitely densely, the sum over the discrete list $\sum_n$ becomes an integral over the continuous variable $\omega$, and the sequence of coefficients $b_n$, one number per allowed frequency, becomes a function $b(\omega)$ defined for every real $\omega$. A decomposition into countably many sines becomes a decomposition into a continuum of them:

$$f(x) = \int_0^{\infty} b(\omega)\sin(\omega x)\,d\omega.$$

The second change is subtler. On the rod $[0,1]$, two near-equal frequencies are nearly the same direction: $\int_0^1 \sin(\pi x)\sin(1.0001\pi x)\,dx = 0.49997$, almost a sine's full overlap with itself. Over the whole line, though, any two distinct frequencies are exactly orthogonal, the integral understood as the limit of the overlap over $[-L, L]$ as $L \to \infty$:

$$\int_{-\infty}^{\infty} \sin(\omega_1 x)\sin(\omega_2 x)\,dx = 0 \quad \text{whenever } \omega_1 \neq \omega_2.$$

The intuition is that on the infinite interval, as long as $\omega_1 \neq \omega_2$, the faster sine gradually slides relative to the slower one, spending exactly as much of the line moving in step with it as moving against it. The matched and opposed stretches are equal, so the overlap cancels to zero.

Extracting the amount of a given frequency is still the same overlap integral as on the rod, now run over the whole line for every $\omega$ rather than each integer $n$. This decomposition into a continuum of orthogonal oscillations is the **Fourier transform**. In sine and cosine form it reads

$$f(x) = \int_0^{\infty} \big[\,a(\omega)\cos(\omega x) + b(\omega)\sin(\omega x)\,\big]\,d\omega,$$

with the weights recovered by the same overlap as before, now integrated over the whole line:

$$a(\omega) = \frac{1}{\pi}\int_{-\infty}^{\infty} f(x)\cos(\omega x)\,dx, \qquad b(\omega) = \frac{1}{\pi}\int_{-\infty}^{\infty} f(x)\sin(\omega x)\,dx.$$

---

## Solving the Heat Equation: Unbounded Case

Earlier we solved the heat equation on a rod with fixed ends. Now take a rod with no boundaries at all and a concrete starting profile: a quantity of heat $Q$ deposited at the single point $x = 0$, a hot needle touched to the origin of an infinite rod, with the rest of the rod cold. Concentrating the heat at a point forces one shift in reading: $u(x,t)$ is now the heat density, the amount of heat per unit length, rather than the temperature.[^density] Written as an initial condition, all the heat sits at one point, $u(x, 0) = Q\,\delta(x)$, where $\delta$ is the spike at the origin that integrates to one, so the total heat on the rod is $Q$. We never need a formula for $\delta$ itself, only the rule that defines it: integrated against any function it returns that function's value at the origin,

$$\int_{-\infty}^{\infty} \delta(x)\,g(x)\,dx = g(0).$$

The plan is the same three steps as the bounded case: decompose the initial profile into frequencies, let each frequency decay at its own rate, reassemble.

Decompose first. The coefficients come from the same projection as on the rod, now run over the whole line. Projecting the point source onto each frequency, the delta collapses the integral to the integrand's value at the origin:

$$a(\omega) = \frac{1}{\pi}\int_{-\infty}^{\infty} Q\,\delta(x)\cos(\omega x)\,dx = \frac{Q}{\pi}\cos(0) = \frac{Q}{\pi}, \qquad b(\omega) = \frac{1}{\pi}\int_{-\infty}^{\infty} Q\,\delta(x)\sin(\omega x)\,dx = \frac{Q}{\pi}\sin(0) = 0.$$

The cosine weight is the same constant at every $\omega$, so the spike projects equally onto all frequencies, and the sine weight is zero everywhere. The profile is built from cosines alone. Each frequency decays as $e^{-\kappa\omega^2 t}$, falling off from $\omega = 0$ as the high frequencies die away fastest. Reassemble by integrating these decayed cosines back into a profile in $x$:

$$u(x,t) = \frac{Q}{\pi}\int_0^{\infty} e^{-\kappa\omega^2 t}\cos(\omega x)\,d\omega.$$

This is a Gaussian integral, and it evaluates in closed form to[^gaussint]

$$u(x,t) = \frac{Q}{2\sqrt{\pi\kappa t}}\,e^{-x^2/(4\kappa t)}.$$

This bell curve is exactly a normal distribution, scaled by the total heat $Q$. Reading its exponent $-x^2/(4\kappa t)$ against the normal's $-x^2/(2\sigma^2)$ sets the variance at $\sigma^2 = 2\kappa t$, growing linearly in time. So the whole evolution collapses to one picture: we begin with all the heat at a single point, the zero-variance limit of a normal, and it spreads as a normal whose variance widens as $2\kappa t$. At small $t$ the curve is tall and narrow, almost all the heat still near $x = 0$; as $t$ grows the variance grows with it, spreading the heat wider while the area underneath, the total heat $Q$, stays fixed.

{{< heat_kernel >}}

This spreading Gaussian is the heat kernel, the exact closed form for a point source on the infinite rod. Had the heat been deposited at $x_0$ rather than the origin, the same spreading would produce the identical Gaussian centered there, $u(x,t) = \frac{Q}{\sqrt{4\pi\kappa t}}\,e^{-(x-x_0)^2/(4\kappa t)}$, since a uniform rod diffuses the same way wherever the source sits.

The point source was a special initial condition, but it unlocks the general solution. Because the heat equation is linear, a general starting profile $f(x)$ can be read as a sum of point sources, one at each location $y$ carrying heat $f(y)$, and the response to the whole is the sum of the responses to each. Every source at $y$ spreads into a kernel centered on $y$, so the solution is the initial profile weighted against the kernel at each point,

$$u(x,t) = \int_{-\infty}^{\infty} f(y)\,\frac{1}{\sqrt{4\pi\kappa t}}\,e^{-(x-y)^2/(4\kappa t)}\,dy.$$

This integral is the general solution on the infinite rod: the heat at each point is the initial profile averaged against a Gaussian of variance $2\kappa t$ centered there, the nearby heat counting most. The point source is the case $f = Q\delta$, where the integral collapses back to the single kernel above.

---
## The Fourier Transform in Complex Exponentials

So far we have built the transform from sine and cosine, the way Fourier originally worked, but the modern representation is almost always written with the complex exponential $e^{i\omega x}$ instead. The two are equivalent by Euler's formula,

$$e^{i\omega x} = \cos(\omega x) + i\sin(\omega x),$$

which rearranges to give the sine and cosine back as exponentials,

$$\cos(\omega x) = \frac{e^{i\omega x} + e^{-i\omega x}}{2}, \qquad \sin(\omega x) = \frac{e^{i\omega x} - e^{-i\omega x}}{2i}.$$

Therefore

$$a(\omega)\cos(\omega x) + b(\omega)\sin(\omega x) = c(\omega)\,e^{i\omega x} + c(-\omega)\,e^{-i\omega x},$$

where $c(\omega) = \tfrac{1}{2}\big(a(\omega) - i\,b(\omega)\big)$ and its partner is the complex conjugate, $c(-\omega) = \tfrac{1}{2}\big(a(\omega) + i\,b(\omega)\big) = \overline{c(\omega)}$ for real $f$.

The sine-and-cosine decomposition and the complex-exponential decomposition therefore hold the same information. In exponential form the transform and its inverse read

$$\hat f(\omega) = \int_{-\infty}^{\infty} f(x)\, e^{-i\omega x}\, dx, \qquad f(x) = \frac{1}{2\pi}\int_{-\infty}^{\infty} \hat f(\omega)\, e^{i\omega x}\, d\omega,$$

where $\hat f(\omega) = 2\pi\,c(\omega)$, and the complex weight $\hat f(\omega)$ carries both the amplitude and the phase of frequency $\omega$. A single frequency is one oscillation described by two numbers, how large it is and where its peaks sit along the axis; a sine and a cosine of the same $\omega$ are this one oscillation in two shifted positions. The complex $\hat f(\omega)$ records both: its size and its angle in the plane are the two numbers that pin down the oscillation.

Two things make this form worth adopting. The first is bookkeeping: sines and cosines are tracked together as one complex exponential, instead of as two separate real families. The second matters more but is not obvious at this point: the complex exponential is the language of the characteristic function, $\phi(\omega) = \mathbb{E}[e^{i\omega X}]$, which is the transform of the density evaluated at $-\omega$ under the sign convention above, and which we will draw on in the discussion of option pricing that follows.


---

## Option Pricing

### Black-Scholes

The Black-Scholes equation is the heat equation in disguise. A change of variables turns it into exactly the diffusion we just solved, so the same machinery prices an option. We carry this out in full in the appendix, reducing the equation to the heat equation and recovering the Black-Scholes formula.

This is included to illustrate the technique on a familiar case, not because the transform is needed here: the Black-Scholes density is known in closed form, lognormal, so the pricing integral can be done directly. Nor is this how Fourier methods are used in practice. Their real value appears when no such closed-form density is available, such as in the Heston stochastic volatility model below.

### Heston

A price is the discounted expected payoff, the payoff $g$ integrated against the density $q$ of the log price $y = \ln S_T$,

$$V = e^{-rT}\int_{-\infty}^{\infty} g(y)\, q(y)\, dy.$$

Stochastic volatility models, Heston among them, leave $q$ with no closed form, so this integral is stuck.

The Fourier transform offers a way around this. Write the payoff as its own inverse transform, $g(y) = \frac{1}{2\pi}\int \hat g(\omega)\, e^{i\omega y}\, d\omega$, and exchange the order of integration. The density now meets the exponential, and that pairing is the characteristic function,

$$\int_{-\infty}^{\infty} e^{i\omega y}\, q(y)\, dy = \mathbb{E}\!\left[e^{i\omega y}\right] = \phi(\omega),$$

the Fourier transform of the density. So the price turns into an integral against $\phi$ rather than $q$,

$$V = \frac{e^{-rT}}{2\pi}\int_{-\infty}^{\infty} \hat g(\omega)\,\phi(\omega)\, d\omega.$$

Heston withholds the density but provides $\phi$ in closed form, so this last integral is computable and the density is never reconstructed. The Heston characteristic function and the mechanics of this inversion are the subject of a future article on stochastic volatility.

---

## Reflection

Working through this, there are three things I find really fascinating.

The first is Fourier's original leap. To take an arbitrary profile, a temperature with a corner in it, a pulse with outright jumps, and assert that it is the sum of infinitely many smooth, endless sine and cosine waves is not an obvious thing to believe. His contemporaries did not believe it; the claim that discontinuous shapes could be assembled from perfectly smooth ones met skepticism from mathematicians as formidable as Lagrange. That Fourier saw it anyway, and was essentially right, is a conviction I find hard to reconstruct after the fact. How he came to see a function as a spectrum is the part I cannot quite explain, only admire.

The second is what happens to orthogonality in the limit. On the finite rod, two waves of unrelated frequency are not orthogonal: their overlap over the interval is some bounded number that doesn't vanish. Yet let the rod grow without bound and that overlap collapses to exactly zero, and distinct frequencies become perfectly orthogonal. Passing to a limit can lead to a new property, qualitatively different from the finite case.

The third is the one I keep returning to. The general solution of the heat equation on the infinite line, which we built by decomposing, evolving, and reassembling, can also be written with no Fourier machinery at all, as the expected value of the initial profile sampled at the endpoint of a diffusing particle. Two pictures that share no obvious vocabulary, a deterministic field smoothing itself and a single random walker wandering, produce the identical formula, because the heat kernel is both the spreading weight of diffusion and the transition density of Brownian motion. Feynman-Kac is the statement that these are the same object, and it is what lets an option price be read either as a PDE to solve or an expectation to take. That a hard analytic problem can always be recast as an average over random paths still strikes me as quietly powerful.

A last point on where this method applies, since it is not universal. The transform needs three things: a linear equation, constant coefficients, and an unbounded domain. The heat equation has all three, which is why it yielded cleanly; Black-Scholes gains them once the log substitution makes its coefficients constant. Many real pricing problems do not, and there the clean decomposition into independent frequencies breaks down. That is where the finite difference method takes over, handling what the transform cannot at the price of approximating rather than solving in closed form.

---

## Appendix: Solving Black-Scholes with the Fourier Transform

The value $V(S,t)$ of a European call satisfies

$$\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} + rS\frac{\partial V}{\partial S} - rV = 0,$$

with the terminal condition that at expiry the option is worth its payoff, $V(S,T) = \max(S - K, 0)$.

We apply the following transformations to reduce it to the heat equation.

1. **Log price, $x = \ln S$.** This makes the variable coefficients constant: the powers of $S$ in $\tfrac{1}{2}\sigma^2 S^2$ and $rS$ cancel the powers produced when a derivative hits a log, since $S\,\partial_S = \partial_x$ and $S^2\partial_{SS} = \partial_{xx} - \partial_x$.
2. **Backward time, $\tau = T - t$.** An option is anchored by its payoff at expiry, while the heat equation runs from an initial condition, so we reverse time and let $\tau$ measure time to expiry. This flips the sign, $\partial_t = -\partial_\tau$.
3. **Exponential factor, $V = e^{ax + b\tau}\,u$.** After the first two steps the equation still carries a first-derivative (drift) term and the $-rV$ term. This factor clears both: the spatial part $e^{ax}$ removes the drift, and the time part $e^{b\tau}$ removes the discount term together with the constant that killing the drift leaves behind, so $b$ is more than the bare rate. The required values are

$$a = -\frac{r - \tfrac{1}{2}\sigma^2}{\sigma^2}, \qquad b = -r - \frac{(r - \tfrac{1}{2}\sigma^2)^2}{2\sigma^2}.$$

With these the equation becomes exactly the diffusion we solved on the infinite rod,

$$\frac{\partial u}{\partial \tau} = \frac{1}{2}\sigma^2 \frac{\partial^2 u}{\partial x^2}, \qquad \kappa = \tfrac{1}{2}\sigma^2.$$

At $\tau = 0$ the factor is $e^{ax}$ and $V(S,T)$ is the payoff, so the initial condition is $u(x, 0) = e^{-ax}\max(e^{x} - K, 0)$.

This is the unbounded heat equation, so its solution is the initial profile convolved with the heat kernel of variance $2\kappa\tau = \sigma^2\tau$, the Gaussian we assembled frequency by frequency in the body, centered at $x$,

$$u(x,\tau) = \int_{-\infty}^{\infty} e^{-ay}\max(e^{y} - K, 0)\,\frac{1}{\sqrt{2\pi\sigma^2\tau}}\,e^{-(y - x)^2/(2\sigma^2\tau)}\, dy.$$

The payoff vanishes below $y = \ln K$, so the lower limit becomes $\ln K$ and the integrand splits into two exponential-times-Gaussian pieces,

$$u(x,\tau) = \int_{\ln K}^{\infty} \big[\,e^{(1-a)y} - K e^{-ay}\,\big]\,G(y)\,dy, \qquad G(y) = \frac{1}{\sqrt{2\pi\sigma^2\tau}}\,e^{-(y - x)^2/(2\sigma^2\tau)}.$$

Each piece has the form $\int_{\ln K}^{\infty} e^{cy}\, G(y)\,dy$, which completing the square evaluates as

$$\int_{\ln K}^{\infty} e^{cy}\, G(y)\,dy = e^{\,cx + \frac{1}{2}c^2\sigma^2\tau}\,N\!\left(\frac{x - \ln K + c\sigma^2\tau}{\sigma\sqrt\tau}\right),$$

applied with $c = 1 - a$ and $c = -a$. The two arguments simplify, using $(1-a)\sigma^2 = r + \tfrac{1}{2}\sigma^2$ and $-a\sigma^2 = r - \tfrac{1}{2}\sigma^2$, to exactly

$$d_{1,2} = \frac{\ln(S/K) + \big(r \pm \tfrac{1}{2}\sigma^2\big)\tau}{\sigma\sqrt\tau}.$$

Restoring $V = e^{ax + b\tau}u$, the chosen $a, b$ make the exponential constants in front of each term simplify, since $b + \tfrac{1}{2}(1-a)^2\sigma^2 = 0$ and $b + \tfrac{1}{2}a^2\sigma^2 = -r$, leaving the two terms as $S\,N(d_1)$ and $K e^{-r\tau} N(d_2)$. With $\tau = T - t$,

$$V(S,t) = S\,N(d_1) - K e^{-r(T-t)}\,N(d_2).$$

The Black-Scholes call price, recovered by reducing the equation to the heat equation and averaging the payoff against its kernel.

[^density]: All the heat sits at a single point of zero width. Reading $u$ as a temperature would put finite heat into zero width, which is infinite at that point. As a density it is fine: an infinitely tall spike whose integral, the total heat, is still the finite $Q$.

[^gaussint]: The identity $\int_0^\infty e^{-a\omega^2}\cos(b\omega)\,d\omega = \tfrac{1}{2}\sqrt{\pi/a}\,e^{-b^2/(4a)}$ can be derived by differentiation under the integral sign (Feynman's trick): differentiating in $b$ and integrating by parts turns the integral into a first-order differential equation in $b$, which solves to the Gaussian above.