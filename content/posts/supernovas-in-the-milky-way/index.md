---
weight: 4
title: "Supernova rates in the Milky Way"
date: 2025-10-10T09:23:00-03:00
draft: false
author: "f"
authorLink: "https://fsb2210.github.io"
description: ""
images: []
math: true
resources:
- name: "featured-image"
  src: "featured-image.png"

tags: ["physics", "astronomy"]
categories: ["physics", "astronomy"]

lightgallery: true

plotly: true
---

In this article we explore how astronomers derive the rate of supernovae in the Milky Way and compare it to actual
observations of such events throughout the history. We will try to understand why astronomers often cite a rate of 2
supernovae per century for the Milky Way, even though we haven't seen one since 1604.
<!--more-->

# Supernovae rates in the Milky Way

In the vast of the galaxy we happen to live in, the Milky Way, stars don't just simply fade away; some go out with a
monumental bang, and we call them *supernova* (from the Latin word *nova* - meaning "new") because when they appear
in the sky they can briefly shine as brightly as an entire galaxy, looking like a "new star" where none was seen
before.

A supernova marks a dramatic turning point in a star’s life. It can be triggered when a massive star exhausts its
nuclear fuel and its core collapses under its own gravity, or when a dense stellar remnant called a white dwarf grows
unstable after stealing matter from a companion. In either case, the explosion releases enormous amounts of energy,
forges many of the heavy elements found in planets and living beings, and sends shock waves racing through interstellar
space.

Such events are not just theoretical: humans have witnessed supernovae long before modern astronomy existed. There
exist ancient records from different civilizations describing sudden, bright stars appearing in the night sky and
remaining visible for weeks or months.

But how often do these explosions happen in our own galaxy? To answer that, astronomers talk about rates: an estimate
of how frequently supernovae occur over time.

In the sections that follow, we will look at historical observations of supernovae and see how they help us infer the
rate at which the Milky Way produces these extraordinary events.

## Historical records of supernovae

Historically, observations of such events are extremely difficult to perform due to our location in the Galaxy being in
the Galactic disk which produces a lot of obscuration along the line of sight. Despite that fact, there are some
records in the last millenium of naked eye observations of "new stars", which happen to be supernovae.

In fact, there are five "historical" supernovae from the years: 1006, 1054, 1181, 1572 and 1604. The first three of
these events were found in records of civilizations from China and Japan while the latter two were mostly followed by
Europeans. In the case of the supernova from 1572, records comes from Tycho Brahe, and for the 1604 supernovae we have
records from Johannes Kepler.

It is now known that these supernovae are located nearby in the Galaxy, within a few kilo-parsecs away from us.

We can see here that we can actually infer a rate of these events to be close to one every two centuries which is not
that different from the one we will derive later on of around two per century. More on that later!


{{< admonition note "Note" >}}
A thorough description of historical records of supernovae can be found the **Handbook of Supernovae**, from which I
extracted the year of their recordings.
{{< /admonition >}}

## Supernovae rates

In order to make estimations on what we expect to see in the future, that is make predictions to check whether our
theories are correct or need to be changed, we need to compute the frequency of supernovae events in the Milky Way.
In order to accomplish such task, we are going to dive into the world of statistics, specifically
**Bayesian inference** and how it is applied.

### Bayes theorem

At the heart of Bayesian inference lies Bayes’ theorem, a very powerful rule that tells us how to update our beliefs in
light of new information. In astronomy, and particularly when dealing with rare events like supernovae, this framework
is especially valuable because observations are often incomplete, noisy, or indirect.

Mathematically, it is written as

\begin{equation*}
    P(\vartheta \\ | \\ \\{h\\}) = \frac{P(\\{h\\} \\ | \\ \vartheta) P(\vartheta)}{P(\\{h\\})}
\end{equation*}

Bayes’ theorem combines three ingredients: what we initially believe about one or more quantities, the *prior*
$P(\vartheta)$, how compatible the observations are with different possible values of those quantities, the *likelihood*
$P(\\{h\\} \\ | \\ \vartheta)$, and the updated belief after seeing the data, the *posterior* $P(\vartheta \\ | \\ \\{h\\})$.
Finally, $P(\\{h\\})$ is the evidence, acting as a normalization constant.

{{< admonition note "Note" >}}
For supernova rates, this means starting with an initial expectation of how frequently these explosions occur in the
Milky Way and then refining that expectation using historical records and modern observations.
{{< /admonition >}}

### Likelihood

The likelihood encodes how the data relate to the model parameters we want to estimate. It answers the question: if a
given supernova rate were true, how likely would it be to observe the data we have? Thus, it forces us to be explicit
about how supernovae occur in time, how they are detected, and what kinds of observations count as evidence.

In practice, defining the likelihood means translating our physical and observational assumptions into a mathematical
form. For example, we must decide how to model the timing of supernova explosions, how long they remain observable, and
how observational limitations affect what we actually record. Once specified, the likelihood acts as the bridge between
theory and data.

{{< admonition note "Note" >}}
Supernovae can be considered *point events* on time. They are independent: the ocurrence of one does not make the next
more likely. This behaves like a [**Poisson process**](https://en.wikipedia.org/wiki/Poisson_point_process). Thus, if
we define $\lambda$ as the rate of supernovae per century, the probability of observing $k$ events in $t$ centuries is:
\begin{equation*}
    P(k \\ | \\ \lambda, t) = \dfrac{(\lambda t)^k e^{− \lambda t}}{k!}
\end{equation*}
{{< /admonition >}}

However, there is a again the problem of dust hiding supernovae from us, they are *hidden*. Thus, the *observed rate*
isn't $\lambda$: it's $\lambda_{\rm MW} \cdot p$ with $p$ the detection probability. Treating both the actual rate
$\lambda_{\rm MW}$ and the detection probability $p$ as unknowns, we would be accounting for blind spots!

We use the **Poisson distribution** to ask the model: Given our guess for the rate and our limited view, how likely is
it that we saw exactly 7 events? (5 supernovae seen with the naked eye + 2 remnants in other ranges of the
electromagnetic spectrum).

### Priors: encoding scientific constraints

Priors represent our knowledge—or lack of— about a parameter before considering the current data. In the
context of supernova rates, a prior might reflect previous astronomical studies, theoretical expectations from stellar
evolution, or simply a statement that we are initially uncertain over a wide range of possible rates.

Priors are not arbitrary guesses; they are a transparent way of stating assumptions. They can be informative, when
strong external knowledge exists, or weakly informative (or even nearly uninformative) when we want the data to speak
more loudly. An important part of Bayesian analysis is examining how sensitive the results are to the choice of prior,
ensuring that conclusions about supernova rates are driven by observations rather than hidden assumptions.

In this case, we will assume that the rate is likely around 2, but let the model explore anything from 0 to 5, that is,
we use a **uniform distribution**.
\begin{equation*}
    P(\lambda) = \dfrac{1}{5}
\end{equation*}

Moreover, since $p$ is a probability, it must be bounded between 0 and 1. Thus, we use the **beta distribution**.
\begin{equation*}
    P(p) = B(\alpha, \beta)
\end{equation*}

{{< admonition tip "Tip" >}}
By setting $\alpha = 4$ and $\beta = 2$, we create a "belief" that leans toward higher detection for remnants (which
our data relies on), while still allowing for the possibility that many are missed.
{{< /admonition >}}

## Rates model

Now we can condense what was previosuly said into code to estimate rates. For that we can create a script using the
*python* language.

First we need to import some packages and define some values:

```python
import numpy as np
import emcee
import scipy.stats as stats

# data
k_obs = 7  # 5 historical + 2 remnants
t_obs = 10
```

Now we are ready to create the logarithm of the probabilities:

```python
def log_prior(theta):
    _, p = theta
    return stats.beta.logpdf(p, 2, 8)

def log_likelihood(theta):
    lam, p = theta
    return stats.poisson.logpmf(k_obs, lam * p * t_obs)

def log_probability(theta):
    lam, p = theta
    # boundary constraints
    if not (0 < lam < 5 and 0.01 < p < 1): return -np.inf
    lp = log_prior(theta)
    if not np.isfinite(lp): return -np.inf
    ll = log_likelihood(theta)
    return lp + ll
```

Finally, we are ready to compute the posterior probabilities. In this case we are going to use the [emcee]()
package that allows us to run an Markov Chain of MonteCarlo in few lines of code:

```python
# MCMC Setup
initial = [2.0, 0.5]
nwalkers, ndim = 32, 2
pos = initial + 1e-4 * np.random.randn(nwalkers, ndim)
nsteps = 1000
burnin = 200

with Pool() as pool:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, pool=pool)
    sampler.run_mcmc(pos, nsteps, progress=True)
samples = sampler.get_chain(discard=burnin, flat=True)
```

We can now make some graphics to understand what we computed

{{< plotly json="/json/mcmc_corner.json" >}}

The graphic is clear: the estimated number of supernovae in the Galaxy is not a sharply peaked towards 2 but actually
is broader given the historical records. The explanation is as follows: seeing 7 events through a "dusty window"
implies there were likely many more supernovae that we missed. The 400 year gap of no observed supernovae with the
naked eye is not physical, it ss just a statistical lull.

## Forecasting supernovae in the Galaxy

We can use our posteriors on the $\lambda_{\rm MW}$ to forecast the chances of a supernova happening *in the Galaxy* in
the next 10, 25 and 50 years:

```python
# extract samples from the sampler
# samples[:, 0] is mu (global average)
# samples[:, 1] is lam_mw (Milky Way specific)
param_names = ["lambda", "p"]
lam_mw_samples = samples[:, 0]
p_samples = samples[:, 1]

# forecast
for yrs in [10, 25, 50]:
    t = yrs / 100
    # probability of at least one event: 1 - P(zero events)
    prob_dist = 1 - np.exp(-lam_mw_samples * t)
    mean_prob = np.mean(prob_dist) * 100
    print(f"- chance of a supernova in the next {yrs} years: {mean_prob:.1f}%")
```

Using the above we find that:

- chance of a supernova in the next 10 years: 27.1%
- chance of a supernova in the next 25 years: 53.6%
- chance of a supernova in the next 50 years: 76.9%

If we are lucky, we can maybe see one in the near future! (Betelgeuse, maybe?)

## Final thoughts

In this post we have seen that the estimated number of supernovae in the Galaxy should not be only guided by what we
see but what we cannot see; the Milky Way is huge, dusty and sure knows how to hide exploding stars.

Reaching these results was possible given the "humble" way of estimating rate of events in the context of Bayesian
statistics.

### On how to improve the previous estimates: hierarchical modeling

Actually, to be more scientifically accurate there are two problems that we have not mentioned earlier:

1. there is only one Milky Way (small sample size) and,
2. dust blocks a lot of the optical view.

To solve this, astronomers typically assume that the rate of supernovae in the Milky Way ($\lambda_{\rm MW}$) is drawn
from a broader distribution of rates ($\mu$) found in similar spiral galaxies and use them as contraints in the *prior*
distributions.

Given the complexity of such approach, which requires to investigate about such events in different galaxies, we leave
it for astronomers to compute rates using this method.
