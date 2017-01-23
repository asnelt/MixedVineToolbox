Mixed Vine Toolbox for Matlab
=============================

Toolbox for canonical vine copula trees with mixed continuous and discrete
margins. If you use this toolbox, then please cite:
A. Onken and S. Panzeri (2016). Mixed vine copulas as joint models of
spike counts and local field potentials. In D. D. Lee, M. Sugiyama,
U. V. Luxburg, I. Guyon and R. Garnett, editors, Advances in Neural
Information Processing Systems 29 (NIPS 2016), pages 1325–1333.


Description
-----------

In this toolbox, we implemented a complete framework based on canonical
vine copulas for modelling multivariate data that are partly discrete and
partly continuous. The resulting multivariate distributions are flexible
with rich dependence structures and arbitrary margins. For continuous
margins, we provide implementations of the normal and the gamma
distributions. For discrete margins, we provide the Poisson, binomial and
negative binomial distributions. As bivariate copula building blocks, we
provide the Gaussian, student and Clayton families as well as rotation
transformed Clayton families. The toolbox includes methods for sampling,
likelihood calculation and inference, all of which have quadratic
complexity. These procedures are combined to estimate entropy and mutual
information by means of Monte Carlo integration.


Demonstration
-------------

The script `demo.m` demonstrates how to apply the Mixed Vine Toolbox. It
constructs a 4D mixed canonical vine with normal, gamma, Poisson and
binomial margins and builds the vine tree from Gaussian, Student, Clayton
and rotated Clayton copula families. It calculates and plots multivariate
marginal probability densities, samples from the distribution, estimates
the model from the samples and calculates entropy.


Functions
---------

* mixedvinefit - Mixed copula vine estimates
* mixedvinepdf - Mixed copula vine probability density function
* mixedvinernd - Mixed copula vine random numbers
* mixedvineentropy - Mixed copula vine entropy estimate
* mixedvineinfo - Mixed copula vine mutual information estimate
* mixedgaussfit - Mixed copula vine estimates with Gaussian copula
* marginfit - Univariate margin estimates
* marginpdf - Univariate margin probability density function
* margincdf - Univariate margin cumulative distribution function
* margininv - Inverse of univariate margin CDF
* copulafit - Copula parameter estimates
* copulapdf - Copula probability density function
* copulacdf - Copula cumulative distribution function
* copulaccdf - Copula conditional cumulative distribution function
* copulaccdfinv - Inverse of copula conditional CDF
* copularnd - Copula random numbers


License
-------

```text
Copyright (C) 2016 Arno Onken

This file is part of the Mixed Vine Toolbox.

The Mixed Vine Toolbox is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as published
by the Free Software Foundation; either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, see <http://www.gnu.org/licenses/>.
```
