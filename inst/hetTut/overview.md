# Overview

*Steve E. Bellan, 2012-2016*  
*Some Rights Reserved*  
*[CC BY-NC 3.0](http://creativecommons.org/licenses/by-nc/3.0/)*

This exercise will help you build intuition for how heterogeneity in contact mixing patterns affects infectious disease dynamics.  We have built a continuous time stochastic *SIR* model for you to explore.  Importantly, in this model not all individuals are the same.  Each individual has their own contact rate, with which they meet other people in the population.  Some individuals meet other people more often, and some meet people less often.

If we think of the total distribution of contact mixing patterns, then we can characterize it using terms such as its mean, its variance, how skewed it is (*i.e.* how long is its tail).  While it is not so much important how we have built this distribution (we use a gamma distribution if you are interested), we are giving you the opportunity to change the mean and variance of individuals' contact rates. Using the code below, play with the outbreak model and see if you can understand how heterogeneity affects the disease dynamics, and specifically the distribution of outbreak sizes.
