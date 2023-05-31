# Part 1: Low Variance

First, let's run the model with very low variance so we can see what happens in a situation that is very similar to the homogenous population model.  This is like the compartmental box models that we have been using in differential equation and discrete time versions of the *SIR* model.

This is the code running behind the scenes:

```r
mxdst <- het.population(n = 100, beta.mean = 2, beta.var = 0.001)
het.hist(mxdst, beta.mean)
het.epidemic(mxdst, runs = 5, end.time = 10, gmma = 1)
```

To add another sample series, click the button above the graphs.

## EXPLANATION OF VARIABLES

 - `beta.mean` is the mean of your contact rate distribution
 - `beta.var` is its variance
 - `end.time` is how long to let the simulations go before cutting them off (if the outbreak hasn't already died)
 - `pop.size` is the total population size
 - `gmma` is the recovery rate (or 1/mean duration of infectiousness)
 - `rho` is the waning rate (or 1/mean duration of immunity)
 
## OUTPUT

The panel on the left is a histogram of the contact rate distribution. The middle shows the epidemic time series. The right panel shows the distribution of outbreak final size (*i.e.* total number of people infected before outbreak is over) from the runs.
