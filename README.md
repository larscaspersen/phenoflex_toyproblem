
<!-- README.md is generated from README.Rmd. Please edit that file -->

# An example of the current fitting routine using PhenoFlex

This is an example with flowering data of apple cultivar “Blanquina” for
the period 1987 - 1997 and 2004 - 2020 in Asturias, Spain. Hourly
weather data was generated from daily temperature observations using the
chillR function `stack_hourly_temps()`. Phenology and temperature data
can be found in the “data” folder.

In the “code” folder you can find a script called “example_fit.R”. It
contains the usual routine we use to prepare the data and then carry out
the first (of usually more than 10) round of fitting using the function
`phenologyFitter()` fromt the chillR package. The function requires
several arguments:

-   **par.guess** contains the initial guess for the parameters. The
    order should be (yc, zc, s1, Tu, E0, E1, A0, A1, Tf, Tc, Tb, slope)
-   **modelfn**, this is the wrapper function for the phenoflex model.
    It can also be customised, so that some parameters are for example
    kept constant.
-   **bloomJDays** is a vector of the Juliand days of the flowering
-   **SeasonList** is a list with the hourly temperature data organized
    in the same order as the provided phenology data
-   **lower** contains the lower bounds of the model parameters
-   **upper** contains the upper bounds of the model parameters
-   **seed** to make results reproducible
-   **control**, a list of further variables directly passed to the
    simulation annealing function currently used in phenoflex
    (`GenSA::GenSA()`). Most interesting is “maxit” which defines the
    maximum number of iterations. “max.time” could be also of interest,
    because it specifies the time the optimizer should run in seconds,
    which can be used instead of “max.it”

The fitting result is checked on RMSE. The temperature response curve
resulting from the new set of parameters is inspected using some helper
functions, which you can find in the file “helper_functions.R”.
