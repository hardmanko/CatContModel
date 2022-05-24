# Examples

In this directory, there are some examples of the use of this package with simulated data sets. There are comments throughout the examples to help you follow the logic of the steps.

In most examples, the data is simulated from a model variant, then that model variant is fit to the data. The data is simulated with known parameter values, so it is possible to examine how well the model is able to recover the true parameter values. These examples are important for testing that the package is working.

In some examples, a model other than the model that generated the data is fit to the data to test whether the incorrect model fits worse than the correct model. This is a good test of model variant specificity.

# Usage

For each example, you should download at least the files with names ending "_analysis.R" and "_data.txt". It may be easiest to download the whole repository to get the whole examples directory (it's not included in the installable package files in releases).

The data are simulated with known parameter values. The parameter values and data are sampled in files with names ending in "_simulateData.R". Those files depend on "DataSimulatingFunctions.R".

The true parameter values that were used to sample the data are stored in files ending in "_parameters.txt".

# List of Examples

The most fully documented example is `betweenItem`, so start there.

## The Four Model Variants

`betweenItem`: Basic examples of fitting the between-item model variant, testing hypotheses, making plots, etc.

`withinItem`: Same as `betweenItem` but for the within-item model variant.

`ZL`: Same as `betweenItem` but for the ZL (Zhang & Luck, 2008) model.

`betweenAndWithin_constrained`: Uses the between-and-within model variant with a constraint on the within-item component of the model that helps to identify the parameters. Basically, `pContWithin` is forced to be 0.5. In general, I do not recommend the between-and-within model variant due to issues with parameter identifiability. This example shows the only way I have found to make the between-and-within variant identifiable.

## Within-Participants Factorial Design

`factorial_betweenItem`: An example with a factorial design (two factors). This example is focused on some of the complex analyses that you can do with factorial designs.

## Between-Participants and Mixed Between/Within Designs

`factorialMixedDesign`: A balanced factorial design with both within-participants and between-participants factors.

`unbalancedFactorialMixedDesign`: Same as `factorialMixedDesign`, but the design is unbalanced. In particular, there are three between-participants groups but with different levels of the within-participants factors for each group.

`betweenParticipants_oneFactor`: A between-participants design with one factor is analyzed, at least in part. This example only shows how to fit the model in a between-participants case and test a main effect of the single factor.

## Parallel Parameter Estimation

`parallel`: Demonstrates running parallel chains using `runParameterEstimation_parallel` and `continueSampling_parallel`. This method takes samples more quickly than running a single chain and makes it easy to do multi-chain convergence diagnostics.

## Linear Data

For the purposes of CatContModel, linear data is bounded on both ends versus circular data which wraps around.

`linear_betweenItem`: An example of using linear data with the between-item model variant. Using  model variant with linear data is straightforward; this example is just focused on the linear data aspect.

## Half-Circle Data

Some WM tasks use a line or bar that rotates around a central point. Participants estimate the angle of the line. 
This angle is in the range [0,180] because the line points in two directions at once. These examples give some ideas for dealing with half-circle data.

`halfCircle_betweenItem`: Contains two methods for fitting half-circle data that contains categories.

`halfCircle_fullCircle`: This example is for a special case in which 1) data from a half-circle task is to be compared to a full-circle task and 2) there is no/negligible categorical responding in the data so the Zhang & Luck (2008) model can be used.