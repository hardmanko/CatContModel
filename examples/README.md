# Examples

In this directory, there are some examples of the use of this package with simulated data sets. There are comments throughout the examples to help you follow the logic of the steps.

In most cases, the data is simulated from the model that is then fit to that data. The data is simulated with known parameter values, so it is possible to examine how well the model is able to recover the true parameter values. I use this for testing the package.

In some cases, in the example, a model other than the model that generated the data is fit to the data to test whether the incorrect model fits worse than the correct model.

# Usage

For each example, you should download at least the files with names ending "_analysis.R" and "_data.txt". You could just download the whole examples directory.

The data are simulated with known parameter values. The parameter values and data are sampled in files with names ending in "_simulateData.R". Those files depend on "DataSimulatingFunctions.R".

The true parameter values are stored in files ending in "_parameters.txt".

# List of Examples

## The Four Model Variants

*`betweenItem`*: Basic examples of fitting the between-item model variant, testing hypotheses, making plots, etc.

*withinItem*: Same as `betweenItem` but for the within-item model variant.

*ZL*: Same as `betweenItem` but for the ZL (Zhang & Luck, 2008) model.

*betweenAndWithin_constrained*: Uses the between-and-within model variant with a constraint on the within-item component of the model that helps to identify the parameters. Basically, `pContWithin` is forced to be 0.5. In general, I do not recommend the between-and-within model variant due to issues with parameter identifiability. This example shows the only way I have found to make the between-and-within variant identifiable.

## Linear Data

*linear_betweenItem*: An example of using linear data with the between-item model variant. Using a different model variant with linear datais straightforward; this example is just focused on the linear data aspect.

## Within-Participants Factorial Design

*factorial_betweenItem*: An example with a factorial design (two factors). This example is focused on some of the complex analyses that you can do with factorial designs.

## Between-Participants and Mixed Designs

*betweenParticipants_oneFactor*: A between-participants design with one factor is analyzed, at least in part. This example only shows how to fit the model in a between-participants case and test a main effect of the single factor.

*factorialMixedDesign*: A balanced factorial design with both within-participants and between-participants factors.

*unbalancedFactorialMixedDesign*: Same as `factorialMixedDesign`, but the design is unbalanced. In particular, there are three between-participants groups but with different levels of the within-participants factors for each group.


