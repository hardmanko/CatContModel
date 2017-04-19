# Examples

In this directory, there are some examples of the use of this package with simulated data sets. There are comments throughout the examples to help you follow the logic of the steps.

In most cases, the data is simulated from the model that is then fit to that data. The data is simulated with known parameter values, so it is possible to examine how well the model is able to recover the true parameter values. I use this for testing the package.

In some cases, in the example, a model other than the model that generated the data is fit to the data to test whether the incorrect model fits worse than the correct model.

# List of Examples

betweenItem: Basic examples of fitting the between-item model variant, testing hypotheses, making plots, etc.

withinItem: Same as betweenItem but for the within-item model variant.

ZL: Same as betweenItem but for the ZL (Zhang & Luck, 2008) model.

linear_betweenItem: An example of using linear data with the between-item model variant. Using a different model variant is straightforward, this is just focused on the linear data aspect.

factorial_betweenItem: An example with a factorial design (two factors).

betweenParticipants_oneFactor: A between-participants design with one factor is analyzed, at least in part. This example only shows how to fit the model in a between-participants case and test a main effect of the single factor.

betweenAndWithin_constrained: Uses the between-and-within model variant with a constraint on the within-item component of the model that helps to identify the parameters. In general, I do not recommend the between-and-within model variant due to issues with parameter identifiability.