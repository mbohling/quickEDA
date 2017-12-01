# quickEDA
R Functions for Exploratory Data Analysis.

Provides functions to quickly summarize assumptions needed in order to determine how to 
proceed in modeling your data using exploratory data analysis.  Can also be used to validate 
your model, making sure that your residuals are in statistical control.

All functions expect a numerical vector.

eda_fourPlot(y)

eda_locationDrift(y)

eda_varianceDrift(y)

eda_randomness(y, plot.it = FALSE) #Plot it = TRUE will plot Lag 1 P/ACF.

eda_normality(y)

These functions can help jumpstart your analysis by providing initial visualizations of your
data, as well as check for basic assumptions underlying a simple univariate regression model.

eda_SummaryStats(y)

Computes and prints summary statistics of your data vector.

eda_plot(y, type = c("run", "lag", "hist", "qqn", "spectrum")

Easily plot your data using the following visualization formats.  Just specify the type.

Download with:
devtools::install_github('minbad/quickEDA') #may need force = TRUE

Still under heavy development and this is my first R package so go easy on me! :)
