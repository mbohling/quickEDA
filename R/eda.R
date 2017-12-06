#' EDA Assumptions Testing - Four Plot
#'
#' This function takes in a vector of numerical data and creates a 4-Plot
#' consisting of:
#' 1. Run Sequence Plot
#' 2. Lag 1 Plot
#' 3. Histogram / Density Plot
#' 4. QQ Normality Plot
#'
#' @param y Numeric vector.
#' @examples
#' eda_fourPlot(rnorm(1000))
#' @export
eda_fourPlot <- function(y) {

  if(!is.numeric(y)) {
    stop("Error: Input vector must be numeric.")
  }

  n = length(y)
  t = 1:n

  par(mfrow = c(2, 2),
      oma = c(0, 0, 2, 0))

  plot(t, y , ylab="Y", xlab="Run Sequence", type="l", main = "Run Sequence")
  plot(y, Lag(y), xlab="Y[i-1]", ylab="Y[i]", main = "Lag Plot")
  h <- hist(y, xlab="Y", main = "Histogram/Density Overlay")
  xfit <- seq(min(y), max(y), length.out = length(y))
  yfit <- dnorm(xfit, mean = mean(y), sd = sd(y))
  yfit <- yfit * diff(h$mids[1:2]) * length(y)
  lines(xfit, yfit, col="blue", lwd=2)
  qqnorm(y, main = "QQ Normality Plot")
  qqline(y)
  mtext("Y: 4-Plot", line = 0.5, outer = TRUE)

  par(mfrow = c(1, 1))
}

#' EDA Assumptions Testing - Drift in Location
#'
#' This function takes in a vector of numerical data and tests drift in location using
#' a simple regression model
#'
#' @param y Numeric vector.
#' @examples
#' eda_locationDrift(rnorm(1000))
#' @export
eda_locationDrift <- function(y) {

  if(!is.numeric(y)) {
    stop("Error: Input vector must be numeric.")
  }

  n = length(y)

  x = c(1:n)
  mod <- lm(y ~ 1 + x)
  ## Critical value to test that the slope is different from zero.
  slopeDifZeroTestResult = if_else(coef(summary(mod))["x","t value"] > qt(.95,n-2), "YES", "NO")

  cat("
          Test Index Fit Straightline for Non-Zero Slope
      ")

  if(slopeDifZeroTestResult == "NO") {
    cat("
        slope coefficient: ", coef(summary(mod))["x","Estimate"],"
        t value for slope: ", coef(summary(mod))["x","t value"], "< Critical Value: ", qt(.95,n-2),"
        Test for Location Drift : PASS

        ")
  }
  else {
    cat("
        slope coefficient: ", coef(summary(mod))["x","Estimate"],"
        t value for slope: ", coef(summary(mod))["x","t value"], "> Critical Value: ", qt(.95,n-2),"
        Test for Location Drift : FAIL

        ")
  }
}

#' EDA Assumptions Testing - Drift in Variation
#'
#' This function takes in a vector of numerical data and tests homoscedasticity
#'
#' @param y Numeric vector.
#' @param testType Variance Homogeniety Test to use.  Defaults to "barlett", other option is "levene".
#' @examples
#' eda_varianceDrift(rnorm(1000))
#' eda_varianceDrift(runif(1000), testType = "levene")
#' @export
eda_varianceDrift <- function(y, testType = "barlett") {

  if(!is.numeric(y)) {
    stop("Error: Input vector must be numeric.")
  }

  n = length(y)

  findGroupSplit <- function(n) {
    for (i in 4:10) {
      if(n %% i == 0) {
        foundFactor <- TRUE
        return(i)
      }
      i <- i + 1
    }
    return(4)
  }

  foundFactor <- FALSE
  kGroups <- 4
  while(!foundFactor) {
    kGroups <- findGroupSplit(n)
    if(!foundFactor) {
      y <- append(mean(y))
      n = length(y)
    }
  }

  if(testType == "levene") {

    int = as.factor(rep(1:kGroups,each=n/kGroups))
    lt <- leveneTest(residuals_rmOL, int, location = "median", bootstrap = TRUE, kruskal.test = TRUE)
    diffVariation <- lt$`F value`[1] > qf(.95, df1 = kGroups - 1, df2 = n - kGroups)
    leveneTestResult = if_else(diffVariation, "YES", "NO")

    print(lt)

    if(leveneTestResult == "YES") {
      cat("
          F-Value: ", lt$`F value`[1], "> Critical Value: ", qf(.95, df1 = kGroups - 1, df2 = n - kGroups),"
          Test for Variation Drift : FAIL
          ")
    } else {
      cat("
          F-Value: ", lt$`F value`[1], "< Critical Value: ", qf(.95, df1 = kGroups - 1, df2 = n - kGroups),"
          Test for Variation Drift : PASS

          ")
    }
  }
  else {
    ## Generate arbitrary interval indicator variable and
    ## run Bartlett's test for drift in variation.

    int = as.factor(rep(1:4, each = length(y)/4))
    bartlett <- bartlett.test(y ~ int)
    bartlettTestResult = if_else(bartlett$statistic < qchisq(.95, bartlett$statistic), "NO", "YES")

    print(bartlett)

    if(bartlettTestResult == "NO") {
      cat("
          K-Squared: ", bartlett$statistic, "< Critical Value: ", qchisq(.95, bartlett$statistic),"
          Test for Variation Drift : PASS

          ")
    } else {
      cat("
          K-Squared: ", bartlett$statistic, "> Critical Value: ", qchisq(.95, bartlett$statistic),"
          Test for Variation Drift : FAIL

          ")
    }
  }
}

#' EDA Assumptions Testing - Randomness
#'
#' This function takes in a vector of numerical data and tests for patterns in the data.
#'
#' @param y Numeric vector.
#' @param plot.it Boolean value that controls whether autocorrelation / partial autocorrelation plots are created. Defaults to TRUE.
#' @examples
#' eda_randomness(rnorm(1000))
#' eda_randomness(rnorm(1000), plot.it = FALSE)
#' @export
eda_randomness <- function(y, plot.it = TRUE) {

  if(!is.numeric(y)) {
    stop("Error: Input vector must be numeric.")
  }

  n = length(y)

  b_AC = acf(y, plot = FALSE)

  # if(plot.it == TRUE) {
  #   confIntervals = c(.90,.95)
  #   plot(b_AC, confIntervals, main="Randomness Check", ylab="Autocorrelation")
  # }

  a = acf(y, lag.max = 100, plot = FALSE)
  ac = round(a$acf[2],5)

  if(plot.it == TRUE) {
    ## Plot the lag 1 P/ACF.
    acf2(y, max.lag = 100)
  }

  # Generate Significance Value
  corr <- acf(y, lag.max=n/4,ci=c(.95,.99),main="", plot = FALSE)
  sig_level <- qnorm((1 + 0.95)/2)/sqrt(corr$n.used)
  autoCorrelationTestResult = if_else(abs(ac) < sig_level,"YES", "NO")

  cat("
          Test Autocorrelation [2] below Significant Value
      ")

  if(autoCorrelationTestResult == "YES") {
    cat("
Autocorrelation: ", abs(ac), "< Significant Value: ", sig_level,"
Test for Random Data : PASS

        ")
  } else {
    cat("
Autocorrelation: ", abs(ac), "> Significant Value: ", sig_level,"
Test for Random Data : FAIL

        ")
  }

  ## Load the lawstat library and perform runs test for Randomness
  rtest <- runs.test(y)
  runsTestResult = if_else(abs(rtest$statistic) < pnorm(.95),"YES", "NO")

  print(rtest)

  if(runsTestResult == "YES") {
    cat("
        Statistic: ", abs(rtest$statistic), "< Significant Value: ", pnorm(.95),"
        Test for Random Data : PASS

        ")
  } else {
    cat("
        Statistic: ", abs(rtest$statistic), "> Significant Value: ", pnorm(.95),"
        Test for Random Data : FAIL

        ")
  }
}

#' EDA Assumptions Testing - Normality
#'
#' This function takes in a vector of numerical data and tests if the data
#' came from a normal distribution.
#'
#' @param y Numeric vector.
#' @examples
#' eda_normality(rnorm(1000))
#' @export
eda_normality <- function(y){

  if(!is.numeric(y)) {
    stop("Error: Input vector must be numeric.")
  }

  n = length(y)

  x = qqnorm(y, plot.it = FALSE)
  ppcc <- cor(x$x,y)
  ppccTestResult = if_else(ppcc < .987, "NO", "YES")

  cat("
        Probabilty Plot Correlation Coefficients (PPCC) Normality Test
      ")

  if(ppccTestResult == "YES") {
    cat("
        PPCC: ", ppcc, "> Significant Value: ", .987,"
        Test for Random Data : PASS

        ")
  } else {
    cat("
        PPCC: ", ppcc, "< Significant Value: ", .987,"
        Test for Random Data : FAIL

        ")
  }

  # Run Anderson-Darling test for Normality
  ADCriticals <- data_frame(Alpha = c(.1, .05, .025, .01),
                            Normal_log = c(.631, .752, .873, 1.035),
                            Weib_Gumb_exp = c(.637, .757, .877, 1.038))

  ad <- ad.test(y)
  andersonDarlingTestResult = if_else(abs(ad$statistic) < ADCriticals[2,2], "YES", "NO")

  print(ad)

  if(andersonDarlingTestResult == "YES") {
    cat("
        Test Statistic: ", abs(ad$statistic), "< Significant Value: ", unlist(ADCriticals[2,2]),"
        Test for Random Data : PASS

        ")
  } else {
    cat("
        Test Statistic: ", abs(ad$statistic), "> Significant Value: ", unlist(ADCriticals[2,2]),"
        Test for Random Data : FAIL

        ")
  }
}

#' Summary Statistics
#'
#' This function takes in a vector of numerical data and computes
#' summary statistics
#'
#' @param y Numeric vector.
#' @examples
#' eda_SummaryStats(rnorm(1000))
#' @export
eda_SummaryStats <- function(y) {

  if(!is.numeric(y)) {
    stop("Error: Input vector must be numeric.")
  }

  ## Compute summary statistics.
  ybar = round(mean(y),5)
  std = round(sd(y),5)
  n = round(length(y),0)
  stderr = round(std/sqrt(n),5)
  v = round(var(y),5)

  # Compute the five number summary.
  # min, lower hinge, Median, upper hinge, max
  z = fivenum(y)
  lhinge = round(z[2],5)
  uhinge = round(z[4],5)
  adjustedMin <- ifelse(min(y) < 0, abs(min(y)), min(y))
  rany = round((max(y) + adjustedMin),5)

  ## Compute the inter-quartile range.
  iqry = round(IQR(y),5)

  Values = c(n,ybar,std,stderr,v,rany,lhinge,uhinge,iqry)
  Statistics = c("Number of Observations ", "Mean", "Std. Dev.",
                 "Std. Dev. of Mean", "Variance", "Range",
                 "Lower Hinge", "Upper Hinge", "Inter-Quartile Range")
  cat("
      Summary Statistics:

      ")

  df <- data_frame(Statistics, Values)
  print(as_tibble(df))
}

#' ggPlot Diagnostic Plots
#'
#' This function takes in a vector of numerical data and a type of ggplot
#' in order to inspect with more detail.
#'
#' @param y Numeric vector.
#' @param type Which plot to run.  Default is "run".
#' @examples
#' eda_ggPlot(rnorm(500))
#' eda_ggPlot(rnorm(500), type = "hist")
#' @export
eda_plot <- function(y, type = c("run", "lag", "hist", "qqn", "spectrum"), binWidth = NULL) {

  if(!is.numeric(y)) {
    stop("Error: Input vector must be numeric.")
  }

  if(type == "run") {
    .plotRun(y)
  }
  else if(type == "lag") {
    .plotLag(y)
  }
  else if(type == "hist") {
    .plotHist(y)
  }
  else if(type == "qqn") {
    .plotQQNorm(y)
  }
  else if(type == "spectrum") {
    .plotSpectrum(y)
  }
  else {
    .plotRun(y)
  }
}

.plotRun <- function(y) {

  par(mfrow = c(1, 1))

  t = 1:length(y)

  ggplot(data = NULL, aes(x = t, y = y)) +
    ggtitle("Run Sequence") +
    xlab("Run") +
    ylab("Y") +
    geom_line() +
    theme_light()
}

.plotLag <- function(y) {

  par(mfrow = c(1, 1))

  t = 1:length(y)

  ggplot(data = NULL, aes(x = y, y = Lag(y))) +
    ggtitle("Lag Plot") +
    xlab("Y[i-1]") +
    ylab("Y[i]") +
    geom_point() +
    theme_light()
}

.plotHist <- function(y, binWidth = 30) {

  par(mfrow = c(1, 1))

  t = 1:length(y)

  ggplot(data = NULL, aes(x = y)) +
    ggtitle("Histogram/Density Overlay") +
    xlab(paste0("Data : binWidth = ", binWidth)) +
    ylab("Frequency") +
    geom_histogram(binwidth = binWidth) +
    suppressWarnings(geom_density(col = "blue", aes(binWidth = binWidth, y = binWidth * ..count..))) +
    theme_light()
}



.plotQQNorm <- function(y) {

  qqnorm(y,main="QQ Normality Plot")
  qqline(y)
}

.plotSpectrum <- function(y) {
  spectrum(y, method = "ar",main = "Spectral density function",
           xlab = "Frequency" ,ylab = "Spectrum")
}
