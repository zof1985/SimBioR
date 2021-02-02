\name{residuals_analysis}

\alias{residuals_analysis}

\title{Winter's Residuals Analysis}

\usage{
residuals_analysis(
  x,
  sampling_rate = 1,
  frequencies   = seq(0.01, 0.49, 0.01) * sampling_rate,
  nodes         = 2,
  keep_node     = nodes,
  filter_fun    = NULL,
  metric        = c("rmse", "sse", "mse"),
  return_ggplot = T
  )
}

\description{
Generate an object of class "ResidualsAnalysis" containing all the data used to obtain an optimal cut-off frequency for the provided signal according to the used filter.

The residuals analysis approach \insertCite{Winter2009}{SimBioR} is used to determine the optimal cut-off frequency for a given filter applied to a signal. This is a geometrical technique mainly used to obtain the cut-off frequency of low-pass filters and generally it is not completely automatic.
Brefly, the technique consists in calculating, for each tested frequency, the residuals between the original and the filtered signal using the selected metric. Afterwards, the user must decide when the residuals vs frequency curve becomes linear. Then, the intercept of the line best fitting the linear region of this curve is used as threshold for defining the optimal cutoff frequency.

To make the process completely automatic, a grid-search approach \insertCite{Lerman1980}{SimBioR} is used. This approach consists in fitting the residuals/frequency curve by means of all the possible combinations of n + 1 segments (where n = number of nodes). The combination of nodes minimizing the residuals is then selected. Finally, the optimal frequency is defined as the frequency corresponding to the node selected by the user (keep_node).

This procedure is computationally intensive, therefore increasing either the number of tested frequencies or the number of nodes will result in a significantly longer computation time.
}

\arguments{

\item{x}{A numeric vector.}

\item{sampling_rate}{A numeric value representing the sampling frequency (in Hz) used to collect the signal x.}

\item{frequencies}{The list of frequencies tested by the function. By default, 200 frequencies are tested within the (0, Nyquist frequency) range.}

\item{nodes}{The number of nodes used to define the segments to be fitted during the grid search step. A value > 1 has to be used.}

\item{keep_node}{The node to be selected as resulting optimal cut-off frequency.}

\item{filter_fun}{The filter function to be used. By default (NULL), a 4th order, phase-corrected, Butterworth low-pass filter is used. Alternatively, a function can be passed accepting two arguments: "f" the tested frequency, and "x" the signal. The function should return a numeric vector with the same length of "x".}

\item{metric}{The metric to be used for calculating the residuals. The following measures are currently supported: "rmse" (Root Mean Square Error), "mse" (Mean square Error), "sse" (Sum of Square Error).}

\item{return_ggplot}{If TRUE (default) a ggplot object is generated describing the outcomes of the analysis. Otherwise, only the optimal frequency is returned.}

}

\value{

If return_ggplot is TRUE (default), a list with two keys:

 \item{opt}{the optimal frequency detected.}
 \item{plot}{a ggplot object ready to be rendered summarizing the outcomes of the performed residuals analysis.}

If return_ggplot is FALSE:

\item{opt}{the optimal frequency detected.}

}


\references{

\insertAllCited{}

}

\examples{


#### generate a noisy signal ####


# sampling frequency
fs <- 20

# signal frequency
f <- 2.4

# sampling time
t <- 5

# samples
x <- seq(0, t, 1 / fs)

# example signal with added white noise
y <- sin(2 * pi * f / t * x) + 0.25 * rnorm(length(x))



#### perform the analysis ####


# obtain just the optimal frequency without plot
opt <- residuals_analysis(y, sampling_rate = fs, return_ggplot = FALSE)


# obtain the optimal cut-off frequency using a custom frequency range
freqs <- seq(0.1, fs / 2, 0.1)
ra <- residuals_analysis(y, sampling_rate = fs, frequencies = freqs)
ra$opt  # optimal frequency
ra$plot # plot of the analysis
}
