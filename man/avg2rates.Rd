% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AffineMortality.R
\name{avg2rates}
\alias{avg2rates}
\title{avg2rates}
\usage{
avg2rates(mu_bar)
}
\arguments{
\item{mu_bar}{the table of average mortality rates}
}
\value{
Returns a matrix of the same dimension as mu_bar with mortality rates
}
\description{
Converts a table with average mortality rates into mortality rates
}
\examples{
mu_bar_rates <- avg2rates(mu_bar)
}
