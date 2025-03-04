% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FSBLUP.R
\name{FSBLUP.stas.cal}
\alias{FSBLUP.stas.cal}
\title{statistical indicator}
\usage{
FSBLUP.stas.cal(x1, x2, type = c("cor", "auc", "rmse", "other"), fn = NULL)
}
\arguments{
\item{x1}{phenotype}

\item{x2}{predicted GEBV}

\item{type}{indicator type, e.g. "cor", "auc", "rmse", "other"}

\item{fn}{self custom function, valid only when "type = other"}
}
\value{
value
}
\description{
Note x1 and x2 must in the same order
}
