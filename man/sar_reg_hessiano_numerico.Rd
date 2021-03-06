% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/script_final_reg_sar.R
\name{sar_reg_hessiano_numerico}
\alias{sar_reg_hessiano_numerico}
\title{Matriz hessiana com derivadas numericas modelo SAR}
\usage{
sar_reg_hessiano_numerico(
  param,
  func,
  method = "Richardson",
  method.args = list(),
  ...
)
}
\arguments{
\item{param}{vetor de parâmetros onde func será avaliada por dirivadas numéricas aproximadas;}

\item{func}{Função implementada da log-verossimilhança ou outra qualquer;}

\item{method}{Método de derivação. Veja \code{\link[numDeriv]{hessian}};}

\item{method.args}{argumentos extras passados para a função.}

\item{...}{Passagem de argumentos extras de funções internas.}
}
\description{
Determina numericamente a matriz hessiana de uma função.
Para mais detalhes ver  \code{\link[numDeriv]{hessian}};
}
