% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/script_final_reg_sar.R
\name{fit_newton}
\alias{fit_newton}
\title{Metodo de Newton com aproximacao numerica do gradiente e hessiana}
\usage{
fit_newton(
  lvero,
  gr,
  hess,
  formula,
  dados,
  A,
  init = NULL,
  verbose = FALSE,
  tol = .Machine$double.eps^(1/2),
  ...
)
}
\arguments{
\item{lvero}{Negativo da função de log-verossimilhança de interesse}

\item{gr}{Função gradiente analítica ou numérica}

\item{hess}{Função da matriz hessiana analítica ou numérica}

\item{formula}{Formula do R estilo y ~ b0 +b1X1 + ... + bnXn. Se não quiser
colocar covariáveis, usar y ~ 1.}

\item{dados}{Base de dados (n x k) com as covariáveis de interesse}

\item{A}{Matrix (n x n) de adjacências do modelo.}

\item{init}{Vetor de chutes nomeado, se não quiser o método implementado}

\item{verbose}{Exibe informações durante a maximização}

\item{tol}{Tolerância de convergência das estimativas.}

\item{...}{Passagem de argumentos extras de funções internas.}
}
\value{
lista com dados da otimização
}
\description{
Função genérica para obter estimativas de ll usando o método de
Newton que tem como base o vetor gradiente e a matrix hessiana. Nessa
implementação, as derivadas são calculadas via aproximação numérica
com apio do pacote numDeriv onde temos uma implementação do método de
Richardson.
}
