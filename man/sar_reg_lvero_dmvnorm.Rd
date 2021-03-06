% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/script_final_reg_sar.R
\name{sar_reg_lvero_dmvnorm}
\alias{sar_reg_lvero_dmvnorm}
\title{Log-verossimilhanca com covariaveis modelo SAR}
\usage{
sar_reg_lvero_dmvnorm(param, formula, dados, A)
}
\arguments{
\item{param}{vetor nomeado de parâmetros a serem estimados}

\item{formula}{Formula do R estilo y ~ b0 +b1X1 + ... + bnXn. Se não quiser
colocar covariáveis, usar y ~ 1.}

\item{dados}{Base de dados (n x k) com as covariáveis de interesse}

\item{A}{Matrix (n x n) de adjacências do modelo.}
}
\description{
Determina a log-verossimilhança com base
em uma realização de uma normal multivariada.
}
