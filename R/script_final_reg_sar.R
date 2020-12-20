## --------------------------------------------------------------
## MI602: Métodos Computacionais em Estatística, Atividade III
## Prof.: Guilherme Ludwig
## Aluno:José Evandeilton lopes (RA264072)
## Data: 14/12/2020
## Objetivo: Programar o modelo SAR dado na atividade III
## --------------------------------------------------------------

## -------------------------------------------------------------
## (1) Implementação da log-vero, gradiente e hessiana numéricos
## -------------------------------------------------------------

#' Log-verossimilhanca com covariaveis modelo SAR
#' @description Determina a log-verossimilhança com base
#' em uma realização de uma normal multivariada.
#' @param param vetor nomeado de parâmetros a serem estimados
#' @param formula Formula do R estilo y ~ b0 +b1X1 + ... + bnXn. Se não quiser
#' colocar covariáveis, usar y ~ 1.
#' @param dados Base de dados (n x k) com as covariáveis de interesse
#' @param A Matrix (n x n) de adjacências do modelo.
#' @export
sar_reg_lvero_dmvnorm <- function(param, formula, dados, A){
  # @importFrom mvtnorm dmvnorm
  # @importFrom stats model.frame model.response model.matrix

  mfx <- model.frame(formula, dados)
  Y   <- model.response(mfx)
  X   <- model.matrix(formula, mfx)
  K   <- c(1:(ncol(X)+2))
  n   <- nrow(X)
  if(length(param) < ncol(X)+2) {
    stop("info: vetor de parametros precisa ser tamanho ncol(X)+2\n")
  }
  pars <- param[1:(ncol(X)+2)]
  names(pars)[K[1:(length(K)-2)]] <- colnames(X)
  names(pars)[K[length(K)-1]] <- "sigma"
  names(pars)[K[length(K)]] <- "rho"

  I  <- diag(nrow(X))
  S  <- solve(I-pars["rho"]*A)
  mu <- X%*%param[1:ncol(X)]
  si <- as.matrix(pars["sigma"]^2*S%*%S)
  ll <- mvtnorm::dmvnorm(Y, mean = mu, sigma = si, log = TRUE)
  return(-ll)
}

#' Vetor gradiente com derivadas numericas modelo SAR
#' @description Determina numericamente ovetor gradiente de uma função.
#' Para mais detalhes ver \code{\link[numDeriv]{grad}};
#' @param param vetor de parâmetros onde func será avaliada por dirivadas numéricas aproximadas;
#' @param func Função implementada da log-verossimilhança ou outra qualquer;
#' @param method Método de derivação. Veja \code{\link[numDeriv]{grad}};
#' @param method.args argumentos extras passados para a função.
#' @param ... Passagem de argumentos extras de funções internas.
#' @importFrom numDeriv grad
#' @export
sar_reg_gradiente_numerico <- function(param, func, method="Richardson", method.args = list(), ...) {
  numDeriv::grad(func = func, x = param, method = method, method.args = method.args, ...)
}

#' Matriz hessiana com derivadas numericas modelo SAR
#' @description Determina numericamente a matriz hessiana de uma função.
#' Para mais detalhes ver  \code{\link[numDeriv]{hessian}};
#' @param param vetor de parâmetros onde func será avaliada por dirivadas numéricas aproximadas;
#' @param func Função implementada da log-verossimilhança ou outra qualquer;
#' @param method Método de derivação. Veja \code{\link[numDeriv]{hessian}};
#' @param method.args argumentos extras passados para a função.
#' @param ... Passagem de argumentos extras de funções internas.
#' @importFrom numDeriv hessian
#' @export
sar_reg_hessiano_numerico <- function(param, func, method="Richardson", method.args = list(), ...) {
  numDeriv::hessian(func = func, x = param, method = method, method.args = method.args, ...)
}

## -------------------------------------------------------------
## (2) Implementação da log-vero, gradiente e hessiana analíticos
## -------------------------------------------------------------

#' Log-verossimilhanca analitica com covariaveis modelo SAR
#' @param param Vetor nomeado de parâmetros a serem estimados
#' @param formula Formula do R estilo y ~ b0 +b1X1 + ... + bnXn. Se não quiser
#' colocar covariáveis, usar y ~ 1.
#' @param dados Base de dados (n x k) com as covariáveis de interesse
#' @param A Matrix (n x n) de adjacências do modelo.
#' @importFrom stats model.frame model.response model.matrix
#' @export
sar_reg_lvero_analitico <- function(param, formula, dados, A){
  # f(y) =(0.5)*log(det((I-rho*A)'*(I-rho*A))) -(n/2)*log(2*pi)-(n/2)*log(sigma^2)-(1/(2*sigma^2))*((Y-X*beta)')*((I-rho*A)'*(I-rho*A))*(Y-X*beta)
  mfx <- model.frame(formula, dados)
  Y   <- model.response(mfx)
  X   <- model.matrix(formula, mfx)
  n   <- length(Y); I   <- diag(n); K   <- ncol(X)
  beta <- param[1:K]
  sigma <- param[K+1]
  rho <- param[K+2]
  T_0 <- (rho * A)
  T_1 <- (t(I) - T_0)
  T_2 <- (I - T_0)
  t_3 <- (Y - (X)%*%(beta))
  t_4 <- (sigma ** 2)
  t_5 <- (2 * t_4)
  t_6 <- (T_1)%*%((T_2)%*%(t_3))
  t_7 <- (1 / t_5)

  ll <- ((((0.5 * log(det((T_1)%*%(T_2)))) - ((n * log((2 * pi))) / 2)) - ((n * log(t_4)) / 2)) - (t(t_3)%*%(t_6) / t_5))

  return(-sum(ll))
}

#' Vetor gradiente com derivadas analiticas modelo SAR
#' @description Determina analiticamente a matriz hessiana do modelo SAR.
#' @param param Vetor nomeado de parâmetros a serem estimados
#' @param formula Formula do R estilo y ~ b0 +b1X1 + ... + bnXn. Se não quiser
#' colocar covariáveis, usar y ~ 1.
#' @param dados Base de dados (n x k) com as covariáveis de interesse
#' @param A Matrix (n x n) de adjacências do modelo.
#' @param ... Passagem de argumentos extras de funções internas.
#' @importFrom stats model.frame model.response model.matrix
#' @export
sar_reg_gradiente_analitico <- function(param, formula, dados, A, ...){
  #d/dbeta = (0.5)*log(det((I-rho*A)'*(I-rho*A))) -(n/2)*log(2*pi)-(n/2)*log(sigma^2)-(1/(2*sigma^2))*((Y-X*beta)')*((I-rho*A)'*(I-rho*A))*(Y-X*beta) = 1/(2*sigma.^2)*X'*(I'-rho*A)*(I-rho*A)*(Y-X*beta)+1/(2*sigma.^2)*X'*(I'-rho*A)*(I+(-rho*A)')*(Y-X*beta)

  #d/dsigma (0.5)*log(det((I-rho*A)'*(I-rho*A))) -(n/2)*log(2*pi)-(n/2)*log(sigma^2)-(1/(2*sigma^2))*((Y-X*beta)')*((I-rho*A)'*(I-rho*A))*(Y-X*beta) = ((Y-X*beta)'*(I'-rho*A)*(I-rho*A)*(Y-X*beta))/sigma.^3-n/sigma

  # d/drho (0.5)*log(det((I-rho*A)'*(I-rho*A))) -(n/2)*log(2*pi)-(n/2)*log(sigma^2)-(1/(2*sigma^2))*((Y-X*beta)')*((I-rho*A)'*(I-rho*A))*(Y-X*beta) = -(0.5*tr(A*(I-rho*A)*inv((I'-rho*A)*(I-rho*A)))+0.5*tr(A*inv((I'-rho*A)*(I-rho*A))*(I'-rho*A))-(((Y-X*beta)'*A*(I-rho*A)*(Y-X*beta))/(2*sigma.^2)+((Y-X*beta)'*(I'-rho*A)*A*(Y-X*beta))/(2*sigma.^2)))

  mfx <- model.frame(formula, dados)
  Y   <- model.response(mfx)
  X   <- model.matrix(formula, mfx)
  n   <- length(Y)
  I   <- diag(n)
  K   <- ncol(X)
  beta <- param[1:K]
  sigma <- param[K+1]
  rho <- param[K+2]

  Ub <- {
    T_0 = (rho * A)
    T_1 = (t(I) - T_0)
    T_2 = (I - T_0)
    t_3 = (Y - (X)%*%(beta))
    t_4 = (sigma ** 2)
    t_5 = (2 * t_4)
    t_6 = (T_1)%*%((T_2)%*%(t_3))
    t_7 = (1 / t_5)
    ((t_7 * t(X)%*%(t_6)) + (t_7 * t(X)%*%((T_1)%*%(((I + -t(T_0)))%*%(t_3)))))
  }

  Us <- {
    T_0 = (rho * A)
    T_1 = (t(I) - T_0)
    T_2 = (I - T_0)
    t_3 = (Y - (X)%*%(beta))
    t_4 = (sigma ** 2)
    t_5 = t(t_3)%*%((T_1)%*%((T_2)%*%(t_3)))
    ((t_5 / (sigma ** 3)) - (n / sigma))
  }

  Ur <- {
    T_0 = (rho * A)
    T_1 = (t(I) - T_0)
    T_2 = (I - T_0)
    t_3 = (Y - (X)%*%(beta))
    t_4 = (sigma ** 2)
    T_5 = (T_1)%*%(T_2)
    T_6 = solve(T_5)
    t_7 = (T_2)%*%(t_3)
    t_8 = (2 * t_4)
    -(((0.5 * traco(((A)%*%(T_2))%*%(T_6))) + (0.5 * traco(((A)%*%(T_6))%*%(T_1)))) - ((t(t_3)%*%((A)%*%(t_7)) / t_8) + (t(t_3)%*%((T_1)%*%((A)%*%(t_3))) / t_8)))
  }
  U <- c(t(Ub), Us, Ur)
  return(U)
}

#' Matriz hessiana com derivadas analiticas modelo SAR
#' @description Determina a matriz hessiana analítica do modelo SAR.
#' @param param Vetor nomeado de parâmetros a serem estimados
#' @param formula Formula do R estilo y ~ b0 +b1X1 + ... + bnXn. Se não quiser
#' colocar covariáveis, usar y ~ 1.
#' @param dados Base de dados (n x k) com as covariáveis de interesse
#' @param A Matrix (n x n) de adjacências do modelo.
#' @param ... Passagem de argumentos extras de funções internas.
#' @importFrom stats model.frame model.response model.matrix
#' @export
sar_reg_hessiano_analitico <- function(param, formula, dados, A, ...){
  #H11 = d/(dbeta dbeta) 1/(2*sigma.^2)*X'*(I'-rho*A)*(I-rho*A)*(Y-X*beta)+1/(2*sigma.^2)*X'*(I'-rho*A)*(I+(-rho*A)')*(Y-X*beta) = -1/sigma.^2*X'*(I'-rho*A)*(I-rho*A)*X

  #H12 = d/(dbeta dsigma)  1/(2*sigma.^2)*X'*(I'-rho*A)*(I-rho*A)*(Y-X*beta)+1/(2*sigma.^2)*X'*(I'-rho*A)*(I+(-rho*A)')*(Y-X*beta) = -(2*sigma)/(sigma.^2).^2*X'*(I'-rho*A)*(I-rho*A)*(Y-X*beta)

  #H13 = d/(dbeta drho) 1/(2*sigma.^2)*X'*(I'-rho*A)*(I-rho*A)*(Y-X*beta)+1/(2*sigma.^2)*X'*(I'-rho*A)*(I+(-rho*A)')*(Y-X*beta) = -(1/sigma.^2*X'*A*(I-rho*A)*(Y-X*beta)+1/sigma.^2*X'*(I'-rho*A)*A*(Y-X*beta))

  # H21 = H12

  # H22 = d/(dsigma dsigma) ((Y-X*beta)'*(I'-rho*A)*(I-rho*A)*(Y-X*beta))/sigma.^3-n/sigma = n/sigma.^2-(sigma.^2*3*(Y-X*beta)'*(I'-rho*A)*(I-rho*A)*(Y-X*beta))/(sigma.^3).^2

  # H23 = d/(dsigma drho) ((Y-X*beta)'*(I'-rho*A)*(I-rho*A)*(Y-X*beta))/sigma.^3-n/sigma = -(((Y-X*beta)'*A*(I-rho*A)*(Y-X*beta))/sigma.^3+((Y-X*beta)'*(I'-rho*A)*A*(Y-X*beta))/sigma.^3)

  # H31 = H13
  # H32 = H23

  # d/(drho drho) -(0.5*tr(A*(I-rho*A)*inv((I'-rho*A)*(I-rho*A)))+0.5*tr(A*inv((I'-rho*A)*(I-rho*A))*(I'-rho*A))-(((Y-X*beta)'*A*(I-rho*A)*(Y-X*beta))/(2*sigma.^2)+((Y-X*beta)'*(I'-rho*A)*A*(Y-X*beta))/(2*sigma.^2))) = 0.5*tr(A*inv((I'-rho*A)*(I-rho*A))*A)-(0.5*tr(A*(I-rho*A)*inv((I'-rho*A)*(I-rho*A))*A*(I-rho*A)*inv((I'-rho*A)*(I-rho*A)))+0.5*tr(A*inv((I'-rho*A)*(I-rho*A))*A*(I-rho*A)*inv((I'-rho*A)*(I-rho*A))*(I'-rho*A)))-(0.5*tr(A*(I-rho*A)*inv((I'-rho*A)*(I-rho*A))*(I'-rho*A)*A*inv((I'-rho*A)*(I-rho*A)))+0.5*tr(A*inv((I'-rho*A)*(I-rho*A))*(I'-rho*A)*A*inv((I'-rho*A)*(I-rho*A))*(I'-rho*A)))+0.5*tr(A*A*inv((I'-rho*A)*(I-rho*A)))-((Y-X*beta)'*A*A*(Y-X*beta))/sigma.^2

  mfx <- model.frame(formula, dados)
  Y   <- model.response(mfx)
  X   <- model.matrix(formula, mfx)
  n   <- length(Y)
  I   <- diag(n)
  K   <- ncol(X)
  beta <- param[1:K]
  sigma <- param[K+1]
  rho <- param[K+2]

  H11 <- {
    T_0 = (rho * A)
    t_1 = (1 / (sigma ** 2))
    T_2 = (t(I) - T_0)
    T_3 = (I - T_0)
    -(t_1 * ((t(X)%*%(T_2))%*%(T_3))%*%(X))
  }

  H12 <- {
    T_0 = (rho * A)
    t_1 = (sigma ** 2)
    t_2 = t(X)%*%(((t(I) - T_0))%*%(((I - T_0))%*%((Y - (X)%*%(beta)))))
    -(((2 * sigma) / (t_1 ** 2)) * t_2)
  }

  H13 <- {
    T_0 = (rho * A)
    t_1 = (1 / (sigma ** 2))
    t_2 = (Y - (X)%*%(beta))
    t_3 = ((I - T_0))%*%(t_2)
    T_4 = (t(I) - T_0)
    -((t_1 * t(X)%*%((A)%*%(t_3))) + (t_1 * t(X)%*%((T_4)%*%((A)%*%(t_2)))))
  }

  H21 <- H12

  H22 <- {
    T_0 = (rho * A)
    t_1 = (Y - (X)%*%(beta))
    t_2 = t(t_1)%*%(((t(I) - T_0))%*%(((I - T_0))%*%(t_1)))
    t_3 = (sigma ** 3)
    ((n / (sigma ** 2)) - ((((sigma ** 2) * 3) * t_2) / (t_3 ** 2)))
  }

  H23 <- {
    T_0 = (rho * A)
    t_1 = (Y - (X)%*%(beta))
    t_2 = ((I - T_0))%*%(t_1)
    t_3 = (sigma ** 3)
    T_4 = (t(I) - T_0)
    -((t(t_1)%*%((A)%*%(t_2)) / t_3) + (t(t_1)%*%((T_4)%*%((A)%*%(t_1))) / t_3))
  }

  H31 <- H13
  H32 <- H23

  H33 <- {
    T_0 = (rho * A)
    T_1 = (I - T_0)
    T_2 = (t(I) - T_0)
    T_3 = solve((T_2)%*%(T_1))
    t_4 = (Y - (X)%*%(beta))
    t_5 = (sigma ** 2)
    t_6 = (2 * t_5)
    T_7 = (A)%*%(T_3)
    T_8 = ((A)%*%(T_1))%*%(T_3)
    T_9 = (T_7)%*%(T_2)
    t_10 = (A)%*%(t_4)
    (((((0.5 * traco((T_7)%*%(A))) - ((0.5 * traco((((T_8)%*%(A))%*%(T_1))%*%(T_3))) + (0.5 * traco(((((T_7)%*%(A))%*%(T_1))%*%(T_3))%*%(T_2))))) - ((0.5 * traco((((T_8)%*%(T_2))%*%(A))%*%(T_3))) + (0.5 * traco((((T_9)%*%(A))%*%(T_3))%*%(T_2))))) + (0.5 * traco(((A)%*%(A))%*%(T_3)))) - (t(t_4)%*%((A)%*%(t_10)) / t_5))
  }

  H    <- rbind(cbind(H11, H12, H13),
                c(H21, H22, H23),
                c(H31, H32, H33))
  rownames(H) <- colnames(H) <- names(param)
  return(H)
}

## -------------------------------------------------------------
## (3) Implementação do método de Newton
## -------------------------------------------------------------

#' Metodo de Newton com aproximacao numerica do gradiente e hessiana
#' @description Função genérica para obter estimativas de ll usando o método de
#' Newton que tem como base o vetor gradiente e a matrix hessiana. Nessa
#' implementação, as derivadas são calculadas via aproximação numérica
#' com apio do pacote numDeriv onde temos uma implementação do método de
#' Richardson.
#' @param lvero Negativo da função de log-verossimilhança de interesse
#' @param gr Função gradiente analítica ou numérica
#' @param hess Função da matriz hessiana analítica ou numérica
#' @param formula Formula do R estilo y ~ b0 +b1X1 + ... + bnXn. Se não quiser
#' colocar covariáveis, usar y ~ 1.
#' @param dados Base de dados (n x k) com as covariáveis de interesse
#' @param A Matrix (n x n) de adjacências do modelo.
#' @param init Vetor de chutes nomeado, se não quiser o método implementado
#' @param verbose Exibe informações durante a maximização
#' @param tol Tolerância de convergência das estimativas.
#' @param ... Passagem de argumentos extras de funções internas.
#' @return lista com dados da otimização
#' @importFrom bbmle parnames
#' @importFrom stats lm
#' @importFrom dplyr bind_rows
#' @export
fit_newton <- function(lvero, gr, hess, formula, dados, A, init = NULL, verbose = FALSE, tol = .Machine$double.eps^(1/2), ...){
  ## Chute inicial com estimativas de OLS para os betas e sigma
  re <- lm(formula, data = dados)
  if(is.null(init)){
    x0 <- c(coef(re), sigma = summary(re)$sigma, rho = 1/summary(re)$sigma^3)
  } else {
    x0 <- init
  }
  if(is.null(names(x0))){
    stop("INFO: init precisa ser um vetor nomeado")
  }
  parnames(lvero) <- names(x0)
  t_  <- 9; i   <- 1; mle <- l <- crit <- c()
  while(t_ > tol){
    H <- hess(param = x0, func = lvero, formula = formula, dados = dados, A = A, ...)
    U <- gr(param = x0, func = lvero, formula = formula, dados = dados, A = A, ...)
    x1 <- x0 - solve(H, U)
    x0 <- x1
    ml <- lvero(param = x0, formula = formula, dados = dados, A = A, ...)
    if(verbose){
      cat("log: iteracao ", i,
          #"mle: ", paste0(x0, collapse = ","),
          "logLik:", ml,
          "crit U^2:", t_, "\n")
    }
    t_ <- sum(U^2)
    mle[[i]] <- x0
    l[[i]] <- ml
    crit[[i]] <- t_
    i <- i + 1
  }
  pr <- data.frame(dplyr:::bind_rows(mle),
                   loglik = as.numeric(l),
                   crit = as.numeric(crit), check.names = F)

  return(list(coeficientes = x0,
              processamento = pr,
              ll = lvero,
              formula = formula, dados = dados, A = A))
}

#' Re-estimacao dos parametros com os otimos obtidos por Newton
#' @description Essa função serve para obter as estatísticas do modelo e
#' confirmar a convergência usando método de quasi-Newton implementado
#' na função optim. Se der erro no algoritmo BFGS, tentamos novamente
#' com o BFGS com restriçã (L-BFGS-B). O Objeto de saida é um modelo
#' completo onde se pode usar funções genéricas como coef(), summary(),
#' vcov(), obter intervalos de confiança e perfis de verossimilhança
#' @param fit Objeto obtido pela função fit_newton
#' @return Um objeto da classe 'mle2' herdada do pacote bbmle onde
#' podemos obter várias estatísticas
#' @importFrom bbmle mle2
#' @export
fit_quasi_newton <- function(fit){
  fit <- try(bbmle::mle2(minuslogl = fit$ll,
                  start = fit$coeficientes,
                  method="BFGS",
                  optimizer = "optim",
                  data = list(dados = fit$dados, A = fit$A, formula = fit$formula)))
  out <- if(class(fit) != "try-error"){
    fit
  } else {
    try(bbmle::mle2(minuslogl = fit$ll,
             start = abs(fit$coeficientes),
             method="L-BFGS-B",
             optimizer = "optim",
             lower = abs(x0)-0.5,
             upper = abs(x0)+0.5,
             data = list(dados = fit$dados, A = fit$A, formula = fit$formula)))
  }
  return(out)
}

## -----------------------------------------------------------
## (4) Preditos
## -----------------------------------------------------------

#' Determina os preditos do modelo ajustado
#' @description Como E(Y) = mu e mu é constante o vetor de valores esperados
#' é constante. Essa predict não é boa para o desenho do problema proposto.
#' @param fit Objeto obtido através da função fit_quasi_newton().
#' @param idcol Nome da variável que será a chave entre com o vetor de predições Y.
#' @return data.frame com duas colunas (id e pred)
#' @importFrom stats model.frame model.response model.matrix
#' @export
sar_reg_predict <- function(fit, idcol = "CODBAIRRO"){

  mfx <- model.frame(fit@data$formula, fit@data$dados)
  Y   <- model.response(mfx)
  X   <- model.matrix(fit@data$formula, mfx)
  A   <- fit@data$A
  I   <- diag(length(Y))

  coefs <- coef(fit)
  rho   <- coefs[length(coefs)]
  out <- as.numeric(X%*%coefs[1:(length(coefs)-2)])
  return(data.frame(id = as.numeric(fit@data$dados[,idcol]), pred = out))
}

#' Determina os preditos out-of-sampledo modelo ajustado
#' @description Essa é uma abordagem que pode ser útil quando temos apenas
#' uma constante como mu, mas se aplica também ao modelo de regressão
#' com covariáveis.
#' @param fit Objeto obtido através da função fit_quasi_newton().
#' @param idcol Nome da variável que será a chave entre com o vetor de predições Y.
#' @importFrom stats model.frame model.response model.matrix
#' @return data.frame com duas colunas (id e pred)
#' @export
sar_reg_predict_oot <- function(fit, idcol = "CODBAIRRO"){
  mfx <- model.frame(fit@data$formula, fit@data$dados)
  Y   <- model.response(mfx)
  X   <- model.matrix(fit@data$formula, mfx)
  A   <- fit@data$A
  I   <- diag(length(Y))
  coefs <- bbmle::coef(fit)
  beta  <- coefs[1:(length(coefs)-2)]
  sigma <- coefs[length(coefs)-1]
  rho   <- coefs[length(coefs)]
  mu    <- as.numeric(X%*%beta)
  S_    <- sigma^2*solve(I-rho*A)
  S     <- S_%*%S_
  out   <- numeric(length(Y))
  for(i in seq_along(Y)){
    out[i] <- mu[i] + S[i,-i]%*%solve(S[-i,-i])%*%(Y[-i]-mu[i])
  }
  return(data.frame(id = as.numeric(fit@data$dados[,idcol]), pred = out))
}

# Determina o traco de uma matrix m
traco <- function(m){
  sum(diag(m))
}

## -------------------------------------------------------------
## FIM
## -------------------------------------------------------------
