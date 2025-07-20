library(betareg)
library(mgcv)

# --- Funções auxiliares ---
train.beta <- function(x, y) {
  if (is.null(colnames(x))) {
    colnames(x) <- paste0("x", seq_len(ncol(x)))
  }
  df <- data.frame(y = y, x)
  formula <- as.formula(paste("y ~", paste(colnames(x), collapse = "+")))
  betareg(formula, data = df, link = "logit")
}

predict.beta <- function(model, newx) {
  if (is.null(colnames(newx))) {
    colnames(newx) <- names(model$coefficients$mean)[-1]
  }
  newdata <- as.data.frame(newx)
  predict(model, newdata = newdata, type = "response")
}

# Treina um GAM para modelar os resíduos absolutos
mad.train.gam <- function(x, res) {
  df <- data.frame(res = res, x)
  formula <- as.formula(paste("res ~", paste0("s(", colnames(x), ")", collapse = " + ")))
  gam(formula, data = df, family = gaussian())
}

# Faz predição dos resíduos absolutos estimados
mad.predict.gam <- function(model, newx) {
  predict(model, newdata = as.data.frame(newx), type = "response")
}

split.conformal.beta <- function(x, y, x0, 
                                 train.fun, predict.fun,
                                 mad.train.fun,
                                 mad.predict.fun,
                                 alpha = 0.1, rho = 0.5,
                                 link = "logit",
                                 split = NULL, seed = NULL, verbose = FALSE) {
  x <- as.matrix(x)
  y <- as.numeric(y)
  x0 <- as.matrix(x0)
  colnames(x) <- c("x1", "x2")
  colnames(x0) <- c("x1", "x2")
  n <- nrow(x)
  
  if (!all(y > 0 & y < 1)) stop("Os valores de y devem estar no intervalo (0,1)")
  
  if (!is.null(split)) {
    i1 <- split
  } else {
    if (!is.null(seed)) set.seed(seed)
    i1 <- sample(1:n, floor(n * rho))
  }
  i2 <- setdiff(1:n, i1)
  
  model <- train.fun(x[i1, , drop = FALSE], y[i1])
  
  pred_train <- predict.fun(model, x)
  pred_test <- predict.fun(model, x0)
  
  g <- switch(link,
              "logit" = qlogis,
              "probit" = qnorm,
              "cloglog" = function(mu) log(-log(1 - mu)),
              stop("Link não implementado"))
  
  g_inv <- switch(link,
                  "logit" = plogis,
                  "probit" = pnorm,
                  "cloglog" = function(eta) 1 - exp(-exp(eta)))
  
  # Calcular resíduos no espaço transformado
  eta_y <- g(y[i2])
  phi <- coef(model)['(phi)']
  mu_star <- digamma(pred_train[i2] * phi) - digamma((1 - pred_train[i2]) * phi)
  res <- abs(eta_y - mu_star)
  
  # === Escalonamento opcional com MAD ===
  if (!is.null(mad.train.fun) && !is.null(mad.predict.fun)) {
    # Residuals da amostra de treino para treinar o GAM
    eta_y_train <- g(y[i1])
    mu_star_train <- digamma(pred_train[i1] * phi) - digamma((1 - pred_train[i1]) * phi)
    res_train <- abs(eta_y_train - mu_star_train)
    
    # Ajustar o modelo de MAD
    mad_model <- mad.train.fun(x[i1, , drop = FALSE], res_train)
    mad_i2 <- mad.predict.fun(mad_model, x[i2, , drop = FALSE])
    mad_x0 <- mad.predict.fun(mad_model, x0)
    
    # Corrigir casos com previsão MAD negativa
    mad_i2 <- pmax(mad_i2, 1e-6)
    mad_x0 <- pmax(mad_x0, 1e-6)
    
    # Escalonar os resíduos
    res <- res / mad_i2
  } else {
    mad_x0 <- rep(1, nrow(x0))
  }
  q_alpha <- quantile(res, probs = 1 - alpha, type = 1)
  
  eta_pred <- digamma(pred_test * phi) - digamma((1 - pred_test) * phi)
  lo_eta <- eta_pred - q_alpha * mad_x0
  up_eta <- eta_pred + q_alpha * mad_x0
  
  lo <- g_inv(lo_eta)
  up <- g_inv(up_eta)
  
  return(list(
    pred = pred_test,
    lo = lo,
    up = up,
    fit = pred_train,
    i1 = i1,
    i2 = i2
  ))
}

# --- Simulação Monte Carlo com avaliação no conjunto de teste ---
B <- 100
n <- 1000
alpha <- 0.05
phi <- 50
beta <- c(-4, 0.5, 1)

gera_dados <- function(n) {
  x1 <- runif(n)
  x2 <- runif(n)
  X <- cbind(1, x1, x2)
  eta <- X %*% beta
  mu <- plogis(eta)
  y <- rbeta(n, mu * phi, (1 - mu) * phi)
  data.frame(x1 = x1, x2 = x2, y = y)
}

cobertura_total <- 0
amplitude_total <- 0
n_total <- 0

set.seed(123)

for (b in 1:B) {
  df <- gera_dados(n)
  x <- as.matrix(df[, c("x1", "x2")])
  y <- df$y
  
  res <- split.conformal.beta(x, y, x0 = x,
                              train.fun = train.beta,
                              predict.fun = predict.beta,
                              mad.train.fun = mad.train.gam,
                              mad.predict.fun = mad.predict.gam,
                              alpha = alpha,
                              rho = 0.7,
                              link = "logit",
                              seed = b)
  
  i2 <- res$i2
  y_val <- y[i2]
  lo <- res$lo[i2]
  up <- res$up[i2]
  
  # Cobertura
  coberto <- (y_val >= lo) & (y_val <= up)
  cobertura_total <- cobertura_total + sum(coberto)
  
  # Amplitude média
  amplitude_total <- amplitude_total + sum(up - lo)
  
  # Total de pontos de validação
  n_total <- n_total + length(i2)
}

# Resultados
cobertura_empirica <- cobertura_total / n_total
amplitude_media <- amplitude_total / n_total

cat(sprintf("Cobertura empírica: %.2f%% (esperado: %.0f%%)\n", 
            100 * cobertura_empirica, 100 * (1 - alpha)))
cat(sprintf("Amplitude média dos intervalos: %.4f\n", amplitude_media))
