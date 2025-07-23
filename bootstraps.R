rm(list = ls())

library(future)      # motor de paralelismo moderno
library(doFuture)    # integra foreach + future
library(foreach)
library(betareg)
library(pracma)
library(zipfR)

set.seed(497532)

## Detecta e registra o número de núcleos (exemplo usando SLURM ou fallback)
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = future::availableCores()))

plan(multicore, workers = n_cores)  # se Windows, trocar para multisession
registerDoFuture()
options(future.rng.onMisuse = "ignore")  # silenciar warning sobre sementes

#########################
## Funções Auxiliares  ##
#########################

# Função para construir y_modificado dada uma amostra de resíduos
y_modificado <- function(r_tb, mu_estimated_star, v_estimated) {
  exp(mu_estimated_star + r_tb * sqrt(v_estimated)) /
    (1 + exp(mu_estimated_star + r_tb * sqrt(v_estimated)))
}

erro_bootstrap_t <- function(y_calc, mu_estimated_star, v_estimated){
  (log(y_calc/(1-y_calc)) - mu_estimated_star) / (sqrt(v_estimated))
}

T_calc <- function(R_boot_t, n_calc, sigma_star){
  (R_boot_t - mean(R_boot_t)) / (sqrt(1 + n_calc^(-1)) * sigma_star)
}
# Função para calcular mu a partir de beta e x (x é um vetor)
mu_invert <- function(beta_invert, x_invert){
  Xtemp <- cbind(1, x_invert)
  eta   <- Xtemp %*% beta_invert
  mu    <- exp(eta)/(1 + exp(eta))
  return(mu)
}

# Função para calcular mu_star
mu_star_calc <- function(mu_estimated, phi_estimated){
  digamma(mu_estimated * phi_estimated) - 
    digamma((1 - mu_estimated) * phi_estimated)
}

# Função para calcular v
v_calc <- function(mu_estimated, phi_estimated ) {
  trigamma(mu_estimated * phi_estimated) + 
    trigamma((1 - mu_estimated) * phi_estimated)
}

# Erro de predição
pred_error <- function(y_aplus, mustar_plus, v_plus) {
  (log(y_aplus/(1 - y_aplus)) - mustar_plus ) / sqrt(v_plus)
}

# Esperança cúbica (para aceleração)
Esp_lt3 <- function(mu_hat_plus, phi_estimated){
  phi_estimated^3 * (psigamma(mu_hat_plus * phi_estimated , deriv = 2) - 
                       psigamma((1 - mu_hat_plus) * phi_estimated, deriv = 2) )
}

# Variância LT
var_lt <- function(mu_hat_plus, phi_estimated){
  phi_estimated^2 * (digamma(mu_hat_plus * phi_estimated) + 
                       digamma((1 - mu_hat_plus) * phi_estimated) ) 
}

# Função auxiliar para cálculo do alfa til (usada no BCa)
alfa_til_inside <- function(alpha, v0_hat, a_hat, quant = TRUE) {
  if (quant == TRUE) {
    termo <- v0_hat + (
      (v0_hat + qnorm(alpha/2)) /
        (1 - a_hat * (v0_hat + qnorm(alpha/2)))
    )
  } else {
    termo <- v0_hat + (
      (v0_hat + qnorm(1 - alpha/2)) /
        (1 - a_hat * (v0_hat + qnorm(1 - alpha/2)))
    )
  }
  return(termo)
}

# Função para calcular os limites dos intervalos pelo método Percentile
calc_limites_percentile <- function(alfa1, alfa2 , R_ab_sorted, 
                                    mu_chapeu_star_plus, v_chapeu_plus, y_obs) {
  n1 <- length(y_obs)
  y_lower     <- numeric(n1)
  y_upper     <- numeric(n1)
  in_interval <- logical(n1)
  
  for(i in 1:n1) {
    lim_inf_idx <- ceiling(B * alfa1)
    lim_sup_idx <- ceiling(B * alfa2)
    lim_inf_idx <- max(1, min(lim_inf_idx, B))
    lim_sup_idx <- max(1, min(lim_sup_idx, B))
    
    Rinf <- R_ab_sorted[i, lim_inf_idx]
    Rsup <- R_ab_sorted[i, lim_sup_idx]
    
    y_lower[i] <- exp(mu_chapeu_star_plus[i] + Rinf * sqrt(v_chapeu_plus[i])) / 
      (1 + exp(mu_chapeu_star_plus[i] + Rinf * sqrt(v_chapeu_plus[i])))
    y_upper[i] <- exp(mu_chapeu_star_plus[i] + Rsup * sqrt(v_chapeu_plus[i])) / 
      (1 + exp(mu_chapeu_star_plus[i] + Rsup * sqrt(v_chapeu_plus[i])))
    
    in_interval[i] <- (y_obs[i] >= y_lower[i] && y_obs[i] <= y_upper[i])
  }
  
  return(list(
    y_lower     = y_lower,
    y_upper     = y_upper,
    in_interval = in_interval
  ))
}

# Função para calcular os limites dos intervalos BCa
calc_limites_BCA <- function(alfa1, alfa2, R_ab_sorted, 
                             mu_chapeu_star_plus, v_chapeu_plus, y_obs) {
  n1 <- length(y_obs)
  y_lower     <- numeric(n1)
  y_upper     <- numeric(n1)
  in_interval <- logical(n1)
  
  for (i in 1:n1) {
    idx_in <- ceiling(B * alfa1[i])
    idx_out <- ceiling(B * alfa2[i])
    idx_in <- max(1, min(idx_in, B))
    idx_out <- max(1, min(idx_out, B))
    
    Rinf <- R_ab_sorted[i, idx_in]
    Rsup <- R_ab_sorted[i, idx_out]
    
    val_inf <- mu_chapeu_star_plus[i] + Rinf * sqrt(v_chapeu_plus[i])
    y_lower[i] <- exp(val_inf) / (1 + exp(val_inf))
    val_sup <- mu_chapeu_star_plus[i] + Rsup * sqrt(v_chapeu_plus[i])
    y_upper[i] <- exp(val_sup) / (1 + exp(val_sup))
    
    in_interval[i] <- (y_obs[i] >= y_lower[i] && y_obs[i] <= y_upper[i])
  }
  
  return(list(
    y_lower     = y_lower,
    y_upper     = y_upper,
    in_interval = in_interval
  ))
}


calc_limites_bootstrap_t <- function(alfa1, alfa2 , R_mean, sigma_n, mu_chapeu_star_plus, 
                                     n_calc, v_chapeu_plus, y_obs, T_sorted) {
  n1 <- length(y_obs)
  y_lower     <- numeric(n1)
  y_upper     <- numeric(n1)
  in_interval <- logical(n1)
  
  for(i in 1:n1) {
    lim_inf_idx <- R_mean + ceiling(B * alfa1) * sigma_n * sqrt(1 + n_calc^(-1))
    lim_sup_idx <- R_mean + ceiling(B * alfa2) * sigma_n * sqrt(1 + n_calc^(-1))
    lim_inf_idx <- max(1, min(lim_inf_idx, B))
    lim_sup_idx <- max(1, min(lim_sup_idx, B))
    
    Tinf <- T_sorted[i, lim_inf_idx]
    Tsup <- T_sorted[i, lim_sup_idx]
    
    y_lower[i] <- exp(mu_chapeu_star_plus[i] + Tinf * sqrt(v_chapeu_plus[i])) / 
      (1 + exp(mu_chapeu_star_plus[i] + Tinf * sqrt(v_chapeu_plus[i])))
    y_upper[i] <- exp(mu_chapeu_star_plus[i] + Tsup * sqrt(v_chapeu_plus[i])) / 
      (1 + exp(mu_chapeu_star_plus[i] + Tsup * sqrt(v_chapeu_plus[i])))
    
    in_interval[i] <- (y_obs[i] >= y_lower[i] && y_obs[i] <= y_upper[i])
  }
  
  return(list(
    y_lower     = y_lower,
    y_upper     = y_upper,
    in_interval = in_interval
  ))
}






###############################
## Estudo de Simulação MC    ##
###############################

MC <- 10    
B  <- 500
sample_sizes <- c(50)
phi_values   <- c(50, 150, 400)
alpha        <- 0.05

results_list <- list()

start <- Sys.time()

for (n in sample_sizes) {
  n_total <- 50 + n
  n1 <- n_total - n
  
  for (phi in phi_values) {
    
    resultados_mc <- foreach(mc = 1:MC,
                             .combine = 'rbind',
                             .packages = c("betareg", "pracma", "zipfR"),
                             .options.future = list(seed = TRUE)) %dopar% {
                               
                               # Gera uma amostra de covariáveis
                               x1_inicial   <- runif(n_total)
                               x2_inicial   <- runif(n_total)
                               
                               # Conjunto de treinamento
                               X            <- cbind(1, x1_inicial[1:n], x2_inicial[1:n])
                               
                               beta         <- c(-1.5, 1.5, 1)
                               
                               mu           <- exp(X %*% beta) / (1 + exp(X %*% beta))
                               y_train      <- rbeta(n, mu * phi, (1 - mu) * phi)
                               
                               # Conjunto de teste
                               X_plus_vec1       <- x1_inicial[(n+1):n_total]
                               X_plus_vec2       <- x2_inicial[(n+1):n_total]
                               X_plus            <- cbind(1, X_plus_vec1, X_plus_vec2)
                               
                               mu_plus_true <- exp(X_plus %*% beta) / (1 + exp(X_plus %*% beta))
                               y_test       <- rbeta(n1, mu_plus_true * phi, (1 - mu_plus_true) * phi)
                               
                               ## Ajuste do modelo beta com os dados de treinamento
                               data_train    <- data.frame(y = y_train, x1 = X[,2], x2 = X[,3])
                               modelo        <- betareg(y ~ x1 + x2, data = data_train)
                               
                               # Extração dos parâmetros estimados
                               phi_modelo_estimado    <- coef(modelo)["(phi)"]
                               betas_modelo_estimado  <- coef(modelo)[1:length(beta)]
                               mu_hat_modelo          <- modelo$fitted.values
                               
                               # Quantidades auxiliares do treinamento
                               mu_hat_star       <- mu_star_calc(mu_hat_modelo, phi_modelo_estimado)
                               v_hat_modelo      <- v_calc(mu_hat_modelo, phi_modelo_estimado)
                               
                               # Previsão para o conjunto de teste usando os parâmetros estimados
                               mu_hat_plus       <- mu_invert(betas_modelo_estimado, X_plus[,2:3])
                               mu_hat_star_plus  <- mu_star_calc(mu_hat_plus, phi_modelo_estimado)
                               v_hat_plus        <- v_calc(mu_hat_plus, phi_modelo_estimado)
                               
                               # Resíduos do modelo ajustado
                               resid_quantile        <- residuals(modelo, type = "quantile")
                               resid_sweighted2      <- residuals(modelo, type = "sweighted2")
                               
                               # Sigmas estimados para o bootstrap-t
                               sigma_estimated_n_q   <- sqrt(n^(-1) * sum( (resid_quantile - mean(resid_quantile))^2))
                               sigma_estimated_n_sw2 <- sqrt(n^(-1) * sum( (resid_sweighted2 - mean(resid_sweighted2))^2))
                               
                               # Matrizes para armazenar resultados bootstrap para o conjunto de teste
                               R_ab_quantile         <- matrix(NA, nrow = n1, ncol = B)
                               R_ab_sweighted2       <- matrix(NA, nrow = n1, ncol = B)
                               T_ab_quantile         <- matrix(NA, nrow = n1, ncol = B)
                               T_ab_sw2              <- matrix(NA, nrow = n1, ncol = B)
                               
                               # Para os intervalos percentile (opcionalmente pode-se armazenar y_boot, mas aqui usamos R_ab)
                               for (b in 1:B) {
                                 # Reamostragem dos resíduos de treinamento
                                 r_tb_quantile       <- sample(resid_quantile, size = n, replace = TRUE)
                                 r_tb_sweighted2     <- sample(resid_sweighted2, size = n, replace = TRUE)
                                 
                                 # Construção de pseudo-respostas para treinamento
                                 y_b_quantile        <- y_modificado(r_tb_quantile, mu_hat_star, v_hat_modelo)
                                 y_b_sweighted2      <- y_modificado(r_tb_sweighted2, mu_hat_star, v_hat_modelo)
                                 
                                 # Ajuste do modelo bootstrap
                                 data_b_quantile     <- data.frame(y = y_b_quantile, x1 = X[,2], x2 = X[,3])
                                 data_b_sweighted2   <- data.frame(y = y_b_sweighted2, x1 = X[,2], x2 = X[,3])
                                 modelo_b_quantile   <- betareg(y ~ x1 + x2, data = data_b_quantile)
                                 modelo_b_sweighted2 <- betareg(y ~ x1 + x2, data = data_b_sweighted2)
                                 
                                 # Extração dos parâmetros bootstrap
                                 phi_hat_boot_quantile       <- coef(modelo_b_quantile)["(phi)"]
                                 beta_hat_quantile_boot      <- coef(modelo_b_quantile)[1:length(beta)]
                                 mu_hat_plus_boot_quantile   <- mu_invert(beta_hat_quantile_boot, X_plus[,2:3])
                                 Rstar_quantile              <- mu_star_calc(mu_hat_plus_boot_quantile, phi_hat_boot_quantile)
                                 v_plus_boot_quantile        <- v_calc(mu_hat_plus_boot_quantile, phi_hat_boot_quantile)
                                 
                                 phi_hat_boot_sweighted2     <- coef(modelo_b_sweighted2)["(phi)"]
                                 beta_hat_boot_sweighted2    <- coef(modelo_b_sweighted2)[1:length(beta)]
                                 mu_hat_plus_boot_sweighted2 <- mu_invert(beta_hat_boot_sweighted2, X_plus[,2:3])
                                 Rstar_sweighted2            <- mu_star_calc(mu_hat_plus_boot_sweighted2, phi_hat_boot_sweighted2)
                                 v_plus_boot_sweighted2      <- v_calc(mu_hat_plus_boot_sweighted2, phi_hat_boot_sweighted2)
                                 
                                 # Reamostragem dos resíduos para o conjunto de teste
                                 r_aplus_b_quantile          <- sample(resid_quantile, size = n1, replace = TRUE)
                                 r_aplus_b_sweighted2        <- sample(resid_sweighted2, size = n1, replace = TRUE)
                                 
                                 # Constrói as pseudo-respostas para o conjunto de teste (bootstrap)
                                 y_aplus_b_quantile          <- y_modificado(r_aplus_b_quantile, Rstar_quantile, v_plus_boot_quantile)
                                 y_aplus_b_sweighted2        <- y_modificado(r_aplus_b_sweighted2, Rstar_sweighted2, v_plus_boot_sweighted2)
                                 
                                 # Calculo dos erros bootstrap 
                                 R_boot_t_quantile           <- erro_bootstrap_t(y_aplus_b_quantile, Rstar_quantile,v_plus_boot_quantile)
                                 R_boot_t_sw2                <- erro_bootstrap_t(y_aplus_b_sweighted2, Rstar_sweighted2,v_plus_boot_sweighted2 )
                                 
                                 R_boot_t_quantile_mean      <- sum(R_boot_t_quantile) / length(R_boot_t_quantile)
                                 R_boot_t_sw2_mean           <- sum(R_boot_t_sw2) / length(R_boot_t_sw2)
                                 
                                 # Cálculo sigma estrela
                                 sigma_star_quantile         <- sqrt(sum((R_boot_t_quantile - R_boot_t_quantile_mean)^2)  / length(R_boot_t_quantile))
                                 sigma_star_sw2              <- sqrt(sum((R_boot_t_sw2 - R_boot_t_sw2_mean)^2)  / length(R_boot_t_sw2))
                                 
                                 
                                 # Cálculo dos erros de predição bootstrap
                                 R_ab_quantile[,b]    <- pred_error(y_aplus_b_quantile, mu_hat_star_plus, v_hat_plus)
                                 R_ab_sweighted2[,b]  <- pred_error(y_aplus_b_sweighted2, mu_hat_star_plus, v_hat_plus)
                                 
                                 T_ab_quantile[,b]    <- T_calc(R_boot_t_quantile, n, sigma_star_quantile)
                                 T_ab_sw2[,b]         <- T_calc(R_boot_t_sw2, n, sigma_star_sw2)
                               } # fim do loop bootstrap
                               
                               # Ordena os valores bootstrap por linha (para cada observação do teste)
                               R_ab_quantile_sorted        <- t(apply(R_ab_quantile, 1, sort))
                               R_ab_sweighted2_sorted      <- t(apply(R_ab_sweighted2, 1, sort))
                               
                               T_ab_quantile_sorted        <- t(apply(T_ab_quantile, 1, sort))
                               T_ab_sw2_sorted             <- t(apply(T_ab_sw2, 1, sort))
                               
                               ## Intervalos Percentile
                               y_quantile_percentile   <- calc_limites_percentile(alfa1 = 0.025, alfa2 = 0.975,
                                                                                  R_ab_quantile_sorted, mu_hat_star_plus, v_hat_plus, y_obs = y_test)
                               y_sweighted2_percentile <- calc_limites_percentile(alfa1 = 0.025, alfa2 = 0.975,
                                                                                  R_ab_sweighted2_sorted, mu_hat_star_plus, v_hat_plus, y_obs = y_test)
                               
                               ## Intervalos BCa
                               
                               # Cálculo dos quantis dos resíduos de treinamento
                               R_m_quantile <- median(resid_quantile)
                               R_m_sweighted2 <- median(resid_sweighted2)
                               num_R_ab_quantile <- apply(R_ab_quantile_sorted, 1, function(row) sum(row < R_m_quantile))
                               num_R_ab_sweighted2 <- apply(R_ab_sweighted2_sorted, 1, function(row) sum(row < R_m_sweighted2))
                               
                               v0_hat_quantile <- qnorm((num_R_ab_quantile + 0.5)/(B+1))
                               v0_hat_sweighted2 <- qnorm((num_R_ab_sweighted2 + 0.5)/ (B+1))
                               
                               Esp_lt3_quantile <- Esp_lt3(mu_hat_plus, phi_modelo_estimado)
                               Esp_lt3_sweighted2 <- Esp_lt3(mu_hat_plus, phi_modelo_estimado)
                               var_lt_quantile <- var_lt(mu_hat_plus, phi_modelo_estimado)
                               var_lt_sweighted2 <- var_lt(mu_hat_plus, phi_modelo_estimado)
                               
                               var_lt3_2_quantile <- var_lt_quantile^(3/2)
                               var_lt3_2_sweighted2 <- var_lt_sweighted2^(3/2)
                               
                               a_correcao <- (n/n1)^(-1/2)
                               a_hat_quantile <- a_correcao * 1/6 * (Esp_lt3_quantile/var_lt3_2_quantile)
                               a_hat_sweighted2 <- a_correcao * 1/6 * (Esp_lt3_sweighted2/var_lt3_2_sweighted2)
                               
                               alfa_til_1_quantile <- pnorm(alfa_til_inside(alpha = 0.05, v0_hat = v0_hat_quantile,
                                                                            a_hat = a_hat_quantile, quant = TRUE))
                               alfa_til_2_quantile <- pnorm(alfa_til_inside(alpha = 0.05, v0_hat = v0_hat_quantile,
                                                                            a_hat = a_hat_quantile, quant = FALSE))
                               
                               alfa_til_1_sweighted2 <- pnorm(alfa_til_inside(alpha = 0.05, v0_hat = v0_hat_sweighted2,
                                                                              a_hat = a_hat_sweighted2, quant = TRUE))
                               alfa_til_2_sweighted2 <- pnorm(alfa_til_inside(alpha = 0.05, v0_hat = v0_hat_sweighted2,
                                                                              a_hat = a_hat_sweighted2, quant = FALSE))
                               
                               y_quantile_BCa <- calc_limites_BCA(alfa_til_1_quantile, alfa_til_2_quantile,
                                                                  R_ab_quantile_sorted, mu_hat_star_plus, v_hat_plus, y_obs = y_test)
                               y_sweighted2_BCa <- calc_limites_BCA(alfa_til_1_sweighted2, alfa_til_2_sweighted2,
                                                                    R_ab_sweighted2_sorted, mu_hat_star_plus, v_hat_plus, y_obs = y_test)
                               
                               
                               # Calcula os intervalos bootstrap-t
                               y_quantile_boot_t    <- calc_limites_bootstrap_t(alfa1 = 0.025, alfa2 = 0.975,
                                                                                R_mean = mean(resid_quantile), 
                                                                                sigma_n = sigma_estimated_n_q, 
                                                                                mu_chapeu_star_plus = mu_hat_star_plus, 
                                                                                n_calc = n, 
                                                                                v_chapeu_plus = v_hat_plus, y_obs = y_test,
                                                                                T_sorted = T_ab_quantile_sorted)
                               
                               y_sw2_boot_t         <- calc_limites_bootstrap_t(alfa1 = 0.025, alfa2 = 0.975,
                                                                                R_mean = mean(resid_sweighted2), 
                                                                                sigma_n = sigma_estimated_n_sw2, 
                                                                                mu_chapeu_star_plus = mu_hat_star_plus, 
                                                                                n_calc = n, 
                                                                                v_chapeu_plus = v_hat_plus, y_obs = y_test,
                                                                                T_sorted = T_ab_sw2_sorted)
                               
                              
                               
                               ## Armazenamento dos resultados: cobertura e amplitude média dos intervalos
                               
                               c(
                                 mean(y_quantile_percentile$in_interval, na.rm = TRUE),
                                 mean(y_quantile_percentile$y_upper - y_quantile_percentile$y_lower, na.rm = TRUE),
                                 mean(y_quantile_BCa$in_interval, na.rm = TRUE),
                                 mean(y_quantile_BCa$y_upper - y_quantile_BCa$y_lower, na.rm = TRUE),
                                 mean(y_quantile_boot_t$in_interval, na.rm = TRUE),
                                 mean(y_quantile_boot_t$y_upper - y_quantile_boot_t$y_lower, na.rm = TRUE),
                                 mean(y_sweighted2_percentile$in_interval, na.rm = TRUE),
                                 mean(y_sweighted2_percentile$y_upper - y_sweighted2_percentile$y_lower, na.rm = TRUE),
                                 mean(y_sweighted2_BCa$in_interval, na.rm = TRUE),
                                 mean(y_sweighted2_BCa$y_upper - y_sweighted2_BCa$y_lower, na.rm = TRUE),
                                 mean(y_sw2_boot_t$in_interval, na.rm = TRUE),
                                 mean(y_sw2_boot_t$y_upper - y_sw2_boot_t$y_lower, na.rm = TRUE)
                               )
                               
                               
                             }
    results_list[[paste("n", n, "phi", phi, sep = "_")]] <- list(
      quantile_percentile_coverage     = mean(resultados_mc[,1]),
      quantile_percentile_amplitude    = mean(resultados_mc[,2]),
      quantile_BCa_coverage            = mean(resultados_mc[,3]),
      quantile_BCa_amplitude           = mean(resultados_mc[,4]),
      quantile_boot_t_coverage         = mean(resultados_mc[,5]),
      quantile_boot_t_amplitude        = mean(resultados_mc[,6]),
      sweighted2_percentile_coverage   = mean(resultados_mc[,7]),
      sweighted2_percentile_amplitude  = mean(resultados_mc[,8]),
      sweighted2_BCa_coverage          = mean(resultados_mc[,9]),
      sweighted2_BCa_amplitude         = mean(resultados_mc[,10]),
      sweighted2_boot_t_coverage       = mean(resultados_mc[,11]),
      sweighted2_boot_t_amplitude      = mean(resultados_mc[,12])
      
    )
    
    cat("Concluído: n =", n, "phi =", phi, "\n")
  }
}

dur <- Sys.time() - start
print(dur)

###############################
## Data-frame final           ##
###############################
results_df <- data.frame(
  n                                  = integer(),
  phi                                = numeric(),
  quantile_percentile_coverage       = numeric(),
  quantile_percentile_amplitude      = numeric(),
  quantile_BCa_coverage              = numeric(),
  quantile_BCa_amplitude             = numeric(),
  quantile_boot_t_coverage           = numeric(),
  quantile_boot_t_amplitude          = numeric(),
  sweighted2_percentile_coverage     = numeric(),
  sweighted2_percentile_amplitude    = numeric(),
  sweighted2_BCa_coverage            = numeric(),
  sweighted2_BCa_amplitude           = numeric(),
  sweighted2_boot_t_coverage         = numeric(),
  sweighted2_boot_t_amplitude        = numeric(),
  stringsAsFactors = FALSE
)

for (key in names(results_list)) {
  split_key <- strsplit(key, "_")[[1]]
  n_val <- as.numeric(split_key[2])
  phi_val <- as.numeric(split_key[4])
  
  res <- results_list[[key]]
  
  row_df <- data.frame(
    n    = n_val,
    phi  = phi_val,
    quantile_percentile_coverage    = res$quantile_percentile_coverage,
    quantile_percentile_amplitude   = res$quantile_percentile_amplitude,
    quantile_BCa_coverage           = res$quantile_BCa_coverage,
    quantile_BCa_amplitude          = res$quantile_BCa_amplitude,
    quantile_boot_t_coverage        = res$quantile_boot_t_coverage,
    quantile_boot_t_amplitude       = res$quantile_boot_t_amplitude,
    sweighted2_percentile_coverage  = res$sweighted2_percentile_coverage,
    sweighted2_percentile_amplitude = res$sweighted2_percentile_amplitude,
    sweighted2_BCa_coverage         = res$sweighted2_BCa_coverage,
    sweighted2_BCa_amplitude        = res$sweighted2_BCa_amplitude,
    sweighted2_boot_t_coverage      = res$sweighted2_boot_t_coverage,
    sweighted2_boot_t_amplitude     = res$sweighted2_boot_t_amplitude,
    stringsAsFactors = FALSE
  )
  
  results_df <- rbind(results_df, row_df)
}

print(results_df)