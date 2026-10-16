# Bloco J (parte 1) — simulador com verdade conhecida.
#
# População finita em dois estágios: estratos > UPAs > pessoas. Cada UPA tem
# uma variável de desenho Z, NÃO observada pelo analista, que entra na medida
# de tamanho (PPS) e no desfecho. Por isso o desenho é informativo dado X
# quando kappa > 0. Versão compacta do simulador usado nos experimentos (que
# tinha um 3º estágio, domicílio > morador); aqui o 2º estágio sorteia pessoas
# diretamente.

#' Gera uma população finita sintética
#'
#' Modelo: `eta = g(X) + gamma Z (1 + delta X2) + u + alfa_h + e`, com
#' `corr(X1, Z) = rho`, `g(x) = 1,5 sen(1,5 x1) + x2^2 + x3 1(x1 > 0)`.
#' Medida de tamanho da UPA: `N_k exp(kappa Z_k)`.
#'
#' @param H Nº de estratos.
#' @param K_por_estrato UPAs por estrato na população.
#' @param media_pessoas_upa Média (Poisson) de pessoas por UPA.
#' @param p Nº de preditores (>= 3).
#' @param rho,gamma,delta,sd_u,sd_e,sd_estrato Parâmetros do modelo.
#' @param kappa Informatividade do 1º estágio (0 = não informativo).
#' @param desfecho `"continuo"` ou `"binario"`.
#' @param prevalencia Prevalência populacional do binário.
#' @param semente Semente.
#' @return Objeto `populacao_rf` (lista com tabelas de UPAs e pessoas).
#' @export
simular_populacao <- function(H = 8L, K_por_estrato = 200L, media_pessoas_upa = 60,
                              p = 6L, rho = 0.3, gamma = 1, delta = 0.5, sd_u = 0.3,
                              sd_e = 1, sd_estrato = 0.5, kappa = 0.7,
                              desfecho = c("continuo", "binario"), prevalencia = 0.1,
                              semente = NULL) {
  desfecho <- match.arg(desfecho)
  checar(p >= 3, "p deve ser >= 3.")
  com_semente(semente, {
    K <- H * K_por_estrato
    upa <- data.frame(upa = seq_len(K), estrato = rep(seq_len(H), each = K_por_estrato),
                      Z = stats::rnorm(K), u = stats::rnorm(K, 0, sd_u),
                      N = 1L + stats::rpois(K, media_pessoas_upa - 1))
    upa$mos <- upa$N * exp(kappa * upa$Z)
    alfa <- stats::rnorm(H, 0, sd_estrato)
    k <- rep(upa$upa, upa$N)
    n <- length(k)
    X <- matrix(stats::rnorm(n * p), n, p, dimnames = list(NULL, paste0("x", seq_len(p))))
    X[, 1] <- rho * upa$Z[k] + sqrt(1 - rho^2) * X[, 1]
    X[, 3] <- as.numeric(stats::runif(n) < 0.4)
    pessoas <- data.frame(upa = k, estrato = upa$estrato[k], X)
    eta <- g_simulador(X) + gamma * upa$Z[k] * (1 + delta * X[, 2]) + upa$u[k] +
      alfa[pessoas$estrato] + stats::rnorm(n, 0, sd_e)
    limiar <- if (desfecho == "binario") stats::quantile(eta, 1 - prevalencia, names = FALSE) else NA
    pessoas$y <- if (desfecho == "binario") as.numeric(eta > limiar) else eta
  })
  structure(list(upa = upa, pessoas = pessoas, alfa = alfa, limiar = limiar,
                 par = list(H = H, K_por_estrato = K_por_estrato, p = p, rho = rho,
                            gamma = gamma, delta = delta, sd_e = sd_e, kappa = kappa,
                            desfecho = desfecho)),
            class = "populacao_rf")
}

g_simulador <- function(X) 1.5 * sin(1.5 * X[, 1]) + X[, 2]^2 + X[, 3] * (X[, 1] > 0)

#' @export
print.populacao_rf <- function(x, ...) {
  cat("<populacao_rf> N =", nrow(x$pessoas), "| UPAs =", nrow(x$upa), "| estratos =",
      x$par$H, "| kappa =", x$par$kappa, "| desfecho:", x$par$desfecho, "\n")
  invisible(x)
}

#' Nomes dos preditores do simulador (x1..xp e estrato)
#' @param pop Objeto `populacao_rf`.
#' @export
preditores_simulador <- function(pop) c(paste0("x", seq_len(pop$par$p)), "estrato")

#' Verdade na população finita: E_pop\[Y | X = x, estrato\]
#'
#' Mistura exata sobre as UPAs realizadas do estrato, com peso
#' `N_k * densidade(x1 | Z_k)`. No binário, integra o ruído individual
#' (`P(eta > limiar)`). Ignora apenas a flutuação finita do ruído individual.
#'
#' @param pop Objeto `populacao_rf`.
#' @param novos data.frame com x1..xp e estrato.
#' @param bloco Linhas processadas por vez.
#' @return Vetor com a verdade em cada linha.
#' @export
verdade_populacao <- function(pop, novos, bloco = 500L) {
  pr <- pop$par
  X <- as.matrix(novos[, paste0("x", seq_len(pr$p))])
  s <- sqrt(1 - pr$rho^2)
  out <- numeric(nrow(novos))
  for (h in unique(novos$estrato)) {
    up <- pop$upa[pop$upa$estrato == h, ]
    ii <- which(novos$estrato == h)
    for (ini in seq(1L, length(ii), by = bloco)) {
      j <- ii[ini:min(length(ii), ini + bloco - 1L)]
      dens <- stats::dnorm(outer(X[j, 1], pr$rho * up$Z, "-") / s)
      Wt <- sweep(dens, 2, up$N, `*`)
      base <- g_simulador(X[j, , drop = FALSE]) + pop$alfa[h]
      efeito <- outer(1 + pr$delta * X[j, 2], pr$gamma * up$Z) +
        matrix(up$u, length(j), nrow(up), byrow = TRUE)
      if (pr$desfecho == "binario") {
        q <- stats::pnorm((base + efeito - pop$limiar) / pr$sd_e)
        out[j] <- rowSums(Wt * q) / rowSums(Wt)
      } else {
        out[j] <- base + rowSums(Wt * efeito) / rowSums(Wt)
      }
    }
  }
  out
}

#' Probabilidades de inclusão PPS com teto 1 (iterativo)
#' @noRd
inclusao_pps <- function(mos, a) {
  pik <- a * mos / sum(mos)
  repeat {
    acima <- pik >= 1
    if (!any(pik > 1)) break
    pik[acima] <- 1
    resto <- !acima
    pik[resto] <- (a - sum(acima)) * mos[resto] / sum(mos[resto])
  }
  pmin(pik, 1)
}

#' PPS sistemático com ordem aleatória
#' @noRd
pps_sistematico <- function(pik) {
  ord <- sample.int(length(pik))
  cs <- c(0, cumsum(pik[ord]))
  pts <- stats::runif(1) + 0:(round(sum(pik)) - 1)
  ord[findInterval(pts, cs)]
}

#' Sorteia uma amostra em dois estágios (PPS de UPAs, AAS de pessoas)
#'
#' @param pop Objeto `populacao_rf`.
#' @param a_por_estrato UPAs por estrato.
#' @param m_por_upa Pessoas por UPA (todas, se a UPA tiver menos).
#' @param semente Semente.
#' @return Lista com `X`, `y`, `desenho` (objeto `desenho_rf` com fpc) e
#'   `linhas` (índices na população).
#' @export
sortear_amostra <- function(pop, a_por_estrato = 20L, m_por_upa = 12L, semente = NULL) {
  com_semente(semente, {
    up <- pop$upa
    sel <- integer(0); pi1 <- numeric(0)
    for (h in seq_len(pop$par$H)) {
      ids <- which(up$estrato == h)
      pik <- inclusao_pps(up$mos[ids], a_por_estrato)
      esc <- pps_sistematico(pik)
      sel <- c(sel, ids[esc]); pi1 <- c(pi1, pik[esc])
    }
    ini <- c(0L, cumsum(up$N))[sel]
    linhas <- integer(0); w <- numeric(0)
    for (i in seq_along(sel)) {
      Nk <- up$N[sel[i]]; m <- min(m_por_upa, Nk)
      linhas <- c(linhas, ini[i] + sample.int(Nk, m))
      w <- c(w, rep(Nk / (m * pi1[i]), m))
    }
  })
  pes <- pop$pessoas[linhas, ]
  f1 <- a_por_estrato / pop$par$K_por_estrato
  list(X = pes[, preditores_simulador(pop)], y = pes$y, linhas = linhas,
       desenho = desenho_rf(peso = w, estrato = pes$estrato, upa = pes$upa,
                            fpc = rep(f1, length(w))))
}

#' Amostra aleatória simples de pessoas (referência sem desenho)
#' @param pop Objeto `populacao_rf`.
#' @param n Tamanho.
#' @param semente Semente.
#' @return Lista como em [sortear_amostra()].
#' @export
amostra_aas <- function(pop, n, semente = NULL) {
  linhas <- com_semente(semente, sample.int(nrow(pop$pessoas), n))
  pes <- pop$pessoas[linhas, ]
  list(X = pes[, preditores_simulador(pop)], y = pes$y, linhas = linhas,
       desenho = desenho_rf(peso = rep(nrow(pop$pessoas) / n, n)))
}
