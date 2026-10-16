# Bloco H — incerteza de desenho da previsão f(x), por florestas-réplica.
#
# Evidência (simulador, passo 6): a variância entre árvores superestima a
# variância da floresta ~10x; meias-amostras ("little bags") foram ruidosas e
# infladas; florestas-réplica Rao-Wu ficaram em VR = 0,73-0,85 (subestimam
# ~20%). Para a média de Hájek o mesmo Rao-Wu deu VR ~ 0,95, então a falta vem
# da floresta. Os intervalos cobrem E[f(x)], NÃO a função verdadeira: o viés de
# suavização da floresta derrubou a cobertura de m(x) para 0,44-0,67.

#' Incerteza de desenho da previsão por florestas-réplica
#'
#' Refaz o procedimento inteiro (mesma estratégia, motor e controle) em
#' `n_replicas` réplicas do desenho e mede a variância entre as florestas,
#' descontando o ruído de Monte Carlo de cada floresta-réplica e somando o da
#' floresta principal:
#' \deqn{\hat V(x) = \mathrm{Var}_r(F_r(x)) - \overline{s^2_r(x)}/B_r + s^2(x)/B.}
#' Com `fonte = "orgao"`, usa os pesos replicados do desenho e a fórmula de
#' variância replicada informada (escala e rscales), sem desconto de Monte Carlo
#' por réplica além do termo s2/B_r.
#'
#' @param fit Objeto `floresta_desenho`.
#' @param novos data.frame dos pontos x de interesse.
#' @param n_replicas Nº de florestas-réplica (ignorado se `fonte = "orgao"`).
#' @param arvores_por_replica Árvores em cada floresta-réplica.
#' @param fonte `"raowu"` (réplicas físicas de UPAs) ou `"orgao"` (pesos replicados).
#' @param nivel Nível do intervalo normal.
#' @param nucleos Processos paralelos (fork; use 1 no Windows).
#' @param semente Semente.
#' @return data.frame com `previsao`, `var`, `ep`, `inf`, `sup` e `var_negativa`
#'   (o estimador de desenho deu < 0 e foi truncado em 0). Atributo `"nota"`.
#' @export
incerteza_replicas <- function(fit, novos, n_replicas = 50L, arvores_por_replica = 50L,
                               fonte = c("raowu", "orgao"), nivel = 0.95, nucleos = 1L,
                               semente = NULL) {
  fonte <- match.arg(fonte)
  checar(inherits(fit, "floresta_desenho"), "fit deve ser floresta_desenho.")
  novos <- como_preditores(novos)
  d <- fit$desenho
  P <- predict(fit, novos, tipo = "arvores")
  f <- rowMeans(P); s2_mc <- apply(P, 1, stats::var) / ncol(P)
  ctl <- fit$controle; ctl$n_arvores <- as.integer(arvores_por_replica)
  R <- if (fonte == "orgao") {
    checar(!is.null(d$pesos_rep), "fonte = 'orgao' exige pesos replicados no desenho.")
    ncol(d$pesos_rep)
  } else n_replicas
  checar(R >= 2, "S\u00E3o necess\u00E1rias pelo menos 2 r\u00E9plicas.")
  seeds <- com_semente(semente, sample.int(.Machine$integer.max, R))

  uma <- function(r) {
    set.seed(seeds[r])
    rp <- if (fonte == "orgao") replica_orgao(d, r) else replica_raowu(d)
    dr <- subconjunto_desenho(d, rp$idx, peso = rp$peso, upa = rp$upa)
    dr$pesos_rep <- NULL
    est <- fit$estrategia
    if (fonte == "orgao" && est$nome == "replicas_orgao") est <- reamostragem_wb()
    # estrato com n_h = 2 vira 1 UPA na réplica Rao-Wu: colapsa para manter a estratégia
    if (!is.null(est$requer(dr))) dr <- suppressWarnings(colapsar_estratos_unitarios(dr))
    fr <- floresta_desenho(fit$X[rp$idx, , drop = FALSE], fit$y[rp$idx], dr, est,
                           fit$motor, ctl)
    Pr <- predict(fr, novos, tipo = "arvores")
    cbind(rowMeans(Pr), apply(Pr, 1, stats::var))
  }
  res <- if (nucleos > 1L && .Platform$OS.type != "windows")
    parallel::mclapply(seq_len(R), uma, mc.cores = nucleos) else lapply(seq_len(R), uma)
  Fr <- vapply(res, function(m) m[, 1], numeric(nrow(novos)))
  Wr <- vapply(res, function(m) m[, 2], numeric(nrow(novos)))
  Fr <- matrix(Fr, nrow(novos)); Wr <- matrix(Wr, nrow(novos))
  ruido <- rowMeans(Wr) / arvores_por_replica
  v_des <- if (fonte == "orgao") {
    rsc <- d$rep_rscales %||% rep(1, R)
    d$rep_escala * rowSums(sweep((Fr - f)^2, 2, rsc, `*`)) - ruido * d$rep_escala * sum(rsc)
  } else apply(Fr, 1, stats::var) - ruido
  v <- pmax(v_des, 0) + s2_mc
  z <- stats::qnorm(1 - (1 - nivel) / 2)
  out <- data.frame(previsao = f, var = v, ep = sqrt(v), inf = f - z * sqrt(v),
                    sup = f + z * sqrt(v), var_negativa = v_des < 0)
  attr(out, "nota") <- paste(
    "Incerteza de desenho da previs\u00E3o (alvo: E[f(x)]), n\u00E3o intervalo para a fun\u00E7\u00E3o",
    "verdadeira. No simulador, Rao-Wu subestimou a vari\u00E2ncia em ~20%.")
  out
}
