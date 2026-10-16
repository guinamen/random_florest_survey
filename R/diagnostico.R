# Bloco D — diagnóstico do desenho e regras de decisão (portões).

#' Diagnóstico do desenho (e, opcionalmente, do desfecho)
#'
#' Mede o que, nos experimentos, mudou a configuração adequada: UPAs por
#' estrato, fração amostral de UPAs, variabilidade dos pesos, UPAs por célula
#' de calibração, ICC do desfecho e informatividade do desenho.
#'
#' O teste de informatividade é o de DuMouchel e Duncan (1983): compara, por
#' F, `y ~ X` com `y ~ X + w + X:w` (mínimos quadrados sem peso). Um efeito
#' significativo do peso dado X indica que a seleção carrega informação sobre y.
#' O p-valor ignora o conglomerado e tende a ser liberal com ICC alto; use-o
#' como triagem.
#'
#' @param desenho Objeto `desenho_rf`.
#' @param y Desfecho (opcional): habilita ICC e informatividade.
#' @param X Preditores (opcional): habilita a informatividade.
#' @param min_upas_celula Limiar de "poucas UPAs" por célula de calibração.
#' @return Objeto `diagnostico_rf`.
#' @export
diagnosticar <- function(desenho, y = NULL, X = NULL, min_upas_celula = 5L) {
  d <- desenho
  nh <- upas_por_estrato(d)
  w <- d$peso
  out <- list(
    n = d$n, n_estratos = length(nh), n_upas = length(unique(d$upa)),
    upas_por_estrato = c(min = min(nh), mediana = stats::median(nh), max = max(nh)),
    estratos_unitarios = sum(nh < 2L),
    f1_max = if (is.null(d$fpc)) NA_real_ else max(d$fpc),
    cv_peso = stats::sd(w) / mean(w),
    n_eff_kish = n_eff_kish(w), n_eff_upa = n_eff_upa(w, d$upa),
    pesos_rep = if (is.null(d$pesos_rep)) 0L else ncol(d$pesos_rep),
    celulas = NULL, icc = NA_real_, informatividade = NULL)
  if (!is.null(d$celula)) {
    uc <- tapply(d$upa, d$celula, function(v) length(unique(v)))
    out$celulas <- list(n = length(uc), min_upas = min(uc),
                        poucas = sum(uc < min_upas_celula), limiar = min_upas_celula)
  }
  if (!is.null(y)) {
    checar(length(y) == d$n, "y com tamanho errado.")
    out$icc <- icc_anova(y, d$upa)
    if (!is.null(X)) out$informatividade <- teste_informatividade(y, como_preditores(X), w)
  }
  structure(out, class = "diagnostico_rf")
}

#' Teste de DuMouchel-Duncan
#' @noRd
teste_informatividade <- function(y, X, w) {
  Z <- stats::model.matrix(~ ., data = X)[, -1, drop = FALSE]
  wn <- w / mean(w)
  m0 <- stats::lm(y ~ Z)
  m1 <- stats::lm(y ~ Z + wn + Z:wn)
  a <- stats::anova(m0, m1)
  list(F = a$F[2], gl = c(a$Df[2], a$Res.Df[2]), p = a$`Pr(>F)`[2],
       delta_r2 = summary(m1)$r.squared - summary(m0)$r.squared)
}

#' @export
print.diagnostico_rf <- function(x, ...) {
  cat("<diagnostico_rf>\n")
  cat(sprintf("  n = %d | estratos = %d | UPAs = %d\n", x$n, x$n_estratos, x$n_upas))
  cat(sprintf("  UPAs/estrato: m\u00EDn %d, mediana %.1f, m\u00E1x %d | estratos com 1 UPA: %d\n",
              x$upas_por_estrato[["min"]], x$upas_por_estrato[["mediana"]],
              x$upas_por_estrato[["max"]], x$estratos_unitarios))
  if (!is.na(x$f1_max)) cat(sprintf("  fra\u00E7\u00E3o de UPAs (f_h) m\u00E1x: %.4f\n", x$f1_max))
  cat(sprintf("  CV dos pesos: %.3f | n_eff (Kish): %.0f | n_eff em UPAs: %.0f\n",
              x$cv_peso, x$n_eff_kish, x$n_eff_upa))
  if (!is.null(x$celulas))
    cat(sprintf("  c\u00E9lulas de calibra\u00E7\u00E3o: %d | UPAs m\u00EDn: %d | com < %d UPAs: %d\n",
                x$celulas$n, x$celulas$min_upas, x$celulas$limiar, x$celulas$poucas))
  if (x$pesos_rep > 0) cat("  pesos replicados:", x$pesos_rep, "\n")
  if (!is.na(x$icc)) cat(sprintf("  ICC do desfecho (UPA): %.3f\n", x$icc))
  if (!is.null(x$informatividade))
    cat(sprintf("  informatividade (DuMouchel-Duncan): F = %.2f, p = %.3g, \u0394R\u00B2 = %.4f\n",
                x$informatividade$F, x$informatividade$p, x$informatividade$delta_r2))
  invisible(x)
}

#' Recomendações a partir do diagnóstico
#'
#' Transforma as medidas em ações explícitas. Cada linha traz a base da regra:
#' `[medido]` (resultado do simulador), `[inferência]` ou `[construção]`
#' (proposta sem teste).
#'
#' @param diag Objeto `diagnostico_rf`.
#' @param estimando Objeto [estimando()] (opcional).
#' @param limiar_icc,limiar_cv,limiar_f1,alfa Limiares das regras.
#' @return data.frame com `sinal`, `valor`, `acao`, `base`.
#' @export
recomendar <- function(diag, estimando = NULL, limiar_icc = 0.05, limiar_cv = 1,
                       limiar_f1 = 0.1, alfa = 0.05) {
  r <- list()
  add <- function(sinal, valor, acao, base)
    r[[length(r) + 1L]] <<- data.frame(sinal = sinal, valor = valor, acao = acao, base = base)
  if (diag$estratos_unitarios > 0)
    add("estratos com 1 UPA", as.character(diag$estratos_unitarios),
        "colapsar_estratos_unitarios() antes de Rao-Wu", "[requisito do m\u00E9todo]")
  if (diag$n_upas == diag$n && diag$pesos_rep > 0)
    add("UPAs n\u00E3o identificadas, h\u00E1 r\u00E9plicas", as.character(diag$pesos_rep),
        "usar reamostragem_replicas() e incerteza com fonte = 'orgao'", "[constru\u00E7\u00E3o]")
  if (!is.na(diag$f1_max) && diag$f1_max > limiar_f1)
    add("fra\u00E7\u00E3o de UPAs alta", sprintf("%.3f", diag$f1_max),
        "considerar fpc na incerteza; efeito pequeno no simulador", "[medido]")
  if (!is.null(diag$celulas) && diag$celulas$poucas > 0)
    add("c\u00E9lulas de calibra\u00E7\u00E3o com poucas UPAs", as.character(diag$celulas$poucas),
        "comparar peso calibrado e n\u00E3o calibrado (calibra\u00E7\u00E3o elevou o vi\u00E9s residual)", "[medido]")
  if (diag$cv_peso > limiar_cv)
    add("pesos muito vari\u00E1veis", sprintf("CV = %.2f", diag$cv_peso),
        "o WB concentra sorteios em poucas unidades; acompanhe n_eff", "[infer\u00EAncia]")
  if (!is.na(diag$icc) && diag$icc > limiar_icc)
    add("ICC relevante", sprintf("%.3f", diag$icc),
        "OOB por UPA para avaliar; estimativas por dom\u00EDnio com vari\u00E2ncia de desenho",
        "[medido]")
  if (!is.null(diag$informatividade)) {
    inf <- diag$informatividade$p < alfa
    add("informatividade", sprintf("p = %.3g", diag$informatividade$p),
        if (inf) "desenho informativo: pesos necess\u00E1rios para descrever a popula\u00E7\u00E3o"
        else "sem evid\u00EAncia de informatividade: pesos dispens\u00E1veis para ranking",
        "[constru\u00E7\u00E3o]")
  }
  if (!is.null(estimando)) {
    acao <- switch(estimando$tipo,
      ranking = "reamostragem_raowu_wb() + motor_ranger(); avaliar com auc_oob()",
      media_condicional = "reamostragem_raowu_wb() + motor_ranger(); avaliar com erro_oob()",
      dominio = "media_hajek() ou survey; n\u00E3o usar a floresta como estimador")
    add("estimando", estimando$tipo, acao, "[medido]")
  }
  if (!length(r)) return(data.frame(sinal = character(), valor = character(),
                                    acao = character(), base = character()))
  do.call(rbind, r)
}
