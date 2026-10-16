# Bloco E — estratégias de reamostragem.
#
# Contrato: toda estratégia é um objeto `estrategia_rf` cujo `gerar(d, B)`
# devolve
#   m   : matriz n x B (contagens inteiras ou multiplicadores de peso)
#   oob : matriz lógica n x B ("fora da bolsa" na árvore b)
# e cujo campo `saida` diz como o núcleo deve passar a árvore ao motor:
#   "contagens" -> a árvore recebe as contagens e NÃO usa o peso (o peso já
#                  entrou no sorteio). Serve a qualquer motor com "contagens".
#   "pesos"     -> a árvore recebe v = m * peso e precisa de corte ponderado.

#' Construtor genérico de estratégia de reamostragem
#'
#' Use para registrar um esquema novo (por exemplo, o bootstrap de
#' pseudo-população da SDRF) sem alterar o resto do pacote.
#'
#' @param nome Rótulo.
#' @param gerar Função `(d, B) -> list(m, oob)`; ver a seção "Contrato".
#' @param saida `"contagens"` ou `"pesos"`.
#' @param usa_peso A estratégia incorpora o peso amostral?
#' @param oob_por Nível do OOB: `"upa"`, `"unidade"`, `"replica"` ou `"nenhum"`.
#' @param requer Função `(d) -> NULL` ou mensagem de erro (pré-condições).
#' @section Contrato:
#' `m` tem dimensão n x B; para `saida = "contagens"` deve ser inteira e não
#' negativa. `oob[i, b]` só pode ser `TRUE` se `m[i, b] == 0`.
#' @return Objeto `estrategia_rf`.
#' @export
estrategia_reamostragem <- function(nome, gerar, saida = c("contagens", "pesos"),
                                    usa_peso, oob_por = c("upa", "unidade", "replica", "nenhum"),
                                    requer = function(d) NULL) {
  structure(list(nome = nome, gerar = gerar, saida = match.arg(saida),
                 usa_peso = usa_peso, oob_por = match.arg(oob_por), requer = requer),
            class = "estrategia_rf")
}

#' @export
print.estrategia_rf <- function(x, ...) {
  cat("<estrategia_rf>", x$nome, "| sa\u00EDda:", x$saida, "| usa peso:", x$usa_peso,
      "| OOB por:", x$oob_por, "\n")
  invisible(x)
}

#' Executa a estratégia e valida o contrato
#' @noRd
gerar_reamostras <- function(estrategia, d, B) {
  msg <- estrategia$requer(d)
  checar(is.null(msg), msg)
  r <- estrategia$gerar(d, B)
  checar(identical(dim(r$m), c(d$n, as.integer(B))), "Estrat\u00E9gia '", estrategia$nome,
         "' devolveu m com dimens\u00E3o errada.")
  checar(identical(dim(r$oob), dim(r$m)), "oob com dimens\u00E3o errada.")
  checar(!any(r$oob & r$m > 0), "Contrato violado: unidade OOB com m > 0.")
  if (estrategia$saida == "contagens")
    checar(all(r$m >= 0 & r$m == round(r$m)), "Contagens devem ser inteiras e >= 0.")
  r
}

requer_duas_upas <- function(d) {
  nh <- upas_por_estrato(d)
  if (any(nh < 2L))
    return(paste0("Rao-Wu exige >= 2 UPAs por estrato; ", sum(nh < 2L),
                  " estrato(s) com 1 UPA. Use colapsar_estratos_unitarios()."))
  NULL
}

#' Contagens de um sorteio de n com reposição, com probabilidade prob
#' @noRd
sorteio_contagens <- function(n, prob = NULL) tabulate(sample.int(n, n, replace = TRUE, prob = prob), n)

#' Multiplicadores Rao-Wu(-Yue) para uma árvore
#' @noRd
mult_raowu <- function(d, fpc = FALSE) {
  m <- numeric(d$n); oob <- logical(d$n)
  for (h in unique(d$estrato)) {
    lin <- which(d$estrato == h)
    upas <- unique(d$upa[lin]); nh <- length(upas)
    c_k <- tabulate(sample.int(nh, nh - 1L, replace = TRUE), nh)
    pos <- match(d$upa[lin], upas)
    if (fpc && !is.null(d$fpc)) {
      lam <- sqrt(1 - d$fpc[lin][1])
      m[lin] <- 1 - lam + lam * c_k[pos] * nh / (nh - 1)
    } else {
      m[lin] <- c_k[pos] * nh / (nh - 1)
      oob[lin] <- c_k[pos] == 0L
    }
  }
  list(m = m, oob = oob)
}

#' Bootstrap iid por unidade (Breiman)
#'
#' @param peso_na_arvore `FALSE` (padrão): floresta de Breiman, ignora o desenho.
#'   `TRUE`: árvores Hájek ("HJ" de Nalenz et al.), exige motor com corte ponderado.
#' @return Objeto `estrategia_rf`.
#' @export
reamostragem_iid <- function(peso_na_arvore = FALSE) {
  estrategia_reamostragem(
    nome = if (peso_na_arvore) "iid_hajek" else "iid",
    saida = if (peso_na_arvore) "pesos" else "contagens",
    usa_peso = peso_na_arvore, oob_por = "unidade",
    gerar = function(d, B) {
      m <- vapply(seq_len(B), function(b) sorteio_contagens(d$n), integer(d$n))
      m <- matrix(m, d$n, B)
      list(m = m, oob = m == 0)
    })
}

#' Bootstrap ponderado (WB): sorteio com probabilidade proporcional ao peso
#'
#' Árvores sem peso: o peso entra só no sorteio. Nos experimentos sintéticos
#' foi o esquema com menor viés residual (ver README). O OOB é por unidade e
#' tende a ser otimista com conglomerados.
#' @return Objeto `estrategia_rf`.
#' @export
reamostragem_wb <- function() {
  estrategia_reamostragem(
    nome = "wb", saida = "contagens", usa_peso = TRUE, oob_por = "unidade",
    gerar = function(d, B) {
      m <- matrix(vapply(seq_len(B), function(b) sorteio_contagens(d$n, d$peso), integer(d$n)), d$n, B)
      list(m = m, oob = m == 0)
    })
}

#' Rao-Wu(-Yue) por árvore com árvores Hájek
#'
#' Em cada estrato, n_h - 1 UPAs com reposição; multiplicador
#' c_k n_h / (n_h - 1). OOB = unidades de UPAs não sorteadas.
#'
#' @param fpc Usa a correção de população finita de Rao-Wu-Yue (exige `fpc` no
#'   desenho). Com fpc nenhum multiplicador é zero, então não há OOB.
#' @return Objeto `estrategia_rf` (saída `"pesos"`: exige corte ponderado).
#' @export
reamostragem_raowu <- function(fpc = FALSE) {
  estrategia_reamostragem(
    nome = if (fpc) "raowu_fpc" else "raowu", saida = "pesos", usa_peso = TRUE,
    oob_por = if (fpc) "nenhum" else "upa", requer = requer_duas_upas,
    gerar = function(d, B) {
      m <- matrix(0, d$n, B); oob <- matrix(FALSE, d$n, B)
      for (b in seq_len(B)) { r <- mult_raowu(d, fpc); m[, b] <- r$m; oob[, b] <- r$oob }
      list(m = m, oob = oob)
    })
}

#' Rao-Wu nas UPAs + bootstrap ponderado nas pessoas (padrão recomendado)
#'
#' Primeiro sorteia as UPAs por Rao-Wu; depois n sorteios de unidades com
#' probabilidade proporcional a peso x multiplicador. Árvores sem peso, então
#' roda no [motor_ranger()]. OOB por UPA. Nos experimentos: OOB honesto em
#' 18/18 cenários; viés residual entre o do Rao-Wu Hájek e o do WB.
#' @return Objeto `estrategia_rf`.
#' @export
reamostragem_raowu_wb <- function() {
  estrategia_reamostragem(
    nome = "raowu_wb", saida = "contagens", usa_peso = TRUE, oob_por = "upa",
    requer = requer_duas_upas,
    gerar = function(d, B) {
      m <- matrix(0L, d$n, B); oob <- matrix(FALSE, d$n, B)
      for (b in seq_len(B)) {
        r <- mult_raowu(d, FALSE)
        m[, b] <- sorteio_contagens(d$n, d$peso * r$m)
        oob[, b] <- r$oob
      }
      list(m = m, oob = oob)
    })
}

#' Réplicas do órgão produtor + bootstrap ponderado
#'
#' Para desenhos publicados só com pesos replicados (estrato e UPA mascarados).
#' A árvore b usa a réplica ((b - 1) mod R) + 1: n sorteios com probabilidade
#' proporcional ao peso replicado. OOB = unidades com peso replicado zero (em
#' BRR, jackknife e bootstrap isso coincide com UPAs inteiras). Réplicas do tipo
#' Fay não têm zeros e, portanto, não dão OOB. [Não testada no simulador.]
#' @return Objeto `estrategia_rf`.
#' @export
reamostragem_replicas <- function() {
  estrategia_reamostragem(
    nome = "replicas_orgao", saida = "contagens", usa_peso = TRUE, oob_por = "replica",
    requer = function(d) if (is.null(d$pesos_rep)) "O desenho n\u00E3o tem pesos replicados." else NULL,
    gerar = function(d, B) {
      R <- ncol(d$pesos_rep)
      m <- matrix(0L, d$n, B); oob <- matrix(FALSE, d$n, B)
      for (b in seq_len(B)) {
        pr <- d$pesos_rep[, ((b - 1L) %% R) + 1L]
        m[, b] <- sorteio_contagens(d$n, pr)
        oob[, b] <- pr == 0
      }
      list(m = m, oob = oob)
    })
}
