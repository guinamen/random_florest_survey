# Bloco B — objeto canônico do desenho amostral.
#
# O desenho é guardado SEPARADO dos dados de análise, mas alinhado linha a
# linha com eles: a linha i de X/y corresponde ao elemento i de cada vetor do
# desenho. Nenhum outro módulo lê nomes de variáveis da pesquisa; todos falam
# com este objeto.

#' Cria o objeto canônico do desenho amostral
#'
#' @param peso Vetor numérico de pesos amostrais (positivos).
#' @param estrato Vetor de estratos. `NULL` = desenho sem estratificação.
#' @param upa Vetor de unidades primárias de amostragem. `NULL` = cada linha é
#'   a sua própria UPA (amostragem de elementos). Os rótulos podem se repetir
#'   entre estratos: a identificação é sempre feita pelo par estrato:UPA.
#' @param fpc Fração amostral de UPAs no estrato (`f_h` em \[0, 1)), um valor
#'   por linha, constante dentro do estrato. `NULL` = 1º estágio tratado como
#'   com reposição (adequado quando `f_h` é pequeno).
#' @param celula Fator com a célula de calibração de cada linha (opcional).
#' @param totais Vetor nomeado com os totais populacionais por célula (opcional;
#'   exigido se `celula` for dado e se quiser recalibrar).
#' @param pesos_rep Matriz n x R de pesos replicados fornecidos pelo órgão
#'   produtor (opcional).
#' @param rep_tipo Tipo das réplicas (`"bootstrap"`, `"JK1"`, `"JKn"`, `"BRR"`,
#'   `"Fay"`, `"outro"`).
#' @param rep_escala,rep_rscales Constantes da variância replicada:
#'   `V = rep_escala * sum_r rep_rscales[r] * (t_r - t)^2`.
#' @return Objeto de classe `desenho_rf`.
#' @examples
#' d <- desenho_rf(peso = rep(10, 6), estrato = c(1, 1, 1, 2, 2, 2),
#'                 upa = c(1, 1, 2, 1, 2, 2))
#' d
#' @export
desenho_rf <- function(peso, estrato = NULL, upa = NULL, fpc = NULL,
                       celula = NULL, totais = NULL, pesos_rep = NULL,
                       rep_tipo = NULL, rep_escala = 1, rep_rscales = NULL) {
  n <- length(peso)
  checar(n >= 2L, "O desenho precisa de pelo menos 2 linhas.")
  checar(is.numeric(peso) && !anyNA(peso) && all(peso > 0),
         "peso deve ser num\u00E9rico, sem NA e estritamente positivo.")
  estrato <- if (is.null(estrato)) rep(1L, n) else estrato
  checar(length(estrato) == n && !anyNA(estrato), "estrato com tamanho errado ou NA.")
  estrato <- as.integer(factor(estrato))
  upa_orig <- if (is.null(upa)) seq_len(n) else upa
  checar(length(upa_orig) == n && !anyNA(upa_orig), "upa com tamanho errado ou NA.")
  # identificador único estrato:UPA (aceita rótulos de UPA aninhados no estrato)
  upa <- as.integer(factor(paste(estrato, upa_orig, sep = ":")))
  if (!is.null(fpc)) {
    checar(length(fpc) == n && all(fpc >= 0 & fpc < 1), "fpc deve estar em [0, 1).")
    checar(all(tapply(fpc, estrato, function(v) length(unique(v))) == 1L),
           "fpc deve ser constante dentro de cada estrato.")
  }
  if (!is.null(celula)) {
    celula <- as.factor(celula)
    checar(length(celula) == n && !anyNA(celula), "celula com tamanho errado ou NA.")
    if (!is.null(totais))
      checar(all(levels(celula) %in% names(totais)),
             "totais deve ter um valor para cada n\u00EDvel de celula.")
  }
  if (!is.null(pesos_rep)) {
    pesos_rep <- as.matrix(pesos_rep)
    checar(nrow(pesos_rep) == n, "pesos_rep deve ter uma linha por unidade.")
    checar(all(pesos_rep >= 0), "pesos_rep n\u00E3o pode ter valores negativos.")
    rep_tipo <- rep_tipo %||% "outro"
    rep_rscales <- rep_rscales %||% rep(1, ncol(pesos_rep))
  }
  structure(list(n = n, peso = as.numeric(peso), estrato = estrato, upa = upa,
                 upa_original = upa_orig, fpc = fpc, celula = celula, totais = totais,
                 pesos_rep = pesos_rep, rep_tipo = rep_tipo, rep_escala = rep_escala,
                 rep_rscales = rep_rscales),
            class = "desenho_rf")
}

#' Cria o desenho a partir de colunas de um data.frame
#'
#' Atalho para quando os nomes das variáveis de desenho são conhecidos (por
#' exemplo, no dicionário da pesquisa).
#'
#' @param dados data.frame com as variáveis de desenho.
#' @param peso,estrato,upa,fpc,celula Nomes (strings) das colunas. Apenas
#'   `peso` é obrigatório.
#' @param pesos_rep Vetor de nomes das colunas de pesos replicados (opcional).
#' @param ... Repassado a [desenho_rf()] (`totais`, `rep_tipo`, ...).
#' @return Objeto `desenho_rf`.
#' @examples
#' dd <- data.frame(w = c(2, 2, 3, 3), h = c(1, 1, 2, 2), k = c(1, 2, 1, 2))
#' desenho_de_colunas(dd, peso = "w", estrato = "h", upa = "k")
#' @export
desenho_de_colunas <- function(dados, peso, estrato = NULL, upa = NULL, fpc = NULL,
                               celula = NULL, pesos_rep = NULL, ...) {
  col <- function(nm) {
    if (is.null(nm)) return(NULL)
    checar(nm %in% names(dados), "Coluna '", nm, "' n\u00E3o existe nos dados.")
    dados[[nm]]
  }
  pr <- if (is.null(pesos_rep)) NULL else {
    falt <- setdiff(pesos_rep, names(dados))
    checar(!length(falt), "Colunas de pesos replicados ausentes: ", paste(falt, collapse = ", "))
    as.matrix(dados[, pesos_rep, drop = FALSE])
  }
  desenho_rf(peso = col(peso), estrato = col(estrato), upa = col(upa), fpc = col(fpc),
             celula = col(celula), pesos_rep = pr, ...)
}

#' Converte um objeto do pacote survey no desenho canônico
#'
#' Aceita `survey.design2` (criado por `survey::svydesign`) e `svyrep.design`
#' (criado por `survey::svrepdesign` ou `as.svrepdesign`). Usa só o primeiro
#' estágio (conglomerado último), como a variância linearizada usual.
#'
#' @param d Objeto de desenho do pacote survey.
#' @return Objeto `desenho_rf`.
#' @export
desenho_de_survey <- function(d) {
  if (inherits(d, "svyrep.design")) {
    w <- as.numeric(stats::weights(d, type = "sampling"))
    rw <- as.matrix(stats::weights(d, type = "analysis"))
    return(desenho_rf(peso = w, pesos_rep = rw, rep_tipo = d$type,
                      rep_escala = d$scale, rep_rscales = d$rscales))
  }
  checar(inherits(d, "survey.design2"),
         "Esperado survey.design2 ou svyrep.design; recebido ", class(d)[1], ".")
  w <- 1 / d$prob
  checar(all(is.finite(w)), "H\u00E1 probabilidades de inclus\u00E3o nulas no desenho.")
  estr <- d$strata[[1]]
  upa <- d$cluster[[1]]
  fpc <- NULL
  if (!is.null(d$fpc$popsize)) {
    N1 <- d$fpc$popsize[, 1]
    n1 <- d$fpc$sampsize[, 1]
    fpc <- n1 / N1
  }
  desenho_rf(peso = w, estrato = estr, upa = upa, fpc = fpc)
}

#' @export
print.desenho_rf <- function(x, ...) {
  nh <- tapply(x$upa, x$estrato, function(v) length(unique(v)))
  cat("<desenho_rf>", x$n, "linhas,", length(nh), "estratos,",
      length(unique(x$upa)), "UPAs\n")
  cat("  UPAs por estrato: m\u00EDn", min(nh), "| mediana", stats::median(nh),
      "| estratos com 1 UPA:", sum(nh < 2L), "\n")
  cat("  peso: soma", format(round(sum(x$peso)), big.mark = ".", decimal.mark = ","),
      "| CV", round(stats::sd(x$peso) / mean(x$peso), 3), "\n")
  if (!is.null(x$fpc)) cat("  fpc: f_h m\u00E1x", round(max(x$fpc), 4), "\n")
  if (!is.null(x$celula)) cat("  c\u00E9lulas de calibra\u00E7\u00E3o:", nlevels(x$celula), "\n")
  if (!is.null(x$pesos_rep))
    cat("  pesos replicados:", ncol(x$pesos_rep), "(", x$rep_tipo, ")\n")
  invisible(x)
}

#' Número de UPAs distintas por estrato
#' @param d Objeto `desenho_rf`.
#' @return Vetor nomeado (estrato -> nº de UPAs).
#' @export
upas_por_estrato <- function(d) vapply(split(d$upa, d$estrato), function(v) length(unique(v)), integer(1))

#' Colapsa estratos com uma única UPA
#'
#' Rao-Wu e a variância linearizada exigem pelo menos 2 UPAs por estrato. Cada
#' estrato unitário é fundido ao estrato seguinte na ordem dos códigos (ou ao
#' anterior, se for o último) — a mesma ideia de `options(survey.lonely.psu)`
#' com colapso. A fusão é registrada no atributo `"colapsos"`.
#'
#' @param d Objeto `desenho_rf`.
#' @return Objeto `desenho_rf` sem estratos unitários.
#' @export
colapsar_estratos_unitarios <- function(d) {
  log <- character(0)
  repeat {
    nh <- upas_por_estrato(d)
    uni <- as.integer(names(nh)[nh < 2L])
    if (!length(uni) || length(nh) < 2L) break
    h <- uni[1]
    cods <- sort(as.integer(names(nh)))
    pos <- match(h, cods)
    alvo <- if (pos < length(cods)) cods[pos + 1L] else cods[pos - 1L]
    d$estrato[d$estrato == h] <- alvo
    log <- c(log, sprintf("%d -> %d", h, alvo))
  }
  d$estrato <- as.integer(factor(d$estrato))
  if (!is.null(d$fpc) && length(log))
    warning("Estratos colapsados com fpc: o fpc do estrato fundido deixa de ser exato.",
            call. = FALSE)
  attr(d, "colapsos") <- log
  d
}

#' Subconjunto de linhas do desenho (com possíveis repetições)
#'
#' @param d Objeto `desenho_rf`.
#' @param idx Índices das linhas (podem repetir).
#' @param peso Pesos novos (opcional; padrão = pesos originais das linhas).
#' @param upa Identificadores de UPA novos (opcional).
#' @return Objeto `desenho_rf`.
#' @export
subconjunto_desenho <- function(d, idx, peso = NULL, upa = NULL) {
  out <- d
  out$n <- length(idx)
  out$peso <- peso %||% d$peso[idx]
  out$estrato <- d$estrato[idx]
  out$upa <- upa %||% d$upa[idx]
  out$upa_original <- d$upa_original[idx]
  if (!is.null(d$fpc)) out$fpc <- d$fpc[idx]
  if (!is.null(d$celula)) out$celula <- droplevels(d$celula[idx])
  if (!is.null(d$pesos_rep)) out$pesos_rep <- d$pesos_rep[idx, , drop = FALSE]
  out
}

#' Réplica física de Rao-Wu
#'
#' Em cada estrato sorteia n_h - 1 UPAs com reposição. Uma UPA sorteada c vezes
#' entra c vezes, com identificadores distintos (é tratada como c UPAs). Os
#' pesos são multiplicados por n_h/(n_h - 1). Serve para refazer o procedimento
#' inteiro numa réplica (incerteza por florestas-réplica).
#'
#' @param d Objeto `desenho_rf`.
#' @return `list(idx, peso, upa)`: linhas, pesos e UPAs da réplica.
#' @export
replica_raowu <- function(d) {
  idx <- integer(0); fat <- numeric(0); upa_nova <- integer(0)
  prox <- 0L
  for (h in unique(d$estrato)) {
    lin <- which(d$estrato == h)
    upas <- unique(d$upa[lin]); nh <- length(upas)
    if (nh < 2L) {
      idx <- c(idx, lin); fat <- c(fat, rep(1, length(lin)))
      upa_nova <- c(upa_nova, prox + match(d$upa[lin], upas)); prox <- prox + nh
      next
    }
    sorteio <- upas[sample.int(nh, nh - 1L, replace = TRUE)]
    for (c in seq_along(sorteio)) {
      s <- lin[d$upa[lin] == sorteio[c]]
      prox <- prox + 1L
      idx <- c(idx, s); fat <- c(fat, rep(nh / (nh - 1), length(s)))
      upa_nova <- c(upa_nova, rep(prox, length(s)))
    }
  }
  list(idx = idx, peso = d$peso[idx] * fat, upa = upa_nova)
}

#' Réplica a partir dos pesos replicados do órgão
#' @param d Objeto `desenho_rf` com `pesos_rep`.
#' @param r Índice da réplica.
#' @return `list(idx, peso, upa)` com as linhas de peso replicado positivo.
#' @export
replica_orgao <- function(d, r) {
  checar(!is.null(d$pesos_rep), "O desenho n\u00E3o tem pesos replicados.")
  pr <- d$pesos_rep[, r]
  idx <- which(pr > 0)
  list(idx = idx, peso = pr[idx], upa = d$upa[idx])
}
