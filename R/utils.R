# Utilitários internos: checagens, sementes locais e estatísticas ponderadas.

#' Interrompe com mensagem se a condição falhar
#' @noRd
checar <- function(cond, ...) {
  if (!isTRUE(cond)) stop(paste0(...), call. = FALSE)
  invisible(TRUE)
}

#' Avalia `expr` com semente local, restaurando o estado do gerador depois
#'
#' Se `semente` for NULL, apenas avalia. Evita efeitos colaterais no
#' `.Random.seed` global do usuário.
#' @noRd
com_semente <- function(semente, expr) {
  if (is.null(semente)) return(expr)
  tem <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (tem) antigo <- get(".Random.seed", envir = globalenv())
  on.exit({
    if (tem) assign(".Random.seed", antigo, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv())
  })
  set.seed(semente)
  expr
}

#' Média de Hájek
#' @noRd
media_pond <- function(y, w) sum(w * y) / sum(w)

#' Tamanho efetivo de Kish: (sum w)^2 / sum w^2
#' @noRd
n_eff_kish <- function(w) sum(w)^2 / sum(w^2)

#' Tamanho efetivo em UPAs: (sum_k W_k)^2 / sum_k W_k^2
#' @noRd
n_eff_upa <- function(w, upa) {
  W <- tapply(w, upa, sum)
  sum(W)^2 / sum(W^2)
}

#' Variância linearizada (conglomerado último, 1º estágio com reposição)
#'
#' Para um total de escores z (já linearizados), soma por estrato
#' n_h/(n_h-1) * sum_k (z_hk - mean z_h)^2, com fator (1 - f_h) se houver fpc.
#' Estratos com uma UPA contribuem zero (o chamador deve colapsar antes).
#' @noRd
var_linearizada <- function(z, estrato, upa, fpc = NULL) {
  zk <- tapply(z, upa, sum)
  hk <- tapply(estrato, upa, function(v) v[1])
  fk <- if (is.null(fpc)) NULL else tapply(fpc, upa, function(v) v[1])
  tot <- 0
  for (h in unique(hk)) {
    sel <- which(hk == h)
    nh <- length(sel)
    if (nh < 2L) next
    v <- zk[sel]
    fat <- nh / (nh - 1) * (if (is.null(fk)) 1 else (1 - fk[sel][1]))
    tot <- tot + fat * sum((v - mean(v))^2)
  }
  unname(tot)
}

#' ICC por ANOVA de um fator (UPA)
#' @noRd
icc_anova <- function(y, g) {
  g <- as.factor(g)
  k <- nlevels(g); n <- length(y)
  if (k < 2L || n <= k) return(NA_real_)
  ni <- tabulate(g)
  mi <- tapply(y, g, mean)
  msb <- sum(ni * (mi - mean(y))^2) / (k - 1)
  msw <- sum((y - mi[g])^2) / (n - k)
  n0 <- (n - sum(ni^2) / n) / (k - 1)
  s2b <- max(0, (msb - msw) / n0)
  if (s2b + msw == 0) return(0)
  s2b / (s2b + msw)
}

#' Converte preditores em data.frame; checa NA
#' @noRd
como_preditores <- function(X) {
  if (is.matrix(X)) X <- as.data.frame(X)
  checar(is.data.frame(X), "X deve ser data.frame ou matriz.")
  checar(ncol(X) >= 1L, "X precisa de pelo menos uma coluna.")
  checar(!anyNA(X), "X tem valores ausentes; trate-os antes (imputa\u00E7\u00E3o ou exclus\u00E3o).")
  if (is.null(names(X)) || any(names(X) == "")) names(X) <- paste0("x", seq_len(ncol(X)))
  X
}

#' Matriz numérica a partir de preditores (fatores viram códigos inteiros)
#' @noRd
como_matriz_num <- function(X) {
  M <- data.matrix(X)
  storage.mode(M) <- "double"
  M
}

`%||%` <- function(a, b) if (is.null(a)) b else a
