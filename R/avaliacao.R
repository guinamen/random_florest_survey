# Bloco G — avaliação e ajuste de hiperparâmetros com OOB.
#
# Evidência (simulador): o OOB por UPA acompanhou o erro em UPAs novas; o OOB
# por unidade foi otimista mesmo sem desenho informativo. Por isso o ajuste usa
# o OOB da estratégia escolhida e avisa quando ele é por unidade.

#' Erro fora-da-bolsa ponderado pelo peso amostral
#'
#' @param fit Objeto `floresta_desenho`.
#' @param perda `"quadratica"`, `"absoluta"` ou `"logloss"` (y em 0/1).
#' @return Escalar: média ponderada da perda nas unidades com previsão OOB
#'   (`NA` se a estratégia não produz OOB).
#' @export
erro_oob <- function(fit, perda = c("quadratica", "absoluta", "logloss")) {
  perda <- match.arg(perda)
  ok <- !is.na(fit$previsao_oob)
  if (!any(ok)) return(NA_real_)
  y <- fit$y[ok]; p <- fit$previsao_oob[ok]; w <- fit$desenho$peso[ok]
  l <- switch(perda,
    quadratica = (y - p)^2,
    absoluta = abs(y - p),
    logloss = { p <- pmin(pmax(p, 1e-6), 1 - 1e-6); -(y * log(p) + (1 - y) * log(1 - p)) })
  media_pond(l, w)
}

#' AUC ponderada
#'
#' Probabilidade, ponderada pelos pesos, de um caso positivo ter escore maior
#' que um negativo (empates contam 1/2). Com pesos unitários é a AUC usual;
#' com pesos inteiros equivale a replicar as linhas.
#'
#' @param y Vetor 0/1.
#' @param escore Escore (maior = mais risco).
#' @param peso Pesos (padrão: 1).
#' @return Escalar em \[0, 1\].
#' @examples
#' auc_ponderada(c(0, 0, 1, 1), c(0.1, 0.4, 0.35, 0.8))
#' @export
auc_ponderada <- function(y, escore, peso = rep(1, length(y))) {
  checar(all(y %in% c(0, 1)), "y deve ser 0/1.")
  ok <- !is.na(escore)
  y <- y[ok]; escore <- escore[ok]; peso <- peso[ok]
  W1 <- sum(peso[y == 1]); W0 <- sum(peso[y == 0])
  if (W1 == 0 || W0 == 0) return(NA_real_)
  # por grupos de escore empatado: soma de pesos negativos abaixo + metade dos empatados
  g <- factor(escore)
  w0g <- tapply(peso * (y == 0), g, sum); w1g <- tapply(peso * (y == 1), g, sum)
  w0g[is.na(w0g)] <- 0; w1g[is.na(w1g)] <- 0
  abaixo <- cumsum(w0g) - w0g
  sum(w1g * (abaixo + w0g / 2)) / (W1 * W0)
}

#' AUC ponderada fora-da-bolsa
#' @param fit Objeto `floresta_desenho` com y binário.
#' @return Escalar.
#' @export
auc_oob <- function(fit) {
  ok <- !is.na(fit$previsao_oob)
  auc_ponderada(fit$y[ok], fit$previsao_oob[ok], fit$desenho$peso[ok])
}

#' Ajuste de hiperparâmetros por erro OOB
#'
#' Ajusta uma floresta para cada linha de `grade` e devolve o erro OOB.
#' Emite aviso se a estratégia tiver OOB por unidade (tende a ser otimista com
#' conglomerados) e erro se não houver OOB.
#'
#' @param X,y,desenho,estrategia,motor Como em [floresta_desenho()].
#' @param grade data.frame cujas colunas são argumentos de [controle_floresta()]
#'   (por exemplo `min_no`, `mtry`).
#' @param n_arvores Árvores por ajuste.
#' @param perda Perda passada a [erro_oob()]; ou `"auc"` (maior é melhor).
#' @param semente Semente base.
#' @param ... Outros argumentos fixos de [controle_floresta()].
#' @return data.frame: a grade com a coluna `erro_oob` (ou `auc_oob`) e o
#'   atributo `"melhor"` (linha escolhida).
#' @export
ajustar_grade <- function(X, y, desenho, grade, estrategia = reamostragem_raowu_wb(),
                          motor = motor_ranger(), n_arvores = 200L, perda = "quadratica",
                          semente = 1L, ...) {
  checar(is.data.frame(grade) && nrow(grade) >= 1, "grade deve ser um data.frame n\u00E3o vazio.")
  checar(estrategia$oob_por != "nenhum", "A estrat\u00E9gia '", estrategia$nome, "' n\u00E3o produz OOB.")
  if (estrategia$oob_por == "unidade")
    warning("OOB por unidade tende a subestimar o erro com conglomerados; prefira ",
            "reamostragem_raowu_wb() para ajustar.", call. = FALSE)
  res <- vapply(seq_len(nrow(grade)), function(i) {
    ctl <- do.call(controle_floresta, c(as.list(grade[i, , drop = FALSE]),
                                        list(n_arvores = n_arvores), list(...)))
    fl <- floresta_desenho(X, y, desenho, estrategia, motor, ctl, semente = semente + i)
    if (perda == "auc") auc_oob(fl) else erro_oob(fl, perda)
  }, numeric(1))
  out <- grade
  if (perda == "auc") out$auc_oob <- res else out$erro_oob <- res
  attr(out, "melhor") <- out[if (perda == "auc") which.max(res) else which.min(res), , drop = FALSE]
  out
}
