# Núcleo — liga desenho, estratégia de reamostragem e motor de árvores.

#' Ajusta uma floresta aleatória que respeita o desenho amostral
#'
#' Fluxo: (1) valida entradas, estimando e compatibilidade estratégia x motor;
#' (2) a estratégia gera, para cada árvore, as contagens ou os multiplicadores
#' e a máscara OOB; (3) o motor ajusta as árvores; (4) calcula a previsão OOB.
#'
#' @param X data.frame (ou matriz) de preditores, sem NA.
#' @param y Resposta numérica (binária como 0/1).
#' @param desenho Objeto [desenho_rf()] alinhado às linhas de `X`.
#' @param estrategia Objeto `estrategia_rf`. Padrão: [reamostragem_raowu_wb()].
#' @param motor Objeto `motor_rf`. Padrão: [motor_ranger()].
#' @param controle Lista de [controle_floresta()].
#' @param estimando Objeto [estimando()] (opcional, mas recomendado): recusa
#'   configurações incompatíveis com a pergunta.
#' @param semente Semente (opcional) para reprodutibilidade; não altera o estado
#'   global do gerador.
#' @return Objeto `floresta_desenho`.
#' @examples
#' if (requireNamespace("ranger", quietly = TRUE)) {
#'   pop <- simular_populacao(K_por_estrato = 60, semente = 1)
#'   am <- sortear_amostra(pop, semente = 2)
#'   fl <- floresta_desenho(am$X, am$y, am$desenho,
#'                          controle = controle_floresta(n_arvores = 50),
#'                          estimando = estimando("media_condicional"), semente = 3)
#'   erro_oob(fl)
#' }
#' @export
floresta_desenho <- function(X, y, desenho, estrategia = reamostragem_raowu_wb(),
                             motor = motor_ranger(), controle = controle_floresta(),
                             estimando = NULL, semente = NULL) {
  X <- como_preditores(X)
  checar(inherits(desenho, "desenho_rf"), "desenho deve vir de desenho_rf().")
  checar(inherits(estrategia, "estrategia_rf"), "estrategia inv\u00E1lida.")
  checar(inherits(motor, "motor_rf"), "motor inv\u00E1lido.")
  checar(is.numeric(y) && !anyNA(y), "y deve ser num\u00E9rico e sem NA.")
  checar(length(y) == nrow(X) && nrow(X) == desenho$n,
         "X, y e desenho precisam ter o mesmo n\u00BA de linhas.")
  checar_estimando(estimando, estrategia)
  checar_compatibilidade(estrategia, motor)

  com_semente(semente, {
    r <- gerar_reamostras(estrategia, desenho, controle$n_arvores)
    V <- if (estrategia$saida == "pesos") r$m * desenho$peso else r$m
    modelo <- motor$ajustar(X, y, V, controle)
    P <- motor$prever(modelo, X)
  })
  n_oob <- rowSums(r$oob)
  prev_oob <- ifelse(n_oob > 0, rowSums(P * r$oob) / pmax(n_oob, 1), NA_real_)

  structure(list(modelo = modelo, motor = motor, estrategia = estrategia,
                 controle = controle, estimando = estimando, desenho = desenho,
                 X = X, y = y, oob = r$oob, previsao_oob = prev_oob,
                 chamada = match.call()),
            class = "floresta_desenho")
}

#' Previsão de uma floresta_desenho
#'
#' @param object Objeto `floresta_desenho`.
#' @param novos data.frame com as mesmas colunas de `X` (padrão: dados de treino).
#' @param tipo `"media"` (média das árvores) ou `"arvores"` (matriz n x B).
#' @param ... Ignorado.
#' @return Vetor ou matriz de previsões.
#' @export
predict.floresta_desenho <- function(object, novos = NULL, tipo = c("media", "arvores"), ...) {
  tipo <- match.arg(tipo)
  novos <- if (is.null(novos)) object$X else como_preditores(novos)
  falt <- setdiff(names(object$X), names(novos))
  checar(!length(falt), "Faltam colunas em 'novos': ", paste(falt, collapse = ", "))
  P <- object$motor$prever(object$modelo, novos[, names(object$X), drop = FALSE])
  if (tipo == "arvores") P else rowMeans(P)
}

#' @export
print.floresta_desenho <- function(x, ...) {
  cat("<floresta_desenho>", x$controle$n_arvores, "\u00E1rvores |", x$estrategia$nome, "+",
      x$motor$nome, "| n =", x$desenho$n, "| p =", ncol(x$X), "\n")
  if (!is.null(x$estimando)) cat("  estimando:", x$estimando$tipo, "\n")
  cobertos <- mean(!is.na(x$previsao_oob))
  cat("  OOB por", x$estrategia$oob_por, "| unidades com previs\u00E3o OOB:",
      sprintf("%.1f%%", 100 * cobertos), "\n")
  e <- erro_oob(x)
  if (!is.na(e)) cat("  erro OOB (quadr\u00E1tico, ponderado):", signif(e, 4), "\n")
  invisible(x)
}
