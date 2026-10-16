# Bloco I — saídas tipadas: estimativas por domínio (fora da floresta),
# checagem de consistência, dependência parcial e importância.

#' Média de Hájek com erro-padrão linearizado, por domínio
#'
#' Estimador clássico de pesquisas (não usa a floresta): é o que se deve usar
#' para prevalências e médias por domínio. Variância por linearização no
#' conglomerado último, 1º estágio com reposição (ou com fpc, se o desenho tiver).
#' Coincide com `survey::svymean`/`svyby` no mesmo desenho.
#'
#' @param y Variável numérica (0/1 para prevalência).
#' @param desenho Objeto `desenho_rf`.
#' @param dominio Vetor de domínios (opcional).
#' @return data.frame com `dominio`, `estimativa`, `ep`, `n`.
#' @examples
#' d <- desenho_rf(peso = c(1, 2, 1, 2), estrato = c(1, 1, 2, 2), upa = 1:4)
#' media_hajek(c(1, 0, 1, 1), d)
#' @export
media_hajek <- function(y, desenho, dominio = NULL) {
  checar(length(y) == desenho$n && !anyNA(y), "y com tamanho errado ou NA.")
  dominio <- if (is.null(dominio)) rep("total", desenho$n) else as.character(dominio)
  checar(length(dominio) == desenho$n, "dominio com tamanho errado.")
  w <- desenho$peso
  doms <- sort(unique(dominio))
  out <- lapply(doms, function(g) {
    I <- dominio == g
    Nd <- sum(w[I]); est <- sum(w[I] * y[I]) / Nd
    z <- I * w * (y - est) / Nd
    v <- var_linearizada(z, desenho$estrato, desenho$upa, desenho$fpc)
    data.frame(dominio = g, estimativa = est, ep = sqrt(v), n = sum(I))
  })
  do.call(rbind, out)
}

#' Checagem de consistência: previsões OOB contra a estimativa direta
#'
#' Compara a média ponderada de y com a média ponderada da previsão OOB, no
#' total e por domínio. Diferença grande (|z| > 2) indica que a floresta está
#' deslocando o nível do desfecho naquele domínio. O EP trata a previsão como
#' fixa (aproximação; ignora a variabilidade da própria floresta).
#'
#' @param fit Objeto `floresta_desenho`.
#' @param dominio Vetor de domínios alinhado às linhas de treino (opcional).
#' @return data.frame com `dominio`, `media_y`, `media_prev`, `diferenca`, `ep`, `z`.
#' @export
checar_consistencia <- function(fit, dominio = NULL) {
  d <- fit$desenho
  ok <- !is.na(fit$previsao_oob)
  checar(any(ok), "A floresta n\u00E3o tem previs\u00F5es OOB.")
  dominio <- if (is.null(dominio)) rep("total", d$n) else as.character(dominio)
  w <- d$peso
  doms <- sort(unique(dominio[ok]))
  out <- lapply(doms, function(g) {
    I <- ok & dominio == g
    e <- ifelse(I, fit$y - fit$previsao_oob, 0)
    Nd <- sum(w[I]); D <- sum(w * e) / Nd
    z <- I * w * (e - D) / Nd
    ep <- sqrt(var_linearizada(z, d$estrato, d$upa, d$fpc))
    data.frame(dominio = g, media_y = media_pond(fit$y[I], w[I]),
               media_prev = media_pond(fit$previsao_oob[I], w[I]),
               diferenca = D, ep = ep, z = D / ep)
  })
  do.call(rbind, out)
}

#' Dependência parcial ponderada pela população
#'
#' Para cada valor v da variável, fixa `X[, variavel] = v` em todas as linhas e
#' tira a média das previsões ponderada pelos pesos amostrais: a curva estima a
#' média POPULACIONAL de f, não a da amostra. \[Construção não validada no
#' simulador.\]
#'
#' @param fit Objeto `floresta_desenho`.
#' @param variavel Nome da coluna.
#' @param valores Valores a avaliar (padrão: quantis 5\% a 95\% ou níveis do fator).
#' @param n_valores Nº de valores quando `valores` é NULL.
#' @param ponderada Usar os pesos amostrais (padrão `TRUE`).
#' @return data.frame com `valor` e `previsao_media`.
#' @export
dependencia_parcial <- function(fit, variavel, valores = NULL, n_valores = 20L,
                                ponderada = TRUE) {
  X <- fit$X
  checar(variavel %in% names(X), "Vari\u00E1vel '", variavel, "' n\u00E3o est\u00E1 em X.")
  x <- X[[variavel]]
  if (is.null(valores))
    valores <- if (is.factor(x)) levels(x) else
      unique(stats::quantile(x, seq(0.05, 0.95, length.out = n_valores), names = FALSE))
  w <- if (ponderada) fit$desenho$peso else rep(1, nrow(X))
  pm <- vapply(seq_along(valores), function(i) {
    Xc <- X
    Xc[[variavel]] <- if (is.factor(x)) factor(rep(valores[i], nrow(X)), levels = levels(x))
                      else rep(valores[i], nrow(X))
    media_pond(predict(fit, Xc), w)
  }, numeric(1))
  data.frame(valor = valores, previsao_media = pm)
}

#' Importância por permutação fora-da-bolsa, ponderada
#'
#' Aumento da perda OOB ponderada quando a variável é permutada. Usa o OOB da
#' estratégia (por UPA em [reamostragem_raowu_wb()]).
#'
#' @param fit Objeto `floresta_desenho`.
#' @param variaveis Colunas a avaliar (padrão: todas).
#' @param n_perm Permutações por variável (média).
#' @param semente Semente.
#' @return data.frame ordenado com `variavel` e `importancia`.
#' @export
importancia_oob <- function(fit, variaveis = names(fit$X), n_perm = 1L, semente = NULL) {
  ok <- rowSums(fit$oob) > 0
  checar(any(ok), "A floresta n\u00E3o tem OOB.")
  base <- erro_oob(fit)
  w <- fit$desenho$peso
  perda_perm <- function(Xp) {
    P <- predict(fit, Xp, tipo = "arvores")
    pr <- rowSums(P * fit$oob)[ok] / rowSums(fit$oob)[ok]
    media_pond((fit$y[ok] - pr)^2, w[ok])
  }
  imp <- com_semente(semente, vapply(variaveis, function(v) {
    mean(vapply(seq_len(n_perm), function(k) {
      Xp <- fit$X; Xp[[v]] <- Xp[[v]][sample.int(nrow(Xp))]
      perda_perm(Xp)
    }, numeric(1))) - base
  }, numeric(1)))
  out <- data.frame(variavel = variaveis, importancia = unname(imp))
  out[order(-out$importancia), , drop = FALSE]
}
