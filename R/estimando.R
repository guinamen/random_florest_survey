# Bloco A — o estimando é declarado antes do ajuste e restringe a configuração.

#' Declara o estimando (a pergunta da análise)
#'
#' A configuração adequada depende da pergunta. Nos experimentos sintéticos:
#' para ranquear risco, os pesos quase não mudaram a AUC; para descrever a
#' população (média ou prevalência condicional a X), a floresta sem pesos teve
#' viés grande (+59\% a +88\% na prevalência) e os pesos removeram a maior parte.
#'
#' @param tipo Um de:
#'   \describe{
#'     \item{`"ranking"`}{ordenar unidades por risco (triagem, AUC).}
#'     \item{`"media_condicional"`}{estimar E\[Y | X\] na população finita.}
#'     \item{`"dominio"`}{média ou prevalência de um domínio. Não é tarefa da
#'       floresta: use [media_hajek()] ou o pacote survey.}
#'   }
#' @return Objeto `estimando_rf` com as exigências do tipo.
#' @examples
#' estimando("media_condicional")
#' @export
estimando <- function(tipo = c("media_condicional", "ranking", "dominio")) {
  tipo <- match.arg(tipo)
  exig <- switch(tipo,
    ranking = list(peso_obrigatorio = FALSE, metrica = "auc",
                   nota = "Pesos opcionais para ordenar; avalie com OOB por UPA."),
    media_condicional = list(peso_obrigatorio = TRUE, metrica = "quadratica",
                   nota = "Pesos obrigat\u00F3rios; o vi\u00E9s residual n\u00E3o some (rho ~ 0,2-0,3 no simulador)."),
    dominio = list(peso_obrigatorio = TRUE, metrica = NA_character_,
                   nota = "Use media_hajek() ou survey::svymean/svyby; a floresta n\u00E3o \u00E9 o estimador."))
  structure(c(list(tipo = tipo), exig), class = "estimando_rf")
}

#' @export
print.estimando_rf <- function(x, ...) {
  cat("<estimando_rf>", x$tipo, "| pesos obrigat\u00F3rios:", x$peso_obrigatorio, "\n ", x$nota, "\n")
  invisible(x)
}

#' Verifica se a estratégia é compatível com o estimando
#' @noRd
checar_estimando <- function(est, estrategia) {
  if (is.null(est)) return(invisible(TRUE))
  checar(inherits(est, "estimando_rf"), "estimando deve vir de estimando().")
  checar(est$tipo != "dominio",
         "Estimando 'dominio': use media_hajek() ou o pacote survey, n\u00E3o a floresta.")
  if (est$peso_obrigatorio)
    checar(isTRUE(estrategia$usa_peso),
           "O estimando '", est$tipo, "' exige pesos, mas a estrat\u00E9gia '",
           estrategia$nome, "' os ignora.")
  invisible(TRUE)
}
