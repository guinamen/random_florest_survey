# Bloco F — motores de árvores.
#
# Contrato: um `motor_rf` tem
#   capacidades : subconjunto de c("contagens", "corte_ponderado")
#   ajustar(X, y, V, controle) -> modelo   (V: matriz n x B, uma coluna por árvore)
#   prever(modelo, X)          -> matriz n x B (previsão de cada árvore)
# "contagens": aceita V inteira como multiplicidade das linhas na árvore.
# "corte_ponderado": usa V como peso no critério de corte e na folha.

#' Construtor genérico de motor de árvores
#'
#' @param nome Rótulo.
#' @param capacidades Vetor com `"contagens"` e/ou `"corte_ponderado"`.
#' @param ajustar Função `(X, y, V, controle) -> modelo`.
#' @param prever Função `(modelo, X) -> matriz n x B`.
#' @return Objeto `motor_rf`.
#' @export
motor_arvores <- function(nome, capacidades, ajustar, prever) {
  capacidades <- match.arg(capacidades, c("contagens", "corte_ponderado"), several.ok = TRUE)
  structure(list(nome = nome, capacidades = capacidades, ajustar = ajustar, prever = prever),
            class = "motor_rf")
}

#' @export
print.motor_rf <- function(x, ...) {
  cat("<motor_rf>", x$nome, "| capacidades:", paste(x$capacidades, collapse = ", "), "\n")
  invisible(x)
}

#' Parâmetros de controle da floresta
#'
#' @param n_arvores Nº de árvores.
#' @param mtry Variáveis por nó (`NULL` = max(1, floor(p/3))).
#' @param min_no Tamanho mínimo da folha. No `motor_ranger` vira `min.bucket`
#'   e conta repetições do sorteio; no `motor_hajek` conta linhas distintas com
#'   peso > 0. Com pesos variáveis as escalas diferem: no simulador, `min_no = 2k`
#'   no ranger caiu na mesma fronteira viés x RMSE que `min_no = k` no Hájek.
#'   Escolha pelo erro OOB ([ajustar_grade()]), não por valor fixo.
#' @param max_prof Profundidade máxima (`NULL` = sem limite no ranger; 30 no Hájek).
#' @param num_threads Threads do ranger.
#' @return Lista de controle.
#' @export
controle_floresta <- function(n_arvores = 500L, mtry = NULL, min_no = 10L,
                              max_prof = NULL, num_threads = 1L) {
  checar(n_arvores >= 1, "n_arvores deve ser >= 1.")
  checar(min_no >= 1, "min_no deve ser >= 1.")
  list(n_arvores = as.integer(n_arvores), mtry = mtry, min_no = as.integer(min_no),
       max_prof = max_prof, num_threads = as.integer(num_threads))
}

resolver_mtry <- function(controle, p) {
  m <- controle$mtry %||% max(1L, floor(p / 3))
  checar(m >= 1 && m <= p, "mtry deve estar entre 1 e o n\u00BA de preditores (", p, ").")
  as.integer(m)
}

#' Motor ranger (rápido; árvores sem peso)
#'
#' Passa as contagens de cada árvore pelo argumento `inbag` do ranger. O peso
#' amostral, quando houver, já entrou no sorteio da estratégia. Resposta binária
#' 0/1 é tratada como regressão (a previsão é uma probabilidade).
#' @return Objeto `motor_rf`.
#' @export
motor_ranger <- function() {
  checar(requireNamespace("ranger", quietly = TRUE), "Instale o pacote 'ranger'.")
  motor_arvores(
    nome = "ranger", capacidades = "contagens",
    ajustar = function(X, y, V, controle) {
      inbag <- lapply(seq_len(ncol(V)), function(b) as.integer(V[, b]))
      fit <- ranger::ranger(x = X, y = as.numeric(y), num.trees = ncol(V),
                     mtry = resolver_mtry(controle, ncol(X)),
                     min.node.size = 2L * controle$min_no, min.bucket = controle$min_no,
                     max.depth = controle$max_prof, replace = TRUE, inbag = inbag,
                     num.threads = controle$num_threads, oob.error = FALSE,
                     seed = sample.int(.Machine$integer.max, 1L))
      list(ranger = fit, num_threads = controle$num_threads)
    },
    prever = function(modelo, X) {
      p <- stats::predict(modelo$ranger, data = X, predict.all = TRUE,
                          num.threads = modelo$num_threads)$predictions
      matrix(p, nrow(X))
    })
}

#' Motor Hájek (referência em R; corte ponderado)
#'
#' Aceita pesos não inteiros (estratégias Rao-Wu e iid com peso na árvore).
#' Lento para bases grandes: use para validar e para comparar com o ranger.
#' Fatores são convertidos em códigos inteiros (tratados como ordenados).
#' @return Objeto `motor_rf`.
#' @export
motor_hajek <- function() {
  motor_arvores(
    nome = "hajek", capacidades = c("contagens", "corte_ponderado"),
    ajustar = function(X, y, V, controle) {
      Xm <- como_matriz_num(X)
      mtry <- resolver_mtry(controle, ncol(Xm))
      prof <- controle$max_prof %||% 30L
      arvs <- lapply(seq_len(ncol(V)), function(b)
        ajustar_arvore_hajek(Xm, as.numeric(y), V[, b], mtry, controle$min_no, prof))
      list(arvores = arvs, colunas = colnames(Xm))
    },
    prever = function(modelo, X) {
      Xm <- como_matriz_num(X)[, modelo$colunas, drop = FALSE]
      matrix(vapply(modelo$arvores, function(a) prever_arvore_hajek(a, Xm), numeric(nrow(Xm))),
             nrow(Xm))
    })
}

#' Verifica se motor e estratégia combinam
#' @noRd
checar_compatibilidade <- function(estrategia, motor) {
  precisa <- if (estrategia$saida == "pesos") "corte_ponderado" else "contagens"
  checar(precisa %in% motor$capacidades,
         "A estrat\u00E9gia '", estrategia$nome, "' precisa de um motor com '", precisa,
         "', mas o motor '", motor$nome, "' s\u00F3 oferece: ",
         paste(motor$capacidades, collapse = ", "), ".")
  invisible(TRUE)
}
