#' rfdesenho: florestas aleatórias para pesquisas com desenho complexo
#'
#' O pacote organiza o ajuste em blocos independentes, cada um com um contrato
#' explícito:
#'
#' \describe{
#'   \item{A. Estimando}{[estimando()] declara a pergunta (ranking, média
#'     condicional na população finita ou estimativa por domínio). A pergunta
#'     decide se os pesos são obrigatórios.}
#'   \item{B. Desenho}{[desenho_rf()], [desenho_de_colunas()] e
#'     [desenho_de_survey()] produzem um objeto canônico (estrato, UPA, peso,
#'     fpc, células de calibração, pesos replicados) paralelo às linhas dos dados.}
#'   \item{D. Diagnóstico}{[diagnosticar()] mede o desenho e [recomendar()]
#'     transforma as medidas em decisões auditáveis.}
#'   \item{E. Reamostragem}{estratégias trocáveis com a mesma interface:
#'     [reamostragem_iid()], [reamostragem_wb()], [reamostragem_raowu()],
#'     [reamostragem_raowu_wb()], [reamostragem_replicas()].}
#'   \item{F. Motor}{[motor_ranger()] (rápido, árvores sem peso) e
#'     [motor_hajek()] (referência em R, corte ponderado). Cada motor declara
#'     capacidades; combinações inválidas são recusadas.}
#'   \item{Núcleo}{[floresta_desenho()] liga estratégia e motor.}
#'   \item{G. Avaliação}{[erro_oob()], [auc_ponderada()], [ajustar_grade()].}
#'   \item{H. Incerteza}{[incerteza_replicas()]: florestas-réplica.}
#'   \item{I. Saídas}{[dependencia_parcial()], [importancia_oob()],
#'     [media_hajek()], [checar_consistencia()].}
#'   \item{J. Conformidade}{[simular_populacao()] e [suite_conformidade()]:
#'     critérios pré-registrados que qualquer combinação nova deve passar.}
#' }
#'
#' @section Base empírica:
#' As escolhas padrão vêm de experimentos com dados sintéticos (PPS em 3
#' estágios, 8 estratos). As direções dos efeitos são o que mais se transfere;
#' as magnitudes, o que menos. Ver o README.
#'
#' @keywords internal
#' @importFrom stats predict
"_PACKAGE"
