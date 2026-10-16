# =============================================================================
# Exemplo: floresta com desenho na PNS 2019 (morador selecionado)
#
# ATENÇÃO: os nomes de variáveis abaixo são CANDIDATOS e precisam ser
# conferidos no dicionário da PNS 2019 antes de rodar. Nenhum foi verificado.
# Todo nome de variável da pesquisa está concentrado no bloco CONFIG: o resto
# do script fala só com o objeto de desenho e com X/y.
#
# Uso:  Rscript pns2019.R  (a partir de uma pasta com o pacote instalado)
# =============================================================================
library(rfdesenho)

CONFIG <- list(
  arquivo      = "~/dados/pns2019/pns2019.rds",
  estrato      = "V0024",     # estrato                        [CONFIRMAR]
  upa          = "UPA_PNS",   # UPA                             [CONFIRMAR]
  peso         = "V00291",    # peso do morador selecionado     [CONFIRMAR: com ou sem calibração?]
  filtro       = function(d) d[!is.na(d$V00291) & d$V00291 > 0, ],  # só moradores selecionados [CONFIRMAR]
  dominio      = "V0001",     # UF, para estimativas por domínio [CONFIRMAR]
  nucleos      = max(1L, parallel::detectCores() - 1L)
)

# ---- 0. Ingestão --------------------------------------------------------------
pns <- readRDS(path.expand(CONFIG$arquivo))
pns <- CONFIG$filtro(pns)

# ---- 1. Desenho ---------------------------------------------------------------
des <- desenho_de_colunas(pns, peso = CONFIG$peso, estrato = CONFIG$estrato, upa = CONFIG$upa)
des <- colapsar_estratos_unitarios(des)
print(des)

# ---- 2. Desfecho e preditores (DEFINIR com o projeto) --------------------------
# Exemplo de forma, não de conteúdo: y binário 0/1 e preditores sem NA.
# y <- as.numeric(pns$<variavel_DRC> == 1)
# X <- data.frame(idade = pns$<idade>, sexo = factor(pns$<sexo>), ...)
stop("Defina y e X no bloco 2 antes de continuar.")

# ---- 3. Diagnóstico -> recomendações ----------------------------------------------
est <- estimando("media_condicional")      # ou "ranking" para triagem
dg <- diagnosticar(des, y, X)
print(dg)
print(recomendar(dg, est))

# ---- 4. Estimativas por domínio: estimador direto, não a floresta ------------------
print(media_hajek(y, des, pns[[CONFIG$dominio]]))

# ---- 5. Ajuste de hiperparâmetros pelo OOB por UPA ----------------------------------
grade <- expand.grid(min_no = c(10, 20, 40, 80), mtry = c(2, 4))
g <- ajustar_grade(X, y, des, grade, estrategia = reamostragem_raowu_wb(),
                   motor = motor_ranger(), n_arvores = 200, num_threads = CONFIG$nucleos)
print(g); melhor <- attr(g, "melhor")

# ---- 6. Floresta final ---------------------------------------------------------------
fl <- floresta_desenho(X, y, des, reamostragem_raowu_wb(), motor_ranger(),
                       controle_floresta(n_arvores = 300, min_no = melhor$min_no,
                                         mtry = melhor$mtry, num_threads = CONFIG$nucleos),
                       estimando = est, semente = 2026)
print(fl)
if (all(y %in% 0:1)) print(auc_oob(fl))

# ---- 7. Saídas ---------------------------------------------------------------------
print(checar_consistencia(fl, pns[[CONFIG$dominio]]))   # |z| > 2: floresta desloca o nível
print(importancia_oob(fl, semente = 1))
# print(dependencia_parcial(fl, "idade"))

# ---- 8. Incerteza de desenho para perfis de interesse --------------------------------
# perfis <- X[c(...), ]
# print(incerteza_replicas(fl, perfis, n_replicas = 50, arvores_por_replica = 50,
#                          nucleos = CONFIG$nucleos, semente = 1))
