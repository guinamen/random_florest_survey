# Bloco J (parte 2) — suíte de conformidade com critérios pré-registrados.
#
# Qualquer combinação nova de estratégia e motor deve passar aqui antes de ir
# para dados reais. Os critérios foram fixados ANTES de rodar e não devem ser
# ajustados depois de ver o resultado.
#
#   C1 (cenário controle, kappa = 0): o viés do método não difere do viés da
#      floresta de Breiman numa AAS de mesmo tamanho: |dif| <= 2 EP.
#   C2 (cenário informativo): a fração do viés ingênuo que NÃO é removida,
#      rho = (vies_m - vies_aas) / (vies_ing - vies_aas), tem rho + 2 EP < 0,5.
#   C3 (cenário informativo, se a estratégia tem qualquer OOB): o erro OOB não é
#      otimista frente ao erro em pessoas de UPAs fora da amostra:
#      z(erro_oob - erro_upas_novas) > -2.

#' Roda a suíte de conformidade
#'
#' @param estrategia,motor,controle Configuração a avaliar.
#' @param R Réplicas por cenário (use >= 50 para decisões; 20 para triagem).
#' @param n_teste Pessoas da população usadas para medir viés e erro.
#' @param cenarios Lista nomeada; cada elemento é uma lista de argumentos de
#'   [simular_populacao()]. Deve conter `controle` e `informativo`.
#' @param amostra Lista de argumentos de [sortear_amostra()].
#' @param nucleos Processos paralelos (fork).
#' @param semente Semente base.
#' @return Objeto `conformidade_rf` com `replicas`, `criterios` e `aprovado`.
#' @export
suite_conformidade <- function(estrategia = reamostragem_raowu_wb(), motor = motor_ranger(),
                               controle = controle_floresta(n_arvores = 100L),
                               R = 20L, n_teste = 1000L,
                               cenarios = list(controle = list(kappa = 0),
                                               informativo = list(kappa = 0.7)),
                               amostra = list(a_por_estrato = 20L, m_por_upa = 12L),
                               nucleos = 1L, semente = 1L) {
  checar(all(c("controle", "informativo") %in% names(cenarios)),
         "cenarios deve ter 'controle' e 'informativo'.")
  ingenua <- reamostragem_iid(FALSE)
  reps <- list()
  for (nm in names(cenarios)) {
    pop <- do.call(simular_populacao, c(cenarios[[nm]], list(semente = semente)))
    teste_idx <- com_semente(semente + 1L, sample.int(nrow(pop$pessoas), n_teste))
    teste <- pop$pessoas[teste_idx, ]
    Xt <- teste[, preditores_simulador(pop)]
    verdade <- verdade_populacao(pop, teste)
    uma <- function(r) {
      set.seed(semente * 1000L + r)
      am <- do.call(sortear_amostra, c(list(pop), amostra))
      aas <- amostra_aas(pop, length(am$y))
      fm <- floresta_desenho(am$X, am$y, am$desenho, estrategia, motor, controle)
      fi <- floresta_desenho(am$X, am$y, am$desenho, ingenua, motor, controle)
      fa <- floresta_desenho(aas$X, aas$y, aas$desenho, ingenua, motor, controle)
      pm <- predict(fm, Xt)
      nova <- !(teste$upa %in% am$desenho$upa_original)
      data.frame(cenario = nm, replica = r,
                 vies_m = mean(pm - verdade), vies_ing = mean(predict(fi, Xt) - verdade),
                 vies_aas = mean(predict(fa, Xt) - verdade),
                 rmse_m = sqrt(mean((pm - verdade)^2)),
                 gap_oob = erro_oob(fm) - mean((teste$y[nova] - pm[nova])^2))
    }
    res <- if (nucleos > 1L && .Platform$OS.type != "windows")
      parallel::mclapply(seq_len(R), uma, mc.cores = nucleos, mc.set.seed = FALSE)
    else lapply(seq_len(R), uma)
    reps[[nm]] <- do.call(rbind, res)
  }
  rep_df <- do.call(rbind, reps); rownames(rep_df) <- NULL
  crit <- avaliar_criterios(rep_df, tem_oob = estrategia$oob_por != "nenhum")
  structure(list(replicas = rep_df, criterios = crit,
                 aprovado = all(crit$veredito %in% c("passa", "n/a")),
                 configuracao = paste(estrategia$nome, "+", motor$nome)),
            class = "conformidade_rf")
}

#' Avalia C1-C3 a partir das réplicas
#' @noRd
avaliar_criterios <- function(d, tem_oob, n_boot = 500L) {
  ctl <- d[d$cenario == "controle", ]; inf <- d[d$cenario == "informativo", ]
  ep <- function(x) stats::sd(x) / sqrt(length(x))
  # C1
  dif <- ctl$vies_m - ctl$vies_aas
  c1_z <- mean(dif) / ep(dif)
  # C2: rho com EP por bootstrap das réplicas (pareadas)
  rho_de <- function(i) (mean(inf$vies_m[i]) - mean(inf$vies_aas[i])) /
                        (mean(inf$vies_ing[i]) - mean(inf$vies_aas[i]))
  rho <- rho_de(seq_len(nrow(inf)))
  rho_ep <- stats::sd(replicate(n_boot, rho_de(sample.int(nrow(inf), replace = TRUE))))
  # C3
  c3_z <- if (tem_oob && !anyNA(inf$gap_oob)) mean(inf$gap_oob) / ep(inf$gap_oob) else NA_real_
  data.frame(
    criterio = c("C1 controle: vi\u00E9s = vi\u00E9s AAS", "C2 informativo: rho + 2EP < 0,5",
                 "C3 informativo: OOB n\u00E3o otimista"),
    estatistica = c(sprintf("dif = %.4f, z = %.2f", mean(dif), c1_z),
                    sprintf("rho = %.3f, EP = %.3f", rho, rho_ep),
                    if (is.na(c3_z)) "sem OOB" else sprintf("gap = %.4f, z = %.2f", mean(inf$gap_oob), c3_z)),
    veredito = c(if (abs(c1_z) <= 2) "passa" else "falha",
                 if (rho + 2 * rho_ep < 0.5) "passa" else "falha",
                 if (is.na(c3_z)) "n/a" else if (c3_z > -2) "passa" else "falha"))
}

#' @export
print.conformidade_rf <- function(x, ...) {
  cat("<conformidade_rf>", x$configuracao, "| r\u00E9plicas por cen\u00E1rio:",
      max(x$replicas$replica), "\n")
  print(x$criterios, row.names = FALSE)
  cat(if (x$aprovado) "  => APROVADA\n" else "  => REPROVADA\n")
  invisible(x)
}
