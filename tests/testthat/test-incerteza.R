test_that("incerteza_replicas devolve variâncias não negativas e intervalos coerentes", {
  skip_if_not_installed("ranger")
  s <- amostra_simulada()$am
  fl <- floresta_desenho(s$X, s$y, s$desenho, controle = controle_floresta(n_arvores = 30), semente = 1)
  inc <- incerteza_replicas(fl, s$X[1:4, ], n_replicas = 5, arvores_por_replica = 10, semente = 2)
  expect_equal(nrow(inc), 4L)
  expect_true(all(inc$var >= 0))
  expect_true(all(inc$inf <= inc$previsao & inc$previsao <= inc$sup))
  expect_match(attr(inc, "nota"), "E\\[f\\(x\\)\\]")
})

test_that("incerteza com pesos replicados do órgão", {
  skip_if_not_installed("ranger")
  s <- amostra_simulada()$am
  d <- s$desenho
  set.seed(3)
  pr <- vapply(1:6, function(r) { rp <- replica_raowu(d); tabulate(rp$idx, d$n) * d$peso * 20 / 19 },
               numeric(d$n))
  d2 <- desenho_rf(d$peso, d$estrato, d$upa_original, pesos_rep = pr, rep_tipo = "bootstrap",
                   rep_escala = 1 / 5)
  fl <- floresta_desenho(s$X, s$y, d2, reamostragem_replicas(),
                         controle = controle_floresta(n_arvores = 12), semente = 1)
  inc <- incerteza_replicas(fl, s$X[1:3, ], arvores_por_replica = 8, fonte = "orgao", semente = 1)
  expect_equal(nrow(inc), 3L)
  expect_true(all(is.finite(inc$var)))
})
