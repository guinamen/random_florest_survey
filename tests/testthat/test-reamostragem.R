estrategias_contagem <- list(reamostragem_iid(), reamostragem_wb(), reamostragem_raowu_wb())

test_that("todas as estratégias cumprem o contrato", {
  d <- desenho_brinquedo()
  for (e in c(estrategias_contagem, list(reamostragem_iid(TRUE), reamostragem_raowu()))) {
    set.seed(1)
    r <- rfdesenho:::gerar_reamostras(e, d, 7L)
    expect_equal(dim(r$m), c(d$n, 7L))
    expect_false(any(r$oob & r$m > 0))
  }
})

test_that("contagens somam n em cada árvore", {
  d <- desenho_brinquedo()
  for (e in estrategias_contagem) {
    set.seed(2)
    r <- rfdesenho:::gerar_reamostras(e, d, 5L)
    expect_true(all(colSums(r$m) == d$n))
  }
})

test_that("WB sorteia proporcionalmente ao peso", {
  d <- desenho_brinquedo()
  set.seed(3)
  r <- rfdesenho:::gerar_reamostras(reamostragem_wb(), d, 4000L)
  freq <- tapply(rowMeans(r$m), d$peso, mean)
  expect_equal(as.vector(freq / freq[1]), c(1, 2, 3, 4), tolerance = 0.03)
})

test_that("Rao-Wu: OOB é sempre a UPA inteira e o multiplicador tem média 1", {
  d <- desenho_brinquedo()
  set.seed(4)
  r <- rfdesenho:::gerar_reamostras(reamostragem_raowu(), d, 2000L)
  # dentro de cada UPA, todas as unidades têm o mesmo status OOB
  por_upa <- apply(r$oob, 2, function(o) tapply(o, d$upa, function(v) length(unique(v))))
  expect_true(all(por_upa == 1L))
  expect_equal(mean(r$m), 1, tolerance = 0.02)
  # exatamente uma UPA "a menos" sorteada por estrato: n_h - 1 sorteios
  expect_true(all(abs(colSums(r$m * d$peso) / sum(d$peso) - 1) < 1))
})

test_that("Rao-Wu + WB: nenhuma contagem em UPA fora da réplica", {
  d <- desenho_brinquedo()
  set.seed(5)
  r <- rfdesenho:::gerar_reamostras(reamostragem_raowu_wb(), d, 50L)
  expect_true(all(r$m[r$oob] == 0))
  expect_true(any(r$oob))
})

test_that("Rao-Wu recusa estrato com 1 UPA", {
  d <- desenho_rf(peso = rep(1, 4), estrato = c(1, 1, 2, 2), upa = c(1, 2, 1, 1))
  expect_error(rfdesenho:::gerar_reamostras(reamostragem_raowu(), d, 2L), "colapsar")
})

test_that("Rao-Wu com fpc não produz OOB", {
  d <- desenho_rf(peso = rep(1, 8), estrato = rep(1:2, each = 4), upa = rep(1:4, 2),
                  fpc = rep(0.5, 8))
  set.seed(6)
  r <- rfdesenho:::gerar_reamostras(reamostragem_raowu(fpc = TRUE), d, 10L)
  expect_false(any(r$oob))
  expect_true(all(r$m > 0))
})

test_that("réplicas do órgão: OOB onde o peso replicado é zero", {
  d <- desenho_brinquedo()
  pr <- matrix(d$peso, d$n, 3); pr[1:5, 1] <- 0; pr[6:10, 2] <- 0
  d2 <- desenho_rf(peso = d$peso, estrato = d$estrato, upa = d$upa_original, pesos_rep = pr)
  set.seed(7)
  r <- rfdesenho:::gerar_reamostras(reamostragem_replicas(), d2, 6L)
  expect_true(all(r$oob[1:5, c(1, 4)]))
  expect_true(all(r$m[1:5, c(1, 4)] == 0))
  expect_error(rfdesenho:::gerar_reamostras(reamostragem_replicas(), d, 2L), "replicados")
})

test_that("estratégia nova que viola o contrato é barrada", {
  ruim <- estrategia_reamostragem("ruim", saida = "contagens", usa_peso = FALSE,
                                  oob_por = "unidade",
                                  gerar = function(d, B) list(m = matrix(1.5, d$n, B),
                                                              oob = matrix(FALSE, d$n, B)))
  expect_error(rfdesenho:::gerar_reamostras(ruim, desenho_brinquedo(), 2L), "inteiras")
})
