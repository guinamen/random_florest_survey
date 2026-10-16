test_that("verdade_populacao bate com a média empírica na própria população", {
  pop <- simular_populacao(K_por_estrato = 150, semente = 1)
  p <- pop$pessoas[pop$pessoas$estrato %in% 1:2, ]
  v <- verdade_populacao(pop, p)
  # E[Y | X] tem a mesma média que Y e prevê Y com inclinação 1
  expect_lt(abs(mean(v) - mean(p$y)), 3 * sd(p$y - v) / sqrt(nrow(p)))
  expect_equal(unname(coef(lm(p$y ~ v))[2]), 1, tolerance = 0.05)
})

test_that("verdade binária fica em [0,1] e a prevalência é a pedida", {
  pop <- simular_populacao(K_por_estrato = 50, desfecho = "binario", prevalencia = 0.1, semente = 2)
  expect_equal(mean(pop$pessoas$y), 0.1, tolerance = 0.01)
  v <- verdade_populacao(pop, pop$pessoas[1:200, ])
  expect_true(all(v >= 0 & v <= 1))
})

test_that("sortear_amostra: pesos estimam N e o desenho é informativo", {
  pop <- simular_populacao(K_por_estrato = 100, semente = 3)
  tot <- replicate(30, sum(sortear_amostra(pop)$desenho$peso))
  expect_equal(mean(tot), nrow(pop$pessoas), tolerance = 0.05)
  am <- sortear_amostra(pop, semente = 4)
  expect_gt(mean(am$y), mean(pop$pessoas$y))   # kappa > 0 sobre-representa Z alto
})

test_that("inclusao_pps respeita o teto 1 e soma a", {
  pik <- rfdesenho:::inclusao_pps(c(100, 1, 1, 1, 1), 3)
  expect_equal(sum(pik), 3)
  expect_true(all(pik <= 1))
})
