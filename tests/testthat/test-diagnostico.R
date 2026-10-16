test_that("diagnosticar mede o desenho e detecta informatividade", {
  s <- amostra_simulada(kappa = 1.2)$am
  dg <- diagnosticar(s$desenho, s$y, s$X)
  expect_s3_class(dg, "diagnostico_rf")
  expect_equal(dg$n_upas, 80L)
  expect_lt(dg$informatividade$p, 0.01)
  expect_gt(dg$icc, 0.05)
  r <- recomendar(dg, estimando("ranking"))
  expect_true(all(c("sinal", "valor", "acao", "base") %in% names(r)))
  expect_true("informatividade" %in% r$sinal)
})

test_that("ICC é ~0 sem conglomerado e informatividade não acusa sem relação", {
  set.seed(1)
  n <- 2000
  d <- desenho_rf(peso = runif(n, 1, 5), upa = rep(1:200, each = 10))
  y <- rnorm(n); X <- data.frame(x = rnorm(n))
  dg <- diagnosticar(d, y, X)
  expect_lt(dg$icc, 0.03)
  expect_gt(dg$informatividade$p, 0.001)
})

test_that("recomendar aponta estratos unitários e células com poucas UPAs", {
  d <- desenho_rf(peso = rep(1, 6), estrato = c(1, 1, 2, 3, 3, 3), upa = c(1, 2, 1, 1, 2, 3),
                  celula = c("a", "a", "a", "b", "b", "b"))
  r <- recomendar(diagnosticar(d))
  expect_true(any(grepl("1 UPA", r$sinal)))
  expect_true(any(grepl("calibração", r$sinal)))
})
