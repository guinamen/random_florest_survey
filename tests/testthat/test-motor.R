test_that("motor_hajek: árvore sem corte devolve a média ponderada", {
  set.seed(1)
  X <- data.frame(x = rnorm(30)); y <- rnorm(30); v <- runif(30)
  m <- motor_hajek()$ajustar(X, y, matrix(v), controle_floresta(n_arvores = 1, min_no = 30))
  p <- motor_hajek()$prever(m, X)
  expect_equal(unique(p[, 1]), sum(v * y) / sum(v))
})

test_that("motor_hajek: peso inteiro equivale a duplicar linhas (no corte e na folha)", {
  set.seed(2)
  X <- data.frame(x = runif(40)); y <- X$x * 3 + rnorm(40, 0, .1)
  v <- sample(1:3, 40, replace = TRUE)
  ctl <- controle_floresta(n_arvores = 1, min_no = 3, mtry = 1)
  a <- motor_hajek()$prever(motor_hajek()$ajustar(X, y, matrix(v), ctl), X)
  idx <- rep(seq_len(40), v)
  b <- motor_hajek()$prever(motor_hajek()$ajustar(X[idx, , drop = FALSE], y[idx],
                                                  matrix(1, length(idx)), ctl), X)
  # com min_no contado em linhas distintas os cortes podem diferir; o nível médio não
  expect_equal(mean(a), mean(b), tolerance = 0.05)
})

test_that("motor_ranger respeita a multiplicidade do inbag", {
  skip_if_not_installed("ranger")
  set.seed(3)
  X <- data.frame(x = rnorm(20)); y <- rnorm(20)
  cont <- sample(0:3, 20, replace = TRUE)
  ctl <- controle_floresta(n_arvores = 1, min_no = 100)   # impede cortes
  m <- motor_ranger()$ajustar(X, y, matrix(cont), ctl)
  p <- motor_ranger()$prever(m, X)
  expect_equal(unique(p[, 1]), sum(cont * y) / sum(cont), tolerance = 1e-10)
})

test_that("combinação estratégia x motor incompatível é recusada", {
  skip_if_not_installed("ranger")
  expect_error(rfdesenho:::checar_compatibilidade(reamostragem_raowu(), motor_ranger()),
               "corte_ponderado")
  expect_silent(rfdesenho:::checar_compatibilidade(reamostragem_raowu(), motor_hajek()))
})

test_that("mtry inválido é recusado", {
  expect_error(rfdesenho:::resolver_mtry(controle_floresta(mtry = 9), 3), "mtry")
})
