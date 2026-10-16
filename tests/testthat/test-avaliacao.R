test_that("auc_ponderada coincide com a AUC usual e com linhas replicadas", {
  set.seed(1)
  y <- rbinom(200, 1, .3); s <- y + rnorm(200)
  r <- rank(s); n1 <- sum(y); n0 <- sum(1 - y)
  auc_usual <- (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  expect_equal(auc_ponderada(y, s), auc_usual)
  w <- sample(1:3, 200, replace = TRUE); idx <- rep(seq_along(y), w)
  expect_equal(auc_ponderada(y, s, w), auc_ponderada(y[idx], s[idx]))
  expect_equal(auc_ponderada(c(0, 1), c(1, 1)), 0.5)   # empate
})

test_that("erro_oob é a média ponderada da perda", {
  skip_if_not_installed("ranger")
  s <- amostra_simulada()$am
  fl <- floresta_desenho(s$X, s$y, s$desenho, controle = controle_floresta(n_arvores = 20), semente = 1)
  ok <- !is.na(fl$previsao_oob)
  expect_equal(erro_oob(fl), weighted.mean((s$y[ok] - fl$previsao_oob[ok])^2, s$desenho$peso[ok]))
})

test_that("ajustar_grade devolve um erro por linha e avisa sobre OOB por unidade", {
  skip_if_not_installed("ranger")
  s <- amostra_simulada()$am
  g <- ajustar_grade(s$X, s$y, s$desenho, data.frame(min_no = c(5, 30)), n_arvores = 15)
  expect_equal(nrow(g), 2L)
  expect_true(all(is.finite(g$erro_oob)))
  expect_equal(nrow(attr(g, "melhor")), 1L)
  expect_warning(ajustar_grade(s$X, s$y, s$desenho, data.frame(min_no = 10),
                               estrategia = reamostragem_wb(), n_arvores = 5), "unidade")
})
