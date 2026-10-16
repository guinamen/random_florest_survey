test_that("floresta_desenho de ponta a ponta com os dois motores", {
  skip_if_not_installed("ranger")
  s <- amostra_simulada()$am
  fl <- floresta_desenho(s$X, s$y, s$desenho, controle = controle_floresta(n_arvores = 30),
                         semente = 1)
  expect_s3_class(fl, "floresta_desenho")
  expect_length(predict(fl), length(s$y))
  expect_equal(dim(predict(fl, s$X[1:3, ], tipo = "arvores")), c(3L, 30L))
  fh <- floresta_desenho(s$X, s$y, s$desenho, reamostragem_raowu(), motor_hajek(),
                         controle_floresta(n_arvores = 5), semente = 1)
  expect_true(is.finite(erro_oob(fh)))
})

test_that("semente torna o ajuste reprodutível e não mexe no gerador global", {
  skip_if_not_installed("ranger")
  s <- amostra_simulada()$am
  set.seed(99); antes <- runif(1); set.seed(99)
  a <- floresta_desenho(s$X, s$y, s$desenho, controle = controle_floresta(n_arvores = 10), semente = 5)
  depois <- runif(1)
  b <- floresta_desenho(s$X, s$y, s$desenho, controle = controle_floresta(n_arvores = 10), semente = 5)
  expect_equal(predict(a), predict(b))
  expect_equal(antes, depois)
})

test_that("previsão OOB usa só árvores em que a unidade ficou fora", {
  skip_if_not_installed("ranger")
  s <- amostra_simulada()$am
  fl <- floresta_desenho(s$X, s$y, s$desenho, controle = controle_floresta(n_arvores = 20), semente = 2)
  P <- predict(fl, tipo = "arvores")
  i <- which(rowSums(fl$oob) > 0)[1]
  expect_equal(fl$previsao_oob[i], mean(P[i, fl$oob[i, ]]))
})

test_that("estimando recusa configurações incompatíveis", {
  skip_if_not_installed("ranger")
  s <- amostra_simulada()$am
  expect_error(floresta_desenho(s$X, s$y, s$desenho, reamostragem_iid(),
                                estimando = estimando("media_condicional")), "exige pesos")
  expect_error(floresta_desenho(s$X, s$y, s$desenho, estimando = estimando("dominio")),
               "media_hajek")
  expect_silent(floresta_desenho(s$X, s$y, s$desenho, reamostragem_iid(),
                                 controle = controle_floresta(n_arvores = 5),
                                 estimando = estimando("ranking")))
})

test_that("entradas inconsistentes são recusadas", {
  skip_if_not_installed("ranger")
  s <- amostra_simulada()$am
  expect_error(floresta_desenho(s$X[-1, ], s$y, s$desenho), "mesmo nº")
  X <- s$X; X$x1[1] <- NA
  expect_error(floresta_desenho(X, s$y, s$desenho), "ausentes")
})
