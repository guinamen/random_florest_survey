test_that("media_hajek coincide com survey::svymean e svyby", {
  skip_if_not_installed("survey")
  s <- amostra_simulada()$am
  d <- s$desenho
  df <- data.frame(h = d$estrato, k = d$upa_original, w = d$peso, y = s$y,
                   g = s$X$x3)
  sd <- survey::svydesign(ids = ~k, strata = ~h, weights = ~w, data = df, nest = TRUE)
  m <- survey::svymean(~y, sd)
  mh <- media_hajek(s$y, rfdesenho::desenho_rf(d$peso, d$estrato, d$upa_original))
  expect_equal(mh$estimativa, unname(coef(m)), tolerance = 1e-10)
  expect_equal(mh$ep, as.vector(survey::SE(m)), tolerance = 1e-8)
  b <- survey::svyby(~y, ~g, sd, survey::svymean)
  md <- media_hajek(s$y, rfdesenho::desenho_rf(d$peso, d$estrato, d$upa_original), s$X$x3)
  expect_equal(md$estimativa, unname(b$y), tolerance = 1e-10)
  expect_equal(md$ep, unname(b$se), tolerance = 1e-8)
})

test_that("media_hajek com fpc coincide com survey", {
  skip_if_not_installed("survey")
  s <- amostra_simulada()$am
  d <- s$desenho
  df <- data.frame(h = d$estrato, k = d$upa_original, w = d$peso, y = s$y, f = d$fpc)
  sd <- survey::svydesign(ids = ~k, strata = ~h, weights = ~w, fpc = ~f, data = df, nest = TRUE)
  m <- survey::svymean(~y, sd)
  expect_equal(media_hajek(s$y, d)$ep, as.vector(survey::SE(m)), tolerance = 1e-8)
})

test_that("checar_consistencia, dependência parcial e importância funcionam", {
  skip_if_not_installed("ranger")
  s <- amostra_simulada()$am
  fl <- floresta_desenho(s$X, s$y, s$desenho, controle = controle_floresta(n_arvores = 30), semente = 1)
  cc <- checar_consistencia(fl, s$X$x3)
  expect_equal(nrow(cc), 2L)
  expect_true(all(is.finite(cc$z)))
  pd <- dependencia_parcial(fl, "x2", valores = c(-2, 0, 2))
  expect_true(pd$previsao_media[1] > pd$previsao_media[2])   # g contém x2^2
  imp <- importancia_oob(fl, semente = 1)
  expect_equal(imp$variavel[1], "x1")
})
