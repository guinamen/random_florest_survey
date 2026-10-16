test_that("desenho_rf valida entradas", {
  expect_error(desenho_rf(peso = c(1, -1)), "positivo")
  expect_error(desenho_rf(peso = c(1, NA)), "positivo")
  expect_error(desenho_rf(peso = 1:4, estrato = 1:3), "estrato")
  expect_error(desenho_rf(peso = 1:4, fpc = c(.1, .2, .1, .2), estrato = c(1, 1, 2, 2)),
               "constante")
})

test_that("UPA é identificada pelo par estrato:UPA", {
  d <- desenho_brinquedo()
  expect_equal(length(unique(d$upa)), 12L)
  expect_equal(unname(upas_por_estrato(d)), rep(4L, 3))
})

test_that("desenho sem estrato e sem UPA trata cada linha como UPA", {
  d <- desenho_rf(peso = rep(2, 10))
  expect_equal(length(unique(d$upa)), 10L)
  expect_equal(length(unique(d$estrato)), 1L)
})

test_that("desenho_de_colunas reproduz desenho_rf", {
  dd <- data.frame(w = c(2, 2, 3, 3), h = c(1, 1, 2, 2), k = c(1, 2, 1, 2))
  a <- desenho_de_colunas(dd, peso = "w", estrato = "h", upa = "k")
  b <- desenho_rf(peso = dd$w, estrato = dd$h, upa = dd$k)
  expect_identical(a, b)
  expect_error(desenho_de_colunas(dd, peso = "nao_existe"), "não existe")
})

test_that("desenho_de_survey extrai estrato, UPA e pesos", {
  skip_if_not_installed("survey")
  s <- amostra_simulada()$am
  df <- data.frame(h = s$desenho$estrato, k = s$desenho$upa_original, w = s$desenho$peso)
  sd <- survey::svydesign(ids = ~k, strata = ~h, weights = ~w, data = df, nest = TRUE)
  d <- desenho_de_survey(sd)
  expect_equal(d$peso, df$w)
  expect_equal(length(unique(d$upa)), length(unique(paste(df$h, df$k))))
  rep <- survey::as.svrepdesign(sd, type = "bootstrap", replicates = 5)
  dr <- desenho_de_survey(rep)
  expect_equal(ncol(dr$pesos_rep), 5L)
})

test_that("colapsar_estratos_unitarios remove estratos com 1 UPA", {
  d <- desenho_rf(peso = rep(1, 6), estrato = c(1, 1, 2, 3, 3, 3), upa = c(1, 2, 1, 1, 2, 3))
  dc <- colapsar_estratos_unitarios(d)
  expect_true(all(upas_por_estrato(dc) >= 2L))
  expect_length(attr(dc, "colapsos"), 1L)
})

test_that("replica_raowu: n_h - 1 UPAs por estrato e peso total esperado preservado", {
  d <- desenho_brinquedo()
  set.seed(1)
  r <- replica_raowu(d)
  nupa <- tapply(r$upa, d$estrato[r$idx], function(v) length(unique(v)))
  expect_equal(as.vector(nupa), rep(3L, 3))           # cópias contam como UPAs distintas
  tot <- replicate(2000, sum(replica_raowu(d)$peso))
  expect_equal(mean(tot), sum(d$peso), tolerance = 0.02)
})

test_that("subconjunto_desenho mantém alinhamento", {
  d <- desenho_brinquedo()
  s <- subconjunto_desenho(d, c(1, 1, 30))
  expect_equal(s$n, 3L)
  expect_equal(s$peso, d$peso[c(1, 1, 30)])
})
