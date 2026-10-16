test_that("suíte de conformidade roda e devolve os três critérios", {
  skip_if_not_installed("ranger")
  cf <- suite_conformidade(controle = controle_floresta(n_arvores = 15), R = 3,
                           n_teste = 200,
                           cenarios = list(controle = list(kappa = 0, K_por_estrato = 60),
                                           informativo = list(kappa = 0.7, K_por_estrato = 60)),
                           amostra = list(a_por_estrato = 10, m_por_upa = 10))
  expect_s3_class(cf, "conformidade_rf")
  expect_equal(nrow(cf$criterios), 3L)
  expect_equal(nrow(cf$replicas), 6L)
})

test_that("suíte completa: Rao-Wu + WB passa (lento; RFDESENHO_LENTO=1)", {
  skip_if_not_installed("ranger")
  skip_if(Sys.getenv("RFDESENHO_LENTO") != "1", "teste lento desativado")
  cf <- suite_conformidade(R = 30, nucleos = max(1L, parallel::detectCores() - 1L))
  print(cf)
  expect_true(cf$aprovado)
})
