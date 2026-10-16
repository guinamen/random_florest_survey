# Dados pequenos e determinísticos compartilhados pelos testes.
desenho_brinquedo <- function() {
  # 3 estratos x 4 UPAs x 5 pessoas
  estrato <- rep(1:3, each = 20)
  upa <- rep(rep(1:4, each = 5), 3)          # rótulos repetidos entre estratos
  peso <- rep(c(10, 20, 30, 40), each = 5, times = 3)
  desenho_rf(peso = peso, estrato = estrato, upa = upa)
}

amostra_simulada <- function(K = 60, semente = 1, ...) {
  pop <- simular_populacao(K_por_estrato = K, semente = semente, ...)
  list(pop = pop, am = sortear_amostra(pop, a_por_estrato = 10, m_por_upa = 10,
                                       semente = semente + 1))
}
