# rfdesenho

Florestas aleatórias para dados de pesquisas com desenho amostral complexo:
estratos, UPAs, pesos, calibração e pesos replicados.

O pacote separa três coisas que costumam vir misturadas:

- **quem são as unidades e como foram sorteadas** (o desenho);
- **quem entra em cada árvore** (a estratégia de reamostragem);
- **como a árvore é construída** (o motor).

Cada uma é um objeto com contrato explícito, e as combinações inválidas são
recusadas antes de rodar.

## Arquitetura

```
A. estimando()            -> a pergunta decide se o peso é obrigatório
B. desenho_rf()           -> objeto canônico (de colunas, de survey ou manual)
D. diagnosticar()         -> medidas do desenho
   recomendar()           -> regras de decisão auditáveis
E. reamostragem_*()       -> iid | wb | raowu | raowu_wb | replicas
F. motor_*()              -> ranger (contagens) | hajek (corte ponderado)
   floresta_desenho()     -> núcleo
G. erro_oob(), auc_oob(), ajustar_grade()
H. incerteza_replicas()
I. media_hajek(), checar_consistencia(), dependencia_parcial(), importancia_oob()
J. simular_populacao(), suite_conformidade()
```

| Arquivo | Bloco |
|---|---|
| `R/estimando.R` | A |
| `R/desenho.R` | B: construtores, colapso de estratos, réplicas físicas |
| `R/diagnostico.R` | D |
| `R/reamostragem.R` | E |
| `R/motor.R`, `R/arvore_hajek.R` | F |
| `R/floresta.R` | núcleo |
| `R/avaliacao.R` | G |
| `R/incerteza.R` | H |
| `R/saidas.R` | I |
| `R/simulador.R`, `R/conformidade.R` | J |

## Instalação e uso rápido

```r
# na pasta acima de rfdesenho/
install.packages(c("ranger", "survey", "testthat"))   # se faltar
devtools::install("rfdesenho")       # ou: R CMD INSTALL rfdesenho

library(rfdesenho)
pop <- simular_populacao(semente = 1)
am  <- sortear_amostra(pop, semente = 2)
fl  <- floresta_desenho(am$X, am$y, am$desenho,
                        estrategia = reamostragem_raowu_wb(), motor = motor_ranger(),
                        controle = controle_floresta(n_arvores = 300),
                        estimando = estimando("media_condicional"), semente = 3)
fl
checar_consistencia(fl)
```

Para dados reais, veja `inst/exemplos/pns2019.R`. Todos os nomes de variáveis
ficam num bloco `CONFIG` e **precisam ser conferidos no dicionário**.

## Qual configuração usar

| Estimando | Estratégia | Motor | Avaliação |
|---|---|---|---|
| ranking / triagem | `raowu_wb` (ou `iid`) | ranger | `auc_oob()` |
| média condicional na população | `raowu_wb` | ranger | `erro_oob()` |
| prevalência / média por domínio | não usar a floresta | — | `media_hajek()` ou survey |
| só pesos replicados (UPA mascarada) | `replicas` | ranger | `erro_oob()` [não testado] |

## Base empírica

Os números vêm de dados sintéticos, não de dados reais.

**Rodadas deste pacote:** 30 réplicas por cenário, população de 2 estágios com
`kappa = 0,7`. C1–C3 são os critérios da suíte de conformidade, descritos na
próxima seção.

| Configuração | ρ (viés não removido) | C3: OOB não otimista | Suíte |
|---|---|---|---|
| `iid` + ranger | 0,999 | passa (z = −1,25) | reprovada em C2 |
| `wb` + ranger | 0,223 | **falha** (z = −6,63) | reprovada em C3 |
| `raowu` + hajek | 0,319 | passa | aprovada |
| `raowu_wb` + ranger | 0,326–0,334 | passa | aprovada |
| `raowu_wb` + hajek | 0,240 | passa | aprovada |

A diferença de ρ entre os dois motores do `raowu_wb` vem da semântica do
`min_no`. Variando o `min_no`, os dois motores caem na mesma fronteira
viés × RMSE: o `min_no = 2k` no ranger equivale a `min_no = k` no Hájek. Por
isso, escolha o tamanho da folha pelo OOB, não por um número fixo.

**Experimentos anteriores** (simulador de 3 estágios, em `rf_desenho/`):

- os pesos são indispensáveis para descrever a população, mas quase não mudam a AUC;
- a fração amostral pequena não piora o viés relativo;
- a calibração aumentou o viés residual (ρ +0,07);
- a incerteza por florestas-réplica Rao-Wu subestima a variância em ~20%;
- intervalos para a função verdadeira não atingem a cobertura nominal (0,44–0,67).

## Suíte de conformidade

Os critérios foram fixados antes de rodar:

- **C1 (controle não informativo):** o viés do método não difere do da floresta
  de Breiman numa AAS de mesmo tamanho (|z| ≤ 2).
- **C2 (informativo):** ρ = (viés_m − viés_AAS)/(viés_ingênua − viés_AAS)
  satisfaz ρ + 2EP < 0,5.
- **C3 (informativo, se houver OOB):** o erro OOB não é otimista frente ao erro
  em pessoas de UPAs fora da amostra (z > −2).

Toda estratégia ou motor novo deve passar aqui antes de ir para dados reais:

```r
suite_conformidade(minha_estrategia, meu_motor, R = 50, nucleos = 8)
```

## Como estender

- **Estratégia nova** (por exemplo, o bootstrap de pseudo-população da SDRF):
  `estrategia_reamostragem(nome, gerar = function(d, B) list(m, oob), saida, usa_peso, oob_por)`.
  O contrato é validado a cada uso.
- **Motor novo:** `motor_arvores(nome, capacidades, ajustar, prever)`.

## Testes

```bash
Rscript -e 'devtools::test()'                       # ~1 min
RFDESENHO_LENTO=1 Rscript -e 'devtools::test()'     # inclui a suíte completa
```

Os testes cobrem:

- o contrato de cada estratégia: dimensões, OOB por UPA inteira, contagens
  inteiras, WB proporcional ao peso;
- a multiplicidade do `inbag` no ranger;
- a equivalência peso inteiro = linha duplicada;
- a igualdade de `media_hajek()` com `survey::svymean`/`svyby`, com e sem fpc;
- a verdade do simulador e a reprodutibilidade por semente.

## Limitações

- O ranger usa o peso só no sorteio, não no critério de corte. Árvores Hájek
  exigem `motor_hajek()`, em R puro, que é lento acima de alguns milhares de linhas.
- **Memória:** cerca de 3 × n × B × 8 bytes. Com n = 90 mil e B = 300, ~650 MB.
- O teste de informatividade (DuMouchel-Duncan) ignora o conglomerado e é
  liberal com ICC alto: serve como triagem.
- `dependencia_parcial()`, `importancia_oob()` e `reamostragem_replicas()` não
  foram validados no simulador.
- O simulador tem 2 estágios. Desenhos em duas fases, amostragem de Poisson e
  ajuste de não resposta pesado não foram testados.
