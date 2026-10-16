# Árvore de regressão com corte e folhas ponderados (critério Hájek).
#
# Com pesos v_i, o ganho de um corte é a redução da soma de quadrados ponderada
#   SSE_w(pai) - SSE_w(esq) - SSE_w(dir),  SSE_w = sum v (y - ybar_v)^2,
# que é N_hat(pai) x Delta_HJ de Nalenz et al. (mesmo argmax). A folha recebe a
# média de Hájek sum(v y) / sum(v). Para y binário equivale ao Gini ponderado.
# Implementação em R puro: referência correta, adequada até alguns milhares de
# linhas. Para bases grandes, prefira motor_ranger() com estratégias de contagem.

#' Melhor corte de um nó
#' @noRd
melhor_corte_hajek <- function(Xn, yn, vn, feats, min_no) {
  W <- sum(vn); S1 <- sum(vn * yn); S2 <- sum(vn * yn^2)
  sse0 <- S2 - S1^2 / W
  melhor <- list(ganho = 1e-12 * max(sse0, 1e-12))
  nn <- length(yn)
  if (nn < 2L * min_no) return(melhor)
  k <- min_no:(nn - min_no)
  for (j in feats) {
    o <- order(Xn[, j])
    xs <- Xn[o, j]; ws <- vn[o]; ys <- yn[o]
    cw <- cumsum(ws); c1 <- cumsum(ws * ys); c2 <- cumsum(ws * ys^2)
    kk <- k[xs[k] < xs[k + 1L]]
    if (!length(kk)) next
    wl <- cw[kk]; wr <- W - wl
    ok <- wl > 0 & wr > 0
    if (!any(ok)) next
    kk <- kk[ok]; wl <- wl[ok]; wr <- wr[ok]
    g <- sse0 - (c2[kk] - c1[kk]^2 / wl) - ((S2 - c2[kk]) - (S1 - c1[kk])^2 / wr)
    i <- which.max(g)
    if (g[i] > melhor$ganho)
      melhor <- list(ganho = g[i], var = j, corte = (xs[kk[i]] + xs[kk[i] + 1L]) / 2)
  }
  melhor
}

#' Ajusta uma árvore Hájek
#'
#' @param X Matriz numérica.
#' @param y Resposta numérica.
#' @param v Pesos (zero = fora da árvore).
#' @param mtry Nº de variáveis sorteadas por nó.
#' @param min_no Mínimo de linhas distintas com peso > 0 em cada folha.
#' @param max_prof Profundidade máxima.
#' @return Lista com a árvore em vetores (var, corte, esq, dir, val).
#' @noRd
ajustar_arvore_hajek <- function(X, y, v, mtry, min_no, max_prof) {
  ativos <- which(v > 0)
  v <- v / mean(v[ativos])
  cap <- 2L * length(ativos) + 1L
  var <- integer(cap); corte <- numeric(cap); esq <- integer(cap); dir <- integer(cap)
  val <- numeric(cap)
  p <- ncol(X)
  pilha <- list(list(id = 1L, idx = ativos, prof = 0L))
  prox <- 2L
  while (length(pilha)) {
    nd <- pilha[[length(pilha)]]
    pilha[[length(pilha)]] <- NULL
    idx <- nd$idx
    val[nd$id] <- sum(v[idx] * y[idx]) / sum(v[idx])
    if (nd$prof >= max_prof || length(idx) < 2L * min_no) next
    b <- melhor_corte_hajek(X[idx, , drop = FALSE], y[idx], v[idx],
                            sample.int(p, mtry), min_no)
    if (is.null(b$var)) next
    vai_esq <- X[idx, b$var] <= b$corte
    var[nd$id] <- b$var; corte[nd$id] <- b$corte
    esq[nd$id] <- prox; dir[nd$id] <- prox + 1L
    pilha[[length(pilha) + 1L]] <- list(id = prox, idx = idx[vai_esq], prof = nd$prof + 1L)
    pilha[[length(pilha) + 1L]] <- list(id = prox + 1L, idx = idx[!vai_esq], prof = nd$prof + 1L)
    prox <- prox + 2L
  }
  k <- seq_len(prox - 1L)
  list(var = var[k], corte = corte[k], esq = esq[k], dir = dir[k], val = val[k])
}

#' Previsão de uma árvore Hájek
#' @noRd
prever_arvore_hajek <- function(arv, X) {
  no <- rep(1L, nrow(X))
  repeat {
    ii <- which(arv$var[no] > 0L)
    if (!length(ii)) break
    nj <- no[ii]
    vai_esq <- X[cbind(ii, arv$var[nj])] <= arv$corte[nj]
    no[ii] <- ifelse(vai_esq, arv$esq[nj], arv$dir[nj])
  }
  arv$val[no]
}
