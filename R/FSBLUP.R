#' @title solve the mix model and calculate ebvs
#'
#' @param phe phenotype data.frame
#' @param trait_col trait col num in phe
#' @param k kinship matrix
#'
#' @importFrom rrBLUP mixed.solve
FSBLUP.Mix <- function(
    phe, trait_col, k
)
{
    require(rrBLUP)

    y <- colnames(phe)[trait_col]
    res_k <- mixed.solve(phe$y, K = k)

    return(res_k$u)
}

FSBLUP.version <- function() {
cat("#----------------------------------------#\n")
cat("#  _____ ____  ____  _    _   _ ____     #\n")
cat("#  |  ___/ ___|| __ )| |  | | | |  _ \\   #\n")
cat("#  | |_  \\___ \\|  _ \\| |  | | | | |_) |  #\n")
cat("#  |  _|  ___) | |_) | |__| |_| |  __/   #\n")
cat("#  |_|   |____/|____/|_____\\___/|_|      #\n")
cat("#                                        #\n")
cat("#----------------------------------------#\n")
}


#' @title grid search for a given range of points
#'
#' @param fn function to calculate indicator
#' @param x search ranges, e.g. (0,1)
#' @param y search ranges, e.g. (0,1)
#' @param point_num search density on two directions, steps is `sqrt(point_num)`
#' @param phe_ab phenotype data.frame
#' @param trait analysis trait
#' @param M1 pedigree kinship matrix
#' @param M2 genomic matrix
#' @param M3 omics matrix
#' @param a parameter `\alpha`
#' @param b parameter `\beta`
#' @param trait analysis trait
#' @param crv.num fold number of validation
#' @param crv.rep.num repeat times in each fold, also as random seed
#' @param stas.phe.col target trait when calculate the indicator, default as @inheritParams trait
#' @param stas.type indicator when optimize the fusion matrix
#' @param stas.fn custom calculated indicators
#' @param train.id reference population id
#' @param test.id validation population id
#' @param ... passed parameters for other function
#'
#' @return a data.frame induces search points and cor
#' @importFrom purrr map2_dfr
#'
FSBLUP.GridSearch <- function(
    fn, x = c(0, 1), y = c(0, 1), point_num = 25, phe_ab, trait, M1, M2, M3, cov, crv.num, crv.rep.num,
    stas.phe.col, stas.type, stas.fn, train.id = NULL, test.id = NULL, ...
)
{
  start_time1 <- Sys.time()
  require(purrr)
  step <- ceiling(sqrt(point_num))
  para_df <- expand.grid(
    a = seq(x[1], x[2], length.out = step),
    b = seq(y[1], y[2], length.out = step)
  )

  #acc_df <- map2_dfr(para_df$a, para_df$b, ~fn(.x, .y, ...))
  if(length(train.id) == 0 & length(test.id) == 0)
  {
    acc_df <- map2_dfr(para_df$a, para_df$b, ~fn(phe_ab = phe_ab, trait = trait, M1 = M1, M2 = M2, M3 = M3, cov = cov,
                                                 a = .x, b = .y, crv.num = crv.num, crv.rep.num = crv.rep.num,
                                                 stas.phe.col = stas.phe.col, stas.type = stas.type, stas.fn = stas.fn), .progress = T)
  } else if(!(length(train.id) == 0 & length(test.id) == 0))
  {
    acc_df <- map2_dfr(para_df$a, para_df$b, ~fn(phe_ab = phe_ab, trait = trait, M1 = M1, M2 = M2, M3 = M3, cov = cov,
                                                 a = .x, b = .y, RP.id = train.id, VP.id = test.id,
                                                 stas.phe.col = stas.phe.col, stas.type = stas.type, stas.fn = stas.fn), .progress = T)
  } else {
    stop("\n\n parameters error in grid search !!! \n\n")
  }

  cat("\nGrid search finish with ", difftime(Sys.time(), start_time1, units = "mins"), " minutes \n")
  return(acc_df)
}


#' @title bisection optimize points from grid_search returns
#'
#' @param acc_df grid_search returns
#' @param fn function to calculate cor
#' @param threshold stop if reached between two iterations
#' @param phe_ab phenotype data.frame
#' @param trait analysis trait
#' @param M1 pedigree kinship matrix
#' @param M2 genomic matrix
#' @param M3 omics matrix
#' @param crv.num fold number of validation
#' @param crv.rep.num repeat times in each fold, also as random seed
#' @param crv.num fold number of validation
#' @param crv.rep.num repeat times in each fold, also as random seed
#' @param stas.phe.col target trait when calculate the indicator, default as @inheritParams trait
#' @param stas.type indicator when optimize the fusion matrix
#' @param stas.fn custom calculated indicators
#' @param train.id reference population id
#' @param test.id validation population id
#' @param ... passed parameters for other function
#' @param max_iter max bisection iterations
#'
#' @return the best points
#' @importFrom purrr map_dfr
#'

FSBLUP.Bisection <- function(
  acc_df, fn, max_iter = 10, threshold = 1e-4, phe_ab, trait, M1, M2, M3, cov, crv.num, crv.rep.num,
  stas.phe.col, stas.type, stas.fn, train.id = NULL, test.id = NULL, ...
)
{
  start_time1 <- Sys.time()
  require(purrr)

  best_index <- which.max(acc_df$cor)
  best_a <- acc_df$a[best_index]
  best_b <- acc_df$b[best_index]
  best_cor <- acc_df$cor[best_index]

  step_a <- diff(range(acc_df$a)) / (nrow(acc_df) - 1)
  step_b <- diff(range(acc_df$b)) / (nrow(acc_df) - 1)

  cli::cli_progress_bar("Bisection", total = max_iter)
  for (iter_num_run in 1:max_iter) {

    new_points <- expand.grid(a = c(best_a - step_a, best_a + step_a),
                              b = c(best_b - step_b, best_b + step_b))

    if(length(train.id) == 0 & length(test.id) == 0)
    {
      new_cor <- map_dfr(1:nrow(new_points), ~ fn(a = new_points$a[.], b = new_points$b[.],
                                                  phe_ab = phe_ab, trait = trait, M1 = M1, M2 = M2, M3 = M3, cov = cov,
                                                  crv.num = crv.num, crv.rep.num = crv.rep.num,
                                                  stas.phe.col = stas.phe.col, stas.type = stas.type, stas.fn = stas.fn))
      cli::cli_progress_update()
    } else if(!(length(train.id) == 0 & length(test.id) == 0))
    {
      new_cor <- map_dfr(1:nrow(new_points), ~ fn(a = new_points$a[.], b = new_points$b[.],
                                                  phe_ab = phe_ab, trait = trait, M1 = M1, M2 = M2, M3 = M3, cov = cov,
                                                  RP.id = train.id, VP.id = test.id,
                                                  stas.phe.col = stas.phe.col, stas.type = stas.type, stas.fn = stas.fn))
      cli::cli_progress_update()
    } else {
      stop("\n\n parameters error in bisection !!! \n\n")
    }


    new_best_index <- which.max(new_cor$cor)
    new_best_a <- new_points$a[new_best_index]
    new_best_b <- new_points$b[new_best_index]
    new_best_cor <- new_cor$cor[new_best_index]


    if (abs(new_best_cor - best_cor) < threshold) {
      break
    }

    if ((new_best_cor - best_cor) > 0.05) {

      step_a <- step_a / 2
      step_b <- step_b / 2
    } else if ((new_best_cor - best_cor) < 0.05) {
      step_a <- step_a * 0.75
      step_b <- step_b * 0.75
    } else {
      step_a <- step_a / 2
      step_b <- step_b / 2
    }

    best_a <- new_best_a
    best_b <- new_best_b
    best_cor <- new_best_cor
  }

  cli::cli_progress_done()
  cat("\nBisection finish with ", difftime(Sys.time(), start_time1, units = "mins"), " minutes \n")

  return(list(optimal_a = best_a, optimal_b = best_b, optimal_cor = best_cor))
}


#' @title adjust matrix to positive definite
adj_pos <- function(x) {
  if (!is.null(x)) {
    diag(x) <- diag(x) + 1e-06
    B <- try(chol(x), silent = TRUE)
    if (inherits(B, what = "try-error")) {
      x <- Matrix::nearPD(x)$mat %>% as.matrix()
      rm(B)
    }
  }
  return(x)
}



#' Combine matrix
#'
#' @param M1 A pedigree matrix
#' @param M2 G genomic matrix
#' @param M3 O omics matrix
#' @param a weight for M3
#' @param b weight for M2
#' @param phe_c phenotype data.frame
#'
#' @return matrix
#'

FSBLUP.CombineMat <- function(
  phe_c, M1, M2, M3, a, b
)
{
  id_1 <- colnames(M1)
  id_2 <- colnames(M2)
  id_3 <- colnames(M3)

  C.mode <- NULL

  C.mode <- ifelse(length(id_1) == length(id_2) && length(id_2) == length(id_3), "A",
                   ifelse(any(length(id_1) == length(id_2), length(id_1) == length(id_3), length(id_2) == length(id_3)), "B", "C"))

  #cat(dim(M1),dim(M2),dim(M3))
  AA000 <- function(phe_c, M1, M2, M3, a, b){
    M2 <- M2[rownames(M1), colnames(M1)]
    M3 <- M3[rownames(M1), colnames(M1)]
    k <- a*M3+b*M2+(1-a-b)*M1
    return(k)
  }

  BB000 <- function(phe_c, M1, M2, M3, a, b){
    id1 <- id_2
    id2 <- id_3
    id1 <- id1[!id1 %in% id2]
    M2_11 <- M2[id1, id1]
    M2_12 <- M2[id1, id2]
    M2_21 <- M2[id2, id1]
    M2_22 <- M2[id2, id2]
    M2_22 <- adj_pos(M2_22)

    iM2_22 <- solve(M2_22)

    M3 <- M3[id2, id2]

    nind <- nrow(M3)
    avg_sum_1 <- sum(M3) / (nind * nind)
    avg_sum_2 <- sum(M2_22) / (nind * nind)
    avg_diag_1 <- sum(diag(M3)) / nind
    avg_diag_2 <- sum(diag(M2_22)) / nind
    sw <- (avg_sum_2 - avg_diag_2) / (avg_sum_1 - avg_diag_1)
    sm <- avg_diag_1 - sw * avg_diag_2
    M3 <- sw*M3+sm

    M3 <- a*M3+(1-a)*M2_22

    H11 <- M2_11 + M2_12 %*% iM2_22 %*% (M3 - M2_22) %*% iM2_22 %*% M2_21
    H12 <- M2_12 %*% iM2_22 %*% M3
    H21 <- M3 %*% iM2_22 %*% M2_21
    H <- cbind(rbind(H11, H21), rbind(H12, M3))
    H <- H / mean(diag(H))

    M1 <- M1[rownames(H), colnames(H)]
    H <- b*H+(1-b)*M1
    #H <- adj_pos(H)

    return(H)
  }

  CC000 <- function(phe_c, M1, M2, M3, a, b){
    id1 <- id_1
    id2 <- id_2
    id1 <- id1[!id_1 %in% id2]
    M1_11 <- M1[id1, id1]
    M1_12 <- M1[id1, id2]
    M1_21 <- M1[id2, id1]
    M1_22 <- M1[id2, id2]
    M1_22 <- adj_pos(M1_22)

    iM1_22 <- solve(M1_22)

    M2 <- M2[id2, id2]

    nind <- nrow(M2)
    avg_sum_1 <- sum(M2) / (nind * nind)
    avg_sum_2 <- sum(M1_22) / (nind * nind)
    avg_diag_1 <- sum(diag(M2)) / nind
    avg_diag_2 <- sum(diag(M1_22)) / nind
    sw <- (avg_sum_2 - avg_diag_2) / (avg_sum_1 - avg_diag_1)
    sm <- avg_diag_1 - sw * avg_diag_2
    M3 <- sw*M2+sm

    M3 <- b*M2+(1-b)*M1_22

    H11 <- M1_11 + M1_12 %*% iM1_22 %*% (M2 - M1_22) %*% iM1_22 %*% M1_21
    H12 <- M1_12 %*% iM1_22 %*% M2
    H21 <- M2 %*% iM1_22 %*% M1_21
    H <- cbind(rbind(H11, H21), rbind(H12, M2))
    H <- H / mean(diag(H))

    id1 <- id_1
    id2 <- id_3
    id1 <- id1[!id1 %in% id2]
    H_11 <- H[id1, id1]
    H_12 <- H[id1, id2]
    H_21 <- H[id2, id1]
    H_22 <- H[id2, id2]
    H_22 <- adj_pos(H_22)

    iH_22 <- solve(H_22)

    M3 <- M3[id2, id2]

    nind <- nrow(M3)
    avg_sum_1 <- sum(M3) / (nind * nind)
    avg_sum_2 <- sum(H_22) / (nind * nind)
    avg_diag_1 <- sum(diag(M3)) / nind
    avg_diag_2 <- sum(diag(H_22)) / nind
    sw <- (avg_sum_2 - avg_diag_2) / (avg_sum_1 - avg_diag_1)
    sm <- avg_diag_1 - sw * avg_diag_2
    M3 <- sw*M3+sm

    M3 <- a*M3+(1-a)*H_22

    H11 <- H_11 + H_12 %*% iH_22 %*% (M3 - H_22) %*% iH_22 %*% H_21
    H12 <- H_12 %*% iH_22 %*% M3
    H21 <- M3 %*% iH_22 %*% H_21
    H <- cbind(rbind(H11, H21), rbind(H12, M3))
    H <- H / mean(diag(H))
  }

  if (C.mode == "A") {
    kk <- AA000(phe_c = phe_c, M1 = M1, M2 = M2, M3 = M3, a = a, b = b)
    kk <- tryCatch({kk = adj_pos(kk)}, error = function(e) {diag(kk) = diag(M1); return(kk)})
  } else if (C.mode == "B") {
    kk <- BB000(phe_c = phe_c, M1 = M1, M2 = M2, M3 = M3, a = a, b = b)
    kk <- tryCatch({kk = adj_pos(kk)}, error = function(e) {diag(kk) = diag(M1); return(kk)})
  } else if (C.mode == "C") {
    kk <- CC000(phe_c = phe_c, M1 = M1, M2 = M2, M3 = M3, a = a, b = b)
    kk <- tryCatch({kk = adj_pos(kk)}, error = function(e) {diag(kk) = diag(M1); return(kk)})
  } else {
    stop("Invalid C.mode")
  }

  return(kk)
}

#' @title statistical indicator
#' @description
#' Note x1 and x2 must in the same order
#'
#' @param x1 phenotype
#' @param x2 predicted GEBV
#' @param type indicator type, e.g. "cor", "auc", "rmse", "other"
#' @param fn self custom function, valid only when "type = other"
#'
#' @return value
FSBLUP.stas.cal <- function(
    x1, x2, type=c("cor", "auc", "rmse", "other"), fn=NULL
) {
  # x1 phenotype
  # x2 predicted GEBV
  switch(
    match.arg(type),
    cor = {
      val <- cor(x1, x2, use="pairwise.complete.obs")
    },
    auc = {
      if(sum(c(0,1) %in% unique(x1)) != 2) stop("Only two levels(case/1,control/0) are allowed!")
      X <- cbind(x1,x2)
      X <- X[order(X[,2]),]
      N.case <- sum(as.numeric(X[,1]) == 1)
      N.cont <- sum(as.numeric(X[,1]) == 0)
      case.index <- which(as.numeric(X[,1]) == 1)
      case.index.mean <- mean(case.index)
      val <- (1/N.cont)*(case.index.mean-0.5*N.case-0.5)
    },
    rmse = {
      val <- sqrt(mean((x1 - x2)^2))
    },
    other = {
      if(type == "other" & is.null(fn)) stop("When stas type is other, should provide function to statistic")
      val <- fn(x1, x2)
    }
  )
  return(val)
}


#' Cross validation calculation for accuracy
#'
#' @param phe_ab phenotype data.frame
#' @param M1 pedigree kinship matrix
#' @param M2 genomic matrix
#' @param M3 omics matrix
#' @param a parameter `\alpha`
#' @param b parameter `\beta`
#' @param trait analysis trait
#' @param crv.num fold number of validation
#' @param crv.rep.num repeat times in each fold, also as random seed
#' @param stas.phe.col target trait when calculate the indicator, default as @inheritParams trait
#' @param stas.type indicator when optimize the fusion matrix
#' @param stas.fn custom calculated indicators
#' @param ... passed parameters for other function
#'
#' @return list
#'

cv.cal_ab <- function(
    phe_ab, trait, M1, M2, M3, cov, a, b, crv.num, crv.rep.num, stas.phe.col = NULL, stas.type, stas.fn, ...
)
{
  suppressMessages(require(rrBLUP))

  if (missing(phe_ab)) stop("Missing phe_ab")
  if (missing(M1)) stop("Missing M1")
  if (missing(M2)) stop("Missing M2")
  if (missing(M3)) stop("Missing M3")
  if (missing(a)) stop("Missing a")
  if (missing(b)) stop("Missing b")
  if (missing(crv.num)) stop("Missing crv.num")
  if (missing(crv.rep.num)) stop("Missing crv.rep.num")

  colnames(phe_ab)[1] <- "id"

  kk <- FSBLUP.CombineMat(phe_c = phe_ab, M1 = M1, M2 = M2, M3 = M3, a = a, b = b)

  kk <- kk[as.character(phe_ab$id), as.character(phe_ab$id)]

  results <- data.frame()
  for (xi in seq_len(crv.rep.num)) {
    for(xj in seq_len(crv.num)) {
      set.seed(xi)
      phe_ab$xpar <- sample(seq_len(crv.num), size = nrow(phe_ab), replace = TRUE,
                            prob = c(rep((1/crv.num), times = crv.num)))
      phe_ab$y0 <- pull(phe_ab, trait)
      phe_ab$yNA0 <- pull(phe_ab, trait)
      xnas <- phe_ab$xpar == xj
      phe_ab$yNA0[xnas] <- NA

      res_k <- tryCatch(mixed.solve(phe_ab$yNA0, X = cov, K = kk), error = {kk=adj_pos(kk);mixed.solve(phe_ab$yNA0, X = cov, K = kk)})
      res_k <- res_k$u[xnas]

      if(is.null(stas.phe.col)) x1oo00oo <- phe_ab$y0[xnas]
      if(!is.null(stas.phe.col)) x1oo00oo <- pull(phe_ab, stas.phe.col)[xnas]

      results <- rbind(
        results,
        data.frame(
          a = a,
          b = b,
          cv_seed = xi,
          cvnum = xj,
          pearson = FSBLUP.stas.cal(x1oo00oo, res_k, type = stas.type, fn = stas.fn)
          #pearson = cor(phe_ab$y0[xnas], res_k, use = "pairwise.complete.obs")
        )
      )
    }
  }

  acc <- results %>%
    group_by(a, b) %>%
    summarise(cor = mean(pearson, na.rm = T), .groups = "drop") %>%
    ungroup() %>%
    arrange(desc(cor))

  return(list(cor = acc$cor[1], a = acc$a[1], b = acc$b[1]))

}



#' Time order calculation for accuracy
#'
#' @param phe_ab phenotype data.frame
#' @param M1 pedigree kinship matrix
#' @param M2 genomic matrix
#' @param M3 omics matrix
#' @param a parameter `\alpha`
#' @param b parameter `\beta`
#' @param trait analysis trait
#' @param RP.id reference population id
#' @param VP.id validation population id
#' @param stas.phe.col target trait when calculate the indicator, default as @inheritParams trait
#' @param stas.type indicator when optimize the fusion matrix
#' @param stas.fn custom calculated indicators
#' @param ... passed parameters for other function
#'
#' @return list
#'

tv.cal_ab <- function(
    phe_ab, trait, M1, M2, M3, cov, a, b, RP.id = NULL, VP.id = NULL, stas.phe.col = NULL, stas.type, stas.fn, ...
)
{
  suppressMessages(require(rrBLUP))

  if(is.null(RP.id) | is.null(VP.id)) stop(" In time-order validation, should provide both Re(RP) and VP ids!!! ")

  colnames(phe_ab)[1] <- "id"

  kk <- FSBLUP.CombineMat(phe_c = phe_ab, M1 = M1, M2 = M2, M3 = M3, a = a, b = b)

  kk <- kk[as.character(phe_ab$id), as.character(phe_ab$id)]

  results <- data.frame()
  phe_ab$xpar <- NA
  phe_ab$xpar[phe_ab$id %in% RP.id] <- "train"
  phe_ab$xpar[phe_ab$id %in% VP.id] <- "test"
  phe_ab$y0 <- pull(phe_ab, trait)
  phe_ab$yNA0 <- pull(phe_ab, trait)
  xnas <- phe_ab$xpar == "test"
  phe_ab$yNA0[xnas] <- NA

  res_k <- tryCatch(mixed.solve(phe_ab$yNA0, X = cov, K = kk), error = {kk=adj_pos(kk);mixed.solve(phe_ab$yNA0, X = cov, K = kk)})
  res_k <- res_k$u[xnas]

  if(is.null(stas.phe.col)) x1oo00oo <- phe_ab$y0[xnas]
  if(!is.null(stas.phe.col)) x1oo00oo <- pull(phe_ab, stas.phe.col)[xnas]

  results <- rbind(
    results,
    data.frame(
      a = a,
      b = b,
      pearson = FSBLUP.stas.cal(x1oo00oo, res_k, type = stas.type, fn = stas.fn)
      #pearson = cor(phe_ab$y0[xnas], res_k, use = "pairwise.complete.obs")
    )
  )

  acc <- results %>%
    group_by(a, b) %>%
    summarise(cor = mean(pearson, na.rm = T), .groups = "drop") %>%
    ungroup() %>%
    arrange(desc(cor))

  return(list(cor = acc$cor[1], a = acc$a[1], b = acc$b[1]))

}



#' Perform cross-validation in reference population to find best parameters
#'
#' @param phe_cv pheno data.frame
#' @param trait_col analysis trait col number or name
#' @param M1 pedigree kinship matrix
#' @param M2 genomic matrix
#' @param M3 omics matrix
#' @param crv.num fold number of validation
#' @param crv.rep.num repeat times in each fold, also as random seed
#' @param gs.point_num grid search point number, default 25
#' @param bi.max_iter bisection maximum iterations
#' @param bi.threshold bisection stop  threshold between two iterations
#' @param stas.phe.col target trait when calculate the indicator, default as @inheritParams trait_Col
#' @param stas.type indicator when optimize the fusion matrix
#' @param stas.fn custom calculated indicators
#'
#' @return list
#'

FSBLUP.CrossV <- function(
    phe_cv, M1, M2, M3, trait_col, cov, crv.num, crv.rep.num, gs.point_num, bi.max_iter, bi.threshold, stas.phe.col, stas.type, stas.fn
)
{

  wind <- Sys.info()[['sysname']] == 'Windows'
  linux <- Sys.info()[['sysname']] == 'Linux'
  mac <- (!linux) & (!wind)

  if(wind)	cpu <- 1

  if(is.numeric(trait_col)) y <- colnames(phe_cv)[trait_col]
  if(is.character(trait_col)) y <- trait_col

  phe_cv <- filter(phe_cv, !is.na(y))
  colnames(phe_cv)[1] <- "id"

  cat("Start grid search procedure ... \n")
  find_ab <- FSBLUP.GridSearch(fn = cv.cal_ab, x = c(0,1),y = c(0,1), point_num = gs.point_num,
                                   phe_ab = phe_cv, trait = y, M1 = M1, M2 = M2, M3 = M3, cov = cov,
                                  crv.num = crv.num, crv.rep.num = crv.rep.num,
                                  stas.phe.col = stas.phe.col, stas.type = stas.type, stas.fn = stas.fn)
  cat("Start bisection procedure ... \n")
  find_ab <- FSBLUP.Bisection(find_ab, fn = cv.cal_ab, max_iter = bi.max_iter, threshold = bi.threshold,
                                 phe_ab = phe_cv, trait = y, M1 = M1, M2 = M2, M3 = M3, cov = cov,
                                 crv.num = crv.num, crv.rep.num = crv.rep.num,
                                 stas.phe.col = stas.phe.col, stas.type = stas.type, stas.fn = stas.fn)

  return(find_ab)
}


#' Perform time order validation in reference population to find best parameters
#'
#' @param phe_cv pheno data.frame
#' @param trait_col analysis trait col number or name
#' @param M1 pedigree kinship matrix
#' @param M2 genomic matrix
#' @param M3 omics matrix
#' @param gs.point_num grid search point number, default 25
#' @param bi.max_iter bisection maximum iterations
#' @param bi.threshold bisection stop  threshold between two iterations
#' @param stas.phe.col target trait when calculate the indicator, default as @inheritParams trait_col
#' @param stas.type indicator when optimize the fusion matrix
#' @param stas.fn custom calculated indicators
#' @param train.id reference population id
#' @param test.id validation population id
#'
#' @return list
#'

FSBLUP.TimeV <- function(
    phe_cv, M1, M2, M3, trait_col, cov, train.id, test.id, gs.point_num, bi.max_iter, bi.threshold, stas.phe.col, stas.type, stas.fn
)
{

  wind <- Sys.info()[['sysname']] == 'Windows'
  linux <- Sys.info()[['sysname']] == 'Linux'
  mac <- (!linux) & (!wind)

  if(wind)	cpu <- 1

  if(is.numeric(trait_col)) y <- colnames(phe_cv)[trait_col]
  if(is.character(trait_col)) y <- trait_col

  phe_cv <- filter(phe_cv, !is.na(y))
  colnames(phe_cv)[1] <- "id"

  cat("Start grid search with ", gs.point_num, " points ... ")
  find_ab <- FSBLUP.GridSearch(fn = tv.cal_ab, x = c(0,1),y = c(0,1), point_num = gs.point_num,
                                  phe_ab = phe_cv, trait = y, M1 = M1, M2 = M2, M3 = M3, cov = cov,
                                  train.id = train.id, test.id = test.id,
                                  stas.phe.col = stas.phe.col, stas.type = stas.type, stas.fn = stas.fn)
  cat("Start bisection procedures ... ")
  find_ab <- FSBLUP.Bisection(find_ab, fn = tv.cal_ab, max_iter = bi.max_iter, threshold = bi.threshold,
                                 phe_ab = phe_cv, trait = y, M1 = M1, M2 = M2, M3 = M3, cov = cov,
                                 train.id = train.id, test.id = test.id,
                                 stas.phe.col = stas.phe.col, stas.type = stas.type, stas.fn = stas.fn)
  return(find_ab)
}



# calculate matrix

FSBLUP.Amatrix <- AGHmatrix::Amatrix

FSBLUP.Gmatrix <- AGHmatrix::Gmatrix

#' Calculate omics similarity matrix
#'
#' @param data omics data
#' @param missing.value missing value in data, will substitute to 0
#' @param FPKM.qc if omics data is transcriptomic fprk format data, will filter transcripts which fpkm < 0.1 in more than 95% individuals
#'
#' @return matrix
#'

FSBLUP.Omatrix <- function(data = NULL, missing.value = NA, FPKM.qc = FALSE) {
  cat("\nStart calculating omic-drived similarity matrix .. ")
  start_time <- Sys.time()
  if(is.character(data[, 1])) {
    rownames(data) <- data[, 1]
    data <- data[, -1]
  }
  if(!is.numeric(data[,1])){
    stop("\n\n Detected non-numeric format value, please check your omic data!!! \n\n")
  }
  data <- as.matrix(data)
  if(!is.na(missing.value)) {
    data[data==missing.value] <- NA
  }
  data[is.na(data)] <- 0
  if(FPKM.qc) {
    data <- data[, -which((colSums(data < 0.1) / nrow(data)) >= 0.95)]
  }
  Omat <- scale(as.matrix(data))
  Omat <- tcrossprod(Omat)/ncol(Omat)
  cat("\nCompleted! Time = ", difftime(Sys.time(), start_time, units = "mins")," minutes \n")
  return(Omat)
}



#' Fusion similarity matrix calculation
#'
#' @param phe pheno data.frame
#' @param M1 pedigree kinship matrix
#' @param M2 genomic matrix
#' @param M3 omics matrix
#' @param pedi pedigree data.frame, id|sire|dam
#' @param snp snp data.frame, individuals(row) * snps(col)
#' @param omic omic data.frame, individuals(row) * omic features(col)
#' @param trait_col analysis trait in the `phe`, col number or name
#' @param po.crv.num fold number of cross validation in parameter optimization procedure
#' @param po.crv.rep.num repeat times of cross validation in parameter optimization procedure
#' @param po.year year variable col of time order validation in parameter optimization procedure
#' @param po.ngen ngen variable col of time order validation in parameter optimization procedure
#' @param po.gs.point_num grid search point number, default 25
#' @param po.bi.max_iter bisection maximum iterations
#' @param po.bi.threshold bisection stop  threshold between two iterations
#' @param stas.phe.col target trait when calculate the indicator, default as @inheritParams trait_col
#' @param stas.type indicator when optimize the fusion matrix
#' @param stas.fn custom calculated indicators
#' @param return.matrix default `T` return fusion matrix, `F` return predicted EBVs
#' @param FPKM.qc if omic data is transcriptomics FPKM format data, will delete transcript <0.1 in more than 95% individuals.
#' @param fix.col column number of fixed effects in phenotype, e.g. c(2, 3)
#'
#' @return matrix
#' @export
#'

FSBLUP <- function(
    phe = NULL, M1 = NULL, M2 = NULL, M3 = NULL, pedi = NULL, snp = NULL, omic = NULL, FPKM.qc = TRUE, fix.col = NULL,
    trait_col = NULL, po.crv.num = 5, po.crv.rep.num=2, po.year = NULL, po.ngen = NULL,
    po.gs.point_num = 25, po.bi.max_iter = 10, po.bi.threshold = 1e-4,
    stas.type = "cor", stas.fn = NULL, stas.phe.col = NULL, return.matrix = FALSE
)
{
  FSBLUP.version()

  start_time = Sys.time()


  if(is.null(M1) & is.null(M2) & is.null(M3)) {
    if(is.null(pedi) & is.null(snp) & is.null(omic)) {
      stop("\n\n Unknow infomation input, please check parameters !!! \n\n")
    }
    cat("\n\nDetected original data input, start pedigree, genomic, omics matrices calculating !!! \n\n")

    if(is.null(pedi)) stop("\n\n There's no pedigree file input, please ckeck input parameters [pedi] !!! \n\n")
    if(is.null(snp)) stop("\n\n There's no pedigree file input, please ckeck input parameters [snp] !!! \n\n")
    if(is.null(omic)) stop("\n\n There's no pedigree file input, please ckeck input parameters [omic] !!! \n\n")

    cat("\nStart calculating A matrix ..\n")
    M1 <- FSBLUP.Amatrix(pedi)
    cat("\nStart calculating G matrix ..\n")
    M2 <- FSBLUP.Gmatrix(snp)
    M3 <- FSBLUP.Omatrix(omic, FPKM.qc = FPKM.qc)
  }
  phe <- as.data.frame(phe)

  VPindex <- is.na(pull(phe, trait_col))
  VPid <- as.character(pull(phe, 1))[VPindex]
  RFid <- as.character(pull(phe, 1))[!VPindex]

  phe_cv <- phe %>% filter(pull(., 1) %in% RFid)

  M1_cv <- M1[rownames(M1) %in% RFid, colnames(M1) %in% RFid]
  M2_cv <- M2[rownames(M2) %in% RFid, colnames(M2) %in% RFid]
  M3_cv <- M3[rownames(M3) %in% RFid, colnames(M3) %in% RFid]

  if(!is.null(fix.col)) {
    cov <- phe[, fix.col]
    cov <- matrix(cov)
  }

  if(!is.null(po.crv.num))
  {
    cat("\nStart optimizing Fusion Similarity matrix ..\n")
    opt <- FSBLUP.CrossV(phe_cv = phe_cv, trait_col = trait_col, M1 = M1_cv, M2 = M2_cv, M3 = M3_cv, cov = cov,
                            crv.num = po.crv.num, crv.rep.num = po.crv.rep.num, gs.point_num = po.gs.point_num,
                            bi.max_iter = po.bi.max_iter, bi.threshold = po.bi.threshold,
                            stas.type = stas.type, stas.fn = stas.fn, stas.phe.col = stas.phe.col)
  } else if((is.null(po.year) | is.null(po.ngen)) & is.null(po.crv.num))
  {
    if(!is.null(po.year))
    {
      if(is.na(as.numeric(po.year))) stop("\n\n Parameter error: 'po.year' should be a number (e.g. 2024) \n\n")
      if(! "year" %in% colnames(phe)) stop("\n\n Can't find 'year' col in phenotype data !!! \n\n")
      train.id <- phe[,1][phe$year <= as.numeric(po.year)]
      test.id <- phe[,1][phe$year > as.numeric(po.year)]
      if(length(test.id) == 0) stop("\n\n In parameter optimization procedure, after ", po.year," year, all individuals belong to training population, please check input parameters (e.g. po.year) !!! \n\n")

    }
    if(!is.null(po.ngen))
    {
      if(is.na(as.numeric(po.ngen))) stop("\n\n Parameter error: 'po.ngen' should be a number (e.g. 5, 6, 7) \n\n")
      if(! "ngen" %in% colnames(phe)) stop("\n\n Can't find 'ngen' col in phenotype data !!! \n\n")
      train.id <- phe[,1][phe$ngen <= as.numeric(po.ngen)]
      test.id <- phe[,1][phe$ngen > as.numeric(po.ngen)]
      if(length(test.id) == 0) stop("\n\n In parameter optimization procedure, after ", po.ngen," ngen, all individuals belong to training population, please check input parameters (e.g. po.ngen) !!! \n\n")
    }
    if(length(train.id) == 0 | length(test.id) == 0) stop("\n\n In parameter optimization procedure, train or test population contains 0 individuals, please check input parameters (e.g. po.year, po.ngen) !!! \n\n")

    cat("Start optimize Fusion Similarity matrix \n")
    opt <- FSBLUP.TimeV(phe_cv = phe_cv, trait_col = trait_col, M1 = M1_cv, M2 = M2_cv, M3 = M3_cv, cov = cov,
                           train.id = train.id, test.id = test.id, gs.point_num = po.gs.point_num,
                           bi.max_iter = po.bi.max_iter, bi.threshold = po.bi.threshold,
                           stas.type = stas.type, stas.fn = stas.fn, stas.phe.col = stas.phe.col)
  } else {
    stop("\n\n Parameters shouldn't be NULL at same time for: 'po.year', 'po.ngen' and 'po.crv.num' !!! \n\n")
  }

  k <- FSBLUP.CombineMat(phe_c = phe, M1 = M1, M2 = M2, M3 = M3, a = opt$optimal_a, b = opt$optimal_b)
  k <- k[as.character(phe[,1]), as.character(phe[,1])]

  cat("Optimized procedure finish with ", difftime(Sys.time(), start_time, units = "mins")," minutes \n")

  if(return.matrix) {
    return(k)
  } else {
    cat("Predicating ... \n\n")
    res_k <- tryCatch(mixed.solve(phe[, trait_col], K = k), error = {k=adj_pos(k);mixed.solve(phe[, trait_col], K = k)})
    return(res_k$u)
  }
}




