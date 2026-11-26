# its.R — 中斷時間序列（ITS）分段迴歸分析函數
# 所有分析數字都必須來自這支檔案的函數。
# 依賴：nlme、sandwich、lmtest、ggplot2

suppressPackageStartupMessages({
  library(nlme)
  library(sandwich)
  library(lmtest)
  library(ggplot2)
})

# ── 0. 設定 ─────────────────────────────────────────────────────
its_configure <- function() {
  options(digits = 7, contrasts = c("contr.treatment", "contr.poly"))
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(20261118)
  invisible(TRUE)
}

.issue <- function(level, msg) data.frame(level = level, message = msg)

# 解析時間：Date、"YYYY-MM"、"YYYY-MM-DD"、數值（年或小數年）
.parse_time <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (is.numeric(x)) return(x)
  x <- trimws(as.character(x))
  if (all(grepl("^\\d{4}-\\d{1,2}$", x))) return(as.Date(paste0(x, "-01")))
  if (all(grepl("^\\d{4}/\\d{1,2}$", x))) return(as.Date(paste0(gsub("/", "-", x), "-01")))
  d <- suppressWarnings(as.Date(x))
  if (!anyNA(d)) return(d)
  n <- suppressWarnings(as.numeric(x))
  if (!anyNA(n)) return(n)
  stop("無法解析時間欄位；請用 YYYY-MM、YYYY-MM-DD 或數值年份")
}

.month_index <- function(d) {
  as.integer(format(d, "%Y")) * 12L + as.integer(format(d, "%m"))
}

# ── 1. 準備資料與檢查 ─────────────────────────────────────────
# t0：介入「開始」的時間點（介入後第一個觀測），格式同時間欄位
# transition：介入後排除幾期當過渡期（0 = 立即生效）
its_prepare <- function(data, time, y, t0, freq, transition = 0,
                        covariates = NULL) {
  issues <- .issue(character(), character())
  add <- function(level, msg) issues <<- rbind(issues, .issue(level, msg))

  need <- c(time, y, covariates)
  miss_col <- setdiff(need, names(data))
  if (length(miss_col)) stop("找不到欄位：", paste(miss_col, collapse = ", "))
  if (!freq %in% c(1, 4, 7, 12, 52)) stop("freq 只支援 1、4、7、12、52")

  tm <- .parse_time(data[[time]])
  t0p <- .parse_time(t0)
  d <- data.frame(row_id = seq_len(nrow(data)), time = tm,
                  y = suppressWarnings(as.numeric(data[[y]])))
  for (v in covariates) d[[v]] <- suppressWarnings(as.numeric(data[[v]]))

  if (anyNA(d$time)) add("stop", "時間欄位有缺值或無法解析的值")
  if (anyDuplicated(d$time)) add("stop", "時間欄位有重複的時間點")
  if (is.unsorted(d$time)) {
    d <- d[order(d$time), ]
    add("note", "資料原本沒有依時間排序，已依時間重新排序（保留原始 row_id）")
  }

  # 間隔是否相等
  dt <- if (inherits(d$time, "Date")) {
    if (freq %in% c(1, 4, 12)) {
      step <- c(`1` = 12L, `4` = 3L, `12` = 1L)[as.character(freq)]
      diff(.month_index(d$time)) == step
    } else {
      step <- if (freq == 52) 7 else 1
      as.numeric(diff(d$time)) == step
    }
  } else {
    dd <- diff(d$time)
    abs(dd - dd[1]) < 1e-6 * max(1, abs(dd[1]))
  }
  if (length(dt) && !all(dt)) {
    add("stop", sprintf("時間間隔不相等（%d 處不連續），請先補齊或確認頻率",
                        sum(!dt)))
  }

  if (anyNA(d$y)) {
    add("stop", sprintf("結果變數有 %d 個缺值或非數值；本 skill 不自動插補",
                        sum(is.na(d$y))))
  }
  for (v in covariates) if (anyNA(d[[v]])) {
    add("stop", sprintf("共變數 %s 有缺值", v))
  }

  # 介入點
  T0 <- NA_integer_
  if (inherits(d$time, "Date") != inherits(t0p, "Date")) {
    add("stop", "介入時間的格式跟時間欄位不一致")
  } else {
    idx <- which(d$time >= t0p - 1e-9)
    if (!length(idx) || t0p <= min(d$time)) {
      add("stop", "介入時間不在序列範圍內（必須晚於第一個觀測、不晚於最後一個觀測）")
    } else {
      T0 <- idx[1]
      if (abs(as.numeric(d$time[T0]) - as.numeric(t0p)) > 1e-6) {
        add("note", sprintf("序列中沒有剛好等於介入時間的點，以其後第一個觀測（%s）當介入起點",
                            format(d$time[T0])))
      }
    }
  }

  d$T <- seq_len(nrow(d))
  if (!is.na(T0)) {
    d$X <- as.integer(d$T >= T0)
    d$Tpost <- pmax(d$T - T0, 0L)
  } else {
    d$X <- 0L
    d$Tpost <- 0L
  }
  d$excluded <- d$X == 1 & d$Tpost < transition

  # 季節位置
  if (freq > 1) {
    pos <- if (inherits(d$time, "Date") && freq == 12) {
      as.integer(format(d$time, "%m"))
    } else if (inherits(d$time, "Date") && freq == 4) {
      (as.integer(format(d$time, "%m")) - 1L) %/% 3L + 1L
    } else if (inherits(d$time, "Date") && freq == 7) {
      as.integer(format(d$time, "%u"))
    } else if (is.numeric(d$time) && freq %in% c(4, 12) &&
               any(abs(d$time - round(d$time)) > 1e-9)) {
      as.integer(round((d$time - floor(d$time)) * freq)) + 1L
    } else {
      add("note", "季節位置依資料順序推算（第一筆 = 第 1 期）")
      ((d$T - 1L) %% freq) + 1L
    }
    d$season_pos <- pos
    d$season <- factor(pos, levels = sort(unique(pos)))
  }

  n_pre <- sum(d$X == 0)
  n_post <- sum(d$X == 1 & !d$excluded)
  if (!is.na(T0)) {
    if (n_pre < 3 || n_post < 3) {
      add("stop", sprintf("介入前 %d 點、介入後 %d 點（扣除過渡期）；任一段少於 3 點無法估計斜率",
                          n_pre, n_post))
    } else if (n_pre < 12 || n_post < 12) {
      add("warn", sprintf("介入前 %d 點、介入後 %d 點；任一段少於 12 點，檢定力有限",
                          n_pre, n_post))
    }
    if (freq > 1 && n_pre < freq) {
      add("warn", "介入前資料不到一個完整季節週期，季節項估計不穩")
    }
  }
  if (transition > 0) {
    add("note", sprintf("排除介入後前 %d 期當過渡期（不納入配適，圖上以空心點表示）",
                        transition))
  }

  is_count <- all(d$y >= 0, na.rm = TRUE) &&
    all(abs(d$y - round(d$y)) < 1e-9, na.rm = TRUE)
  if (is_count && mean(d$y, na.rm = TRUE) < 10) {
    add("warn", "結果變數是平均小於 10 的計數，線性模型近似不佳，建議以 quasi-Poisson 為主")
  }

  info <- list(time_col = time, y_col = y, t0_input = t0,
               t0_time = if (!is.na(T0)) d$time[T0] else NA,
               T0 = T0, freq = freq, transition = transition,
               covariates = covariates, n = nrow(d),
               n_pre = n_pre, n_post = n_post,
               n_excluded = sum(d$excluded), is_count = is_count,
               y_mean = mean(d$y, na.rm = TRUE))
  status <- if (any(issues$level == "stop")) "stop" else "ok"
  structure(list(data = d, info = info, issues = issues, status = status),
            class = "its_prep")
}

print.its_prep <- function(x, ...) {
  i <- x$info
  cat(sprintf("n = %d｜介入起點 T0 = %s（第 %s 期）｜介入前 %d 點、介入後 %d 點｜頻率 %d\n",
              i$n, format(i$t0_time), i$T0, i$n_pre, i$n_post, i$freq))
  cat("狀態：", if (x$status == "ok") "可以繼續" else "必須停止", "\n")
  if (nrow(x$issues)) {
    lab <- c(stop = "⛔ 停止", warn = "⚠️ 警告", note = "ℹ️ 說明")
    for (k in seq_len(nrow(x$issues))) {
      cat(sprintf("  %s：%s\n", lab[x$issues$level[k]], x$issues$message[k]))
    }
  }
  invisible(x)
}

.check_ok <- function(prep) {
  if (!inherits(prep, "its_prep")) stop("請先用 its_prepare() 準備資料")
  if (prep$status != "ok") stop("資料檢查未通過，不能繼續：\n",
                                paste(prep$issues$message[prep$issues$level == "stop"],
                                      collapse = "\n"))
}

# ── 2. 探索 ─────────────────────────────────────────────────────
its_explore <- function(prep, out_dir = ".", y_label = NULL) {
  .check_ok(prep)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  d <- prep$data
  i <- prep$info
  y_label <- y_label %||% i$y_col
  tx <- as.numeric(d$time)
  pre <- d[d$X == 0, ]

  # 時序圖 + 介入前 loess
  g1 <- ggplot(d, aes(time, y)) +
    geom_line(colour = "grey70") +
    geom_point(aes(shape = ifelse(excluded, "過渡期（排除）", "納入配適")),
               size = 1.3, colour = "grey30") +
    scale_shape_manual(values = c("納入配適" = 16, "過渡期（排除）" = 1),
                       name = NULL) +
    geom_smooth(data = pre, method = "loess", formula = y ~ x, se = FALSE,
                colour = "#1f6feb", linewidth = 0.8) +
    geom_vline(xintercept = d$time[i$T0], linetype = "dashed",
               colour = "#d1242f") +
    labs(x = i$time_col, y = y_label,
         title = "時序圖（紅虛線 = 介入起點；藍線 = 介入前 loess 趨勢）") +
    theme_minimal() + theme(legend.position = "bottom")
  f1 <- file.path(out_dir, "explore_timeseries.png")
  ggsave(f1, g1, width = 8, height = 4.2, dpi = 120)

  files <- c(timeseries = f1)
  season_tab <- NULL
  if (i$freq > 1) {
    season_tab <- aggregate(y ~ season, data = pre, FUN = mean)
    g2 <- ggplot(pre, aes(season, y)) +
      geom_boxplot(fill = "#dbe9ff", outlier.size = 0.8) +
      labs(x = "季節位置（月／季／星期）", y = y_label,
           title = "季節圖（只用介入前資料，避免混入介入效果）") +
      theme_minimal()
    f2 <- file.path(out_dir, "explore_season.png")
    ggsave(f2, g2, width = 8, height = 3.8, dpi = 120)
    files <- c(files, season = f2)
  }

  # 介入前趨勢線性檢查：二次項 F 檢定（有季節就控制季節）
  f_lin <- if (i$freq > 1) y ~ T + season else y ~ T
  m_lin <- lm(f_lin, data = pre)
  m_quad <- update(m_lin, . ~ . + I(T^2))
  p_curve <- anova(m_lin, m_quad)$`Pr(>F)`[2]

  # 季節性 F 檢定（介入前）
  p_season <- NA
  if (i$freq > 1) {
    p_season <- anova(lm(y ~ T, data = pre), lm(y ~ T + season, data = pre))$`Pr(>F)`[2]
  }

  suggestion <- list(
    season = if (i$freq == 1) "none" else if (i$freq %in% c(4, 12)) "dummy" else "fourier",
    family = if (i$is_count && i$y_mean < 10) "quasipoisson" else "gaussian"
  )
  res <- list(files = files, n_pre = i$n_pre, n_post = i$n_post,
              p_pre_curvature = p_curve, p_season = p_season,
              season_means = season_tab, suggestion = suggestion)
  class(res) <- "its_explore"
  res
}

print.its_explore <- function(x, ...) {
  cat(sprintf("介入前 %d 點、介入後 %d 點\n", x$n_pre, x$n_post))
  cat(sprintf("介入前趨勢的二次項檢定 p = %.4g（< .05 表示介入前趨勢明顯彎曲，反事實要小心）\n",
              x$p_pre_curvature))
  if (!is.na(x$p_season)) {
    cat(sprintf("介入前季節性 F 檢定 p = %.4g\n", x$p_season))
  }
  cat(sprintf("建議設定：季節項 = %s；模型家族 = %s\n",
              x$suggestion$season, x$suggestion$family))
  cat("圖檔：", paste(x$files, collapse = "、"), "\n")
  invisible(x)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ── 3. 模型設定 ─────────────────────────────────────────────────
its_spec <- function(prep, season = c("dummy", "fourier", "none"), K = 2,
                     covariates = NULL) {
  .check_ok(prep)
  season <- match.arg(season)
  if (prep$info$freq == 1 && season != "none") {
    stop("年資料沒有季節性，season 請用 \"none\"")
  }
  rhs <- c("T", "X", "Tpost")
  if (season == "dummy") rhs <- c(rhs, "season")
  if (season == "fourier") {
    K <- min(K, floor(prep$info$freq / 2))
    rhs <- c(rhs, paste0(c("S", "C"), rep(seq_len(K), each = 2)))
  }
  covariates <- covariates %||% character()
  bad <- setdiff(covariates, names(prep$data))
  if (length(bad)) stop("共變數不在 its_prepare() 的 covariates 裡：", paste(bad, collapse = ", "))
  rhs <- c(rhs, covariates)
  structure(list(formula = as.formula(paste("y ~", paste(rhs, collapse = " + "))),
                 season = season, K = K, covariates = covariates),
            class = "its_spec")
}

# 加上 Fourier 欄位（需要時）
.design_data <- function(prep, spec) {
  d <- prep$data
  if (spec$season == "fourier") {
    for (k in seq_len(spec$K)) {
      d[[paste0("S", k)]] <- sin(2 * pi * k * d$season_pos / prep$info$freq)
      d[[paste0("C", k)]] <- cos(2 * pi * k * d$season_pos / prep$info$freq)
    }
  }
  d
}

# ── 4. 配適 ─────────────────────────────────────────────────────
# method = "ols"：lm；"gls"：nlme::gls + AR(p) 誤差（REML）
its_fit <- function(prep, spec, method = c("ols", "gls"), p = 0,
                    approved = FALSE, estimation = "REML") {
  if (!isTRUE(approved)) stop("使用者尚未確認模型設定，不能配適")
  .check_ok(prep)
  method <- match.arg(method)
  d <- .design_data(prep, spec)
  dfit <- d[!d$excluded, ]
  f <- spec$formula
  if (method == "ols") {
    m <- lm(f, data = dfit)
  } else {
    cor <- if (p > 0) corARMA(p = p, q = 0, form = ~ T) else NULL
    m <- eval(bquote(gls(.(f), data = dfit, correlation = cor,
                         method = .(estimation))))
  }
  attr(m, "its") <- list(method = method, p = if (method == "ols") 0 else p,
                         spec = spec, info = prep$info)
  m
}

.its_label <- function(m) {
  a <- attr(m, "its")
  if (a$method == "ols") "OLS" else sprintf("GLS AR(%d)", a$p)
}

# ── 5. 殘差診斷 ─────────────────────────────────────────────────
its_diagnose <- function(fit, prep, alpha = 0.05, out_dir = ".", tag = NULL) {
  a <- attr(fit, "its")
  n <- nobs(fit)
  freq <- prep$info$freq
  lb_lag <- min(if (freq > 1) 2 * freq else 10, floor(n / 5))
  bg_order <- if (freq > 1) freq else 4
  out <- data.frame(test = character(), statistic = numeric(),
                    p_value = numeric(), note = character())
  if (a$method == "ols") {
    r <- residuals(fit)
    dw <- dwtest(fit, alternative = "two.sided")
    bg <- bgtest(fit, order = bg_order)
    out <- rbind(out,
      data.frame(test = "Durbin–Watson（一階，雙尾）", statistic = unname(dw$statistic),
                 p_value = dw$p.value, note = "DW ≈ 2 表示沒有一階自相關"),
      data.frame(test = sprintf("Breusch–Godfrey（到 %d 階）", bg_order),
                 statistic = unname(bg$statistic), p_value = bg$p.value,
                 note = "檢定高階自相關"))
  } else {
    r <- residuals(fit, type = "normalized")
  }
  lb <- Box.test(r, lag = lb_lag, type = "Ljung-Box")
  out <- rbind(out, data.frame(
    test = sprintf("Ljung–Box（%d 個落後期）%s", lb_lag,
                   if (a$method == "gls") "，標準化殘差" else ""),
    statistic = unname(lb$statistic), p_value = lb$p.value,
    note = "綜合檢定多個落後期"))
  sw <- shapiro.test(r)
  out <- rbind(out, data.frame(test = "Shapiro–Wilk（常態性）",
                               statistic = unname(sw$statistic),
                               p_value = sw$p.value, note = "參考用；n 大時很敏感"))
  ac <- acf(r, lag.max = min(36, n - 1), plot = FALSE)$acf[-1]
  autocorr <- any(out$p_value[!grepl("Shapiro", out$test)] < alpha)

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  tag <- tag %||% gsub("[^A-Za-z0-9]", "", .its_label(fit))
  fpng <- file.path(out_dir, sprintf("diag_%s.png", tag))
  png(fpng, width = 1100, height = 760, res = 120)
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  plot(r, type = "b", pch = 16, cex = 0.5, xlab = "時間順序", ylab = "殘差",
       main = paste(.its_label(fit), "殘差 vs. 時間"))
  abline(h = 0, lty = 2)
  qqnorm(r, main = "常態 Q-Q 圖", pch = 16, cex = 0.5); qqline(r)
  acf(r, lag.max = min(36, n - 1), main = "殘差 ACF")
  pacf(r, lag.max = min(36, n - 1), main = "殘差 PACF")
  par(op)
  dev.off()

  res <- list(model = .its_label(fit), tests = out, autocorrelation = autocorr,
              acf_lag1 = ac[1], acf_seasonal = if (freq > 1 && length(ac) >= freq) ac[freq] else NA,
              n_acf_outside = sum(abs(ac) > qnorm(0.975) / sqrt(n)),
              n_acf_lags = length(ac), plot = fpng, alpha = alpha)
  class(res) <- "its_diag"
  res
}

print.its_diag <- function(x, ...) {
  cat("──", x$model, "殘差診斷 ──\n")
  tt <- x$tests
  tt$statistic <- round(tt$statistic, 3)
  tt$p_value <- signif(tt$p_value, 3)
  print(tt, row.names = FALSE)
  cat(sprintf("ACF lag1 = %.3f%s；%d / %d 個落後期超出 95%% 界線\n",
              x$acf_lag1,
              if (!is.na(x$acf_seasonal)) sprintf("，季節落後期 = %.3f", x$acf_seasonal) else "",
              x$n_acf_outside, x$n_acf_lags))
  cat("判定：", if (x$autocorrelation) "有自相關（BG 或 Ljung–Box p < α）" else "沒有足夠證據顯示自相關", "\n")
  cat("圖檔：", x$plot, "\n")
  invisible(x)
}

# ── 6. 選 AR 階數 ───────────────────────────────────────────────
its_select_ar <- function(prep, spec, p_max = 4) {
  .check_ok(prep)
  tab <- data.frame(p = 0:p_max, AIC = NA_real_, BIC = NA_real_, converged = FALSE)
  for (p in 0:p_max) {
    m <- tryCatch(its_fit(prep, spec, "gls", p, approved = TRUE, estimation = "ML"),
                  error = function(e) NULL)
    if (!is.null(m)) {
      tab$AIC[p + 1] <- AIC(m)
      tab$BIC[p + 1] <- BIC(m)
      tab$converged[p + 1] <- TRUE
    }
  }
  tab$dAIC <- tab$AIC - min(tab$AIC, na.rm = TRUE)
  best <- tab$p[which.min(tab$AIC)]
  structure(list(table = tab, best_p = best), class = "its_ar")
}

print.its_ar <- function(x, ...) {
  cat("GLS AR(p) 比較（ML 估計，AIC 越小越好）：\n")
  t <- x$table
  t[c("AIC", "BIC", "dAIC")] <- round(t[c("AIC", "BIC", "dAIC")], 1)
  print(t, row.names = FALSE)
  cat(sprintf("→ 建議 p = %d\n", x$best_p))
  invisible(x)
}

# ── 7. 效果估計 ─────────────────────────────────────────────────
.model_matrix_all <- function(fit, prep) {
  a <- attr(fit, "its")
  d <- .design_data(prep, a$spec)
  mm <- model.matrix(delete.response(terms(a$spec$formula)), data = d)
  cf <- coef(fit)
  miss <- setdiff(names(cf), colnames(mm))
  if (length(miss)) stop("設計矩陣欄位不一致：", paste(miss, collapse = ", "))
  mm[, names(cf), drop = FALSE]
}

# vcov_type：model（模型本身）或 newey_west（只限 OLS）
its_vcov <- function(fit, type = c("model", "newey_west"), nw_lag = NULL) {
  type <- match.arg(type)
  if (type == "model") return(vcov(fit))
  if (!inherits(fit, "lm")) stop("Newey–West 只用在 OLS")
  n <- nobs(fit)
  nw_lag <- nw_lag %||% floor(4 * (n / 100)^(2 / 9))
  v <- NeweyWest(fit, lag = nw_lag, prewhite = FALSE, adjust = TRUE)
  attr(v, "nw_lag") <- nw_lag
  v
}

its_effects <- function(fit, prep, k = NULL, level = 0.95, vcov_type = "model",
                        nw_lag = NULL) {
  a <- attr(fit, "its")
  b <- coef(fit)
  V <- its_vcov(fit, vcov_type, nw_lag)
  dfr <- nobs(fit) - length(b)
  q <- qt(1 - (1 - level) / 2, dfr)
  d <- prep$data
  kmax <- max(d$Tpost)
  freq <- prep$info$freq
  k <- k %||% unique(c(0, if (freq > 1 && freq <= kmax) freq, kmax))
  k <- k[k <= kmax]

  vec <- function(...) {
    v <- setNames(numeric(length(b)), names(b))
    z <- list(...)
    v[names(z)] <- unlist(z)
    v
  }
  rows <- list(
    list("介入前斜率 β1（每期）", vec(T = 1)),
    list("水準變化 β2（介入當下）", vec(X = 1)),
    list("斜率變化 β3（每期）", vec(Tpost = 1)),
    list("介入後斜率 β1+β3（每期）", vec(T = 1, Tpost = 1))
  )
  for (kk in k) {
    rows[[length(rows) + 1]] <- list(sprintf("介入後第 %d 期的效果 β2+%d·β3", kk, kk),
                                     vec(X = 1, Tpost = kk))
  }
  tab <- do.call(rbind, lapply(rows, function(r) {
    est <- sum(r[[2]] * b)
    se <- sqrt(drop(t(r[[2]]) %*% V %*% r[[2]]))
    data.frame(quantity = r[[1]], estimate = est, se = se,
               lower = est - q * se, upper = est + q * se,
               p_value = 2 * pt(-abs(est / se), dfr))
  }))

  # 整個介入後期間的相對效果（delta method）
  mm <- .model_matrix_all(fit, prep)
  post <- d$X == 1 & !d$excluded
  A <- colSums(mm[post, , drop = FALSE] * 0)
  A[c("X", "Tpost")] <- c(sum(post), sum(d$Tpost[post]))
  mm_cf <- mm
  mm_cf[, "X"] <- 0
  mm_cf[, "Tpost"] <- 0
  B <- colSums(mm_cf[post, , drop = FALSE])
  num <- sum(A * b)
  den <- sum(B * b)
  ratio <- num / den
  grad <- A / den - num * B / den^2
  se_r <- sqrt(drop(t(grad) %*% V %*% grad))
  rel <- data.frame(quantity = "介入後期間的平均相對變化（效果 ÷ 反事實）",
                    estimate = ratio, se = se_r,
                    lower = ratio - q * se_r, upper = ratio + q * se_r,
                    p_value = 2 * pt(-abs(ratio / se_r), dfr))

  structure(list(table = tab, relative = rel, model = .its_label(fit),
                 vcov_type = vcov_type, nw_lag = attr(V, "nw_lag"),
                 level = level, df = dfr,
                 mean_effect = num / sum(post), mean_counterfactual = den / sum(post)),
            class = "its_effects")
}

print.its_effects <- function(x, ...) {
  cat(sprintf("── 效果估計：%s%s，%d%% CI（t 分配，df = %d）──\n", x$model,
              if (x$vcov_type == "newey_west") sprintf(" + Newey–West(%d)", x$nw_lag) else "",
              round(100 * x$level), x$df))
  t <- rbind(x$table, x$relative)
  t[, 2:5] <- round(t[, 2:5], 4)
  t$p_value <- signif(t$p_value, 3)
  print(t, row.names = FALSE)
  invisible(x)
}

# ── 8. 作圖：實際值、配適線、反事實線 ───────────────────────────
its_plot <- function(fit, prep, file = NULL, y_label = NULL, title = NULL) {
  d <- prep$data
  i <- prep$info
  mm <- .model_matrix_all(fit, prep)
  b <- coef(fit)
  d$fitted <- drop(mm %*% b)
  mm_cf <- mm
  mm_cf[, "X"] <- 0
  mm_cf[, "Tpost"] <- 0
  d$counterfactual <- drop(mm_cf %*% b)
  d$counterfactual[d$X == 0] <- NA
  # 反事實線從介入前最後一點接起
  d$counterfactual[i$T0 - 1] <- d$fitted[i$T0 - 1]
  d$fit_pre <- ifelse(d$X == 0, d$fitted, NA)
  d$fit_post <- ifelse(d$X == 1, d$fitted, NA)
  d$point <- ifelse(d$excluded, "過渡期（排除）", "觀測值")

  g <- ggplot(d, aes(time)) +
    geom_point(aes(y = y, shape = point), colour = "grey40", size = 1.2) +
    scale_shape_manual(values = c("觀測值" = 16, "過渡期（排除）" = 1), name = NULL) +
    geom_line(aes(y = fit_pre, colour = "配適值"), linewidth = 0.7, na.rm = TRUE) +
    geom_line(aes(y = fit_post, colour = "配適值"), linewidth = 0.7, na.rm = TRUE) +
    geom_line(aes(y = counterfactual, colour = "反事實（沒有介入）"),
              linetype = "dashed", linewidth = 0.7, na.rm = TRUE) +
    scale_colour_manual(values = c("配適值" = "#1f6feb", "反事實（沒有介入）" = "#d1242f"),
                        name = NULL) +
    geom_vline(xintercept = d$time[i$T0], linetype = "dotted") +
    labs(x = i$time_col, y = y_label %||% i$y_col,
         title = title %||% sprintf("中斷時間序列：%s", .its_label(fit))) +
    theme_minimal() + theme(legend.position = "bottom")
  if (!is.null(file)) ggsave(file, g, width = 8.5, height = 4.5, dpi = 120)
  invisible(list(plot = g, data = d[, c("time", "y", "fitted", "counterfactual", "excluded")]))
}

# ── 9. 敏感度分析 ───────────────────────────────────────────────
# main：主要模型；which：要做的項目
its_sensitivity <- function(main, prep, which = c("newey_west", "covariates",
                                                  "transition", "pre_window", "quasipoisson"),
                            covariates = NULL, transition = NULL, pre_window = NULL,
                            data = NULL, nw_lag = NULL, level = 0.95) {
  a <- attr(main, "its")
  spec <- a$spec
  pick <- function(e, lab, note = "") {
    t <- e$table
    r2 <- t[grepl("β2（", t$quantity), ]
    r3 <- t[grepl("^斜率變化", t$quantity), ]
    data.frame(model = lab, b2 = r2$estimate, b2_se = r2$se, b2_lower = r2$lower,
               b2_upper = r2$upper, b2_p = r2$p_value,
               b3 = r3$estimate, b3_se = r3$se, b3_lower = r3$lower,
               b3_upper = r3$upper, b3_p = r3$p_value, scale = "原尺度（差值）",
               note = note)
  }
  out <- pick(its_effects(main, prep, level = level),
              paste("主要模型：", .its_label(main)))
  ols <- if (a$method == "ols") main else its_fit(prep, spec, "ols", approved = TRUE)

  if ("newey_west" %in% which) {
    e_ols <- its_effects(ols, prep, level = level)
    e_nw <- its_effects(ols, prep, level = level, vcov_type = "newey_west", nw_lag = nw_lag)
    note <- if (any(e_nw$table$se[2:3] < e_ols$table$se[2:3])) {
      "⚠️ Newey–West 的 SE 比 OLS 還小；後段點數少時 HAC 標準誤不可靠，僅供參考"
    } else ""
    if (a$method != "ols") out <- rbind(out, pick(e_ols, "OLS（未處理自相關）"))
    out <- rbind(out, pick(e_nw, sprintf("OLS + Newey–West(%d)", e_nw$nw_lag), note))
  }
  if ("covariates" %in% which && length(covariates)) {
    sp2 <- its_spec(prep, spec$season, spec$K, covariates = c(spec$covariates, covariates))
    m2 <- its_fit(prep, sp2, a$method, a$p, approved = TRUE)
    out <- rbind(out, pick(its_effects(m2, prep, level = level),
                           sprintf("主要模型 + 共變數（%s）", paste(covariates, collapse = "、"))))
  }
  if ("transition" %in% which && !is.null(transition) && !is.null(data)) {
    i <- prep$info
    p2 <- its_prepare(data, i$time_col, i$y_col, i$t0_input, i$freq,
                      transition = transition, covariates = i$covariates)
    if (p2$status == "ok") {
      m3 <- its_fit(p2, spec, a$method, a$p, approved = TRUE)
      out <- rbind(out, pick(its_effects(m3, p2, level = level),
                             sprintf("排除介入後前 %d 期（過渡期）", transition)))
    } else {
      out <- rbind(out, data.frame(model = sprintf("排除前 %d 期", transition),
                                   b2 = NA, b2_se = NA, b2_lower = NA, b2_upper = NA, b2_p = NA,
                                   b3 = NA, b3_se = NA, b3_lower = NA, b3_upper = NA, b3_p = NA,
                                   scale = "", note = "排除後點數不足，無法估計"))
    }
  }
  if ("pre_window" %in% which && !is.null(pre_window) && !is.null(data)) {
    i <- prep$info
    keep_from <- i$T0 - pre_window
    if (keep_from < 1) {
      out <- rbind(out, .na_row(sprintf("只用介入前 %d 期", pre_window), "介入前資料不足這麼多期"))
    } else {
      o <- order(.parse_time(data[[i$time_col]]))
      d2 <- data[o, , drop = FALSE][keep_from:nrow(data), , drop = FALSE]
      p2 <- its_prepare(d2, i$time_col, i$y_col, i$t0_input, i$freq,
                        transition = i$transition, covariates = i$covariates)
      if (p2$status == "ok") {
        m4 <- its_fit(p2, spec, a$method, a$p, approved = TRUE)
        out <- rbind(out, pick(its_effects(m4, p2, level = level),
                               sprintf("只用介入前最近 %d 期（從 %s 起）", pre_window,
                                       format(p2$data$time[1]))))
      } else {
        out <- rbind(out, .na_row(sprintf("只用介入前 %d 期", pre_window), "資料檢查未通過"))
      }
    }
  }
  if ("quasipoisson" %in% which) {
    if (!prep$info$is_count) {
      out <- rbind(out, data.frame(model = "quasi-Poisson", b2 = NA, b2_se = NA,
                                   b2_lower = NA, b2_upper = NA, b2_p = NA, b3 = NA,
                                   b3_se = NA, b3_lower = NA, b3_upper = NA, b3_p = NA,
                                   scale = "", note = "結果變數不是計數，略過"))
    } else {
      qp <- its_fit_quasipoisson(prep, spec, level = level, nw_lag = nw_lag)
      out <- rbind(out, qp$row)
    }
  }
  rownames(out) <- NULL
  class(out) <- c("its_sens", "data.frame")
  out
}

.na_row <- function(model, note) {
  data.frame(model = model, b2 = NA, b2_se = NA, b2_lower = NA, b2_upper = NA, b2_p = NA,
             b3 = NA, b3_se = NA, b3_lower = NA, b3_upper = NA, b3_p = NA,
             scale = "", note = note)
}

# quasi-Poisson（計數資料）：係數取指數 = 率比；SE 用 Newey–West 處理自相關
its_fit_quasipoisson <- function(prep, spec, level = 0.95, nw_lag = NULL, hac = TRUE) {
  d <- .design_data(prep, spec)
  dfit <- d[!d$excluded, ]
  m <- glm(spec$formula, data = dfit, family = quasipoisson)
  n <- nobs(m)
  nw_lag <- nw_lag %||% floor(4 * (n / 100)^(2 / 9))
  V <- if (hac) NeweyWest(m, lag = nw_lag, prewhite = FALSE, adjust = TRUE) else vcov(m)
  b <- coef(m)
  se <- sqrt(diag(V))
  z <- qnorm(1 - (1 - level) / 2)
  rr <- function(j) c(exp(b[j]), se[j], exp(b[j] - z * se[j]), exp(b[j] + z * se[j]),
                      2 * pnorm(-abs(b[j] / se[j])))
  r2 <- rr("X")
  r3 <- rr("Tpost")
  row <- data.frame(model = sprintf("quasi-Poisson%s", if (hac) sprintf(" + Newey–West(%d)", nw_lag) else ""),
                    b2 = r2[1], b2_se = r2[2], b2_lower = r2[3], b2_upper = r2[4], b2_p = r2[5],
                    b3 = r3[1], b3_se = r3[2], b3_lower = r3[3], b3_upper = r3[4], b3_p = r3[5],
                    scale = "率比（RR）；SE 為 log 尺度",
                    note = sprintf("離散參數 = %.2f", summary(m)$dispersion))
  list(model = m, row = row, dispersion = summary(m)$dispersion)
}

print.its_sens <- function(x, ...) {
  cat("── 敏感度分析：β2（水準變化）與 β3（斜率變化）──\n")
  for (k in seq_len(nrow(x))) {
    r <- x[k, ]
    if (is.na(r$b2)) {
      cat(sprintf("• %s：%s\n", r$model, r$note))
      next
    }
    cat(sprintf("• %s［%s］\n    β2 = %.3f（95%% CI %.3f ~ %.3f，p = %.3g）\n    β3 = %.4f（95%% CI %.4f ~ %.4f，p = %.3g）%s\n",
                r$model, r$scale, r$b2, r$b2_lower, r$b2_upper, r$b2_p,
                r$b3, r$b3_lower, r$b3_upper, r$b3_p,
                if (nzchar(r$note)) paste0("\n    ", r$note) else ""))
  }
  invisible(x)
}

# ── 10. 繁中解釋句 ──────────────────────────────────────────────
its_explain <- function(eff, prep, y_label = NULL, unit = "", alpha = 0.05,
                        period = NULL) {
  i <- prep$info
  y_label <- y_label %||% i$y_col
  period <- period %||% c(`1` = "年", `4` = "季", `7` = "天", `12` = "月", `52` = "週")[as.character(i$freq)]
  t <- eff$table
  f <- function(x) formatC(x, format = "f", digits = 2, big.mark = ",")
  ci <- function(r) sprintf("95%% CI %s 到 %s", f(r$lower), f(r$upper))
  sig <- function(r) r$p_value < alpha
  dirw <- function(x, up = "增加", down = "減少") if (x >= 0) up else down

  r1 <- t[1, ]; r2 <- t[2, ]; r3 <- t[3, ]; r4 <- t[4, ]
  s <- c(
    sprintf("**介入前趨勢**：每%s平均%s %s%s（%s）。", period, dirw(r1$estimate),
            f(abs(r1$estimate)), unit, ci(r1)),
    sprintf("**水準變化（β2）**：介入當下，%s比依介入前趨勢預期的值%s %s%s（%s，%s）。%s",
            y_label, dirw(r2$estimate), f(abs(r2$estimate)), unit, ci(r2),
            .fp(r2$p_value),
            if (sig(r2)) "有證據顯示介入當下水準改變。"
            else "信賴區間包含 0，沒有足夠證據顯示介入當下水準改變（不等於證明沒有改變）。"),
    sprintf("**斜率變化（β3）**：介入後每%s的趨勢比介入前%s %s%s（%s，%s）。%s",
            period, dirw(r3$estimate, "多", "少"), f(abs(r3$estimate)), unit, ci(r3),
            .fp(r3$p_value),
            if (sig(r3)) "有證據顯示趨勢改變。" else "沒有足夠證據顯示趨勢改變。"),
    sprintf("**介入後斜率**：每%s %s%s（%s）。", period, f(r4$estimate), unit, ci(r4))
  )
  rel <- eff$relative
  s <- c(s, sprintf("**整體相對變化**：介入後期間，%s平均比反事實%s %.1f%%；相對變化 = %+.1f%%（95%% CI %+.1f%% 到 %+.1f%%，%s）。%s",
                    y_label, dirw(rel$estimate, "高", "低"), 100 * abs(rel$estimate),
                    100 * rel$estimate, 100 * rel$lower, 100 * rel$upper,
                    .fp(rel$p_value),
                    "這是把整個介入後期間的效果加總後的比例，跟 β2、β3 的單點檢定不同，結論可能不一樣。"))
  s <- c(s, "ITS 是準實驗設計：以上解釋成立的前提是「沒有介入的話，介入前趨勢會延續」，而且介入點附近沒有其他同時發生的事件。")
  paste(s, collapse = "\n\n")
}

.fp <- function(p) ifelse(p < 1e-4, "p < 0.0001", paste0("p = ", formatC(signif(p, 3), format = "fg", digits = 3)))

its_explain_explore <- function(ex, prep) {
  i <- prep$info
  s <- sprintf("介入前有 %d 個觀測、介入後有 %d 個%s。", i$n_pre, i$n_post,
               if (i$n_excluded) sprintf("（另有 %d 期過渡期不納入配適）", i$n_excluded) else "")
  if (min(i$n_pre, i$n_post) < 12) s <- c(s, "其中一段少於 12 點，檢定力有限，結果要保守解讀。")
  if (!is.na(ex$p_season)) {
    s <- c(s, sprintf("介入前資料的季節性 F 檢定 %s，%s",
                      .fp(ex$p_season),
                      if (ex$p_season < .05) "季節性明顯，模型必須控制季節，否則前後兩段涵蓋的月份不同會扭曲效果。"
                      else "季節性不明顯，但頻率 > 1 時仍建議保留季節項。"))
  }
  s <- c(s, sprintf("介入前趨勢的二次項檢定 %s，%s", .fp(ex$p_pre_curvature),
                    if (ex$p_pre_curvature < .05) "介入前趨勢明顯彎曲，用直線延伸的反事實可能不可信，建議做「只用較近的介入前資料」敏感度分析。"
                    else if (ex$p_pre_curvature < .10) "在邊緣（.05–.10），請對照時序圖的 loess 線判斷，並建議做「只用較近的介入前資料」敏感度分析。"
                    else "沒有證據顯示介入前趨勢彎曲，線性趨勢的假設尚可接受。"))
  paste(s, collapse = "")
}

its_explain_diag <- function(dg) {
  tt <- dg$tests
  lines <- sprintf("%s：統計量 = %.3f，%s", tt$test, tt$statistic, .fp(tt$p_value))
  concl <- if (dg$autocorrelation) {
    sprintf("自相關檢定中至少一項 p < %.2f，殘差**有自相關**（lag 1 的 ACF = %.2f）。此時 OLS 的係數仍不偏，但標準誤會被低估，信賴區間太窄、p 值太小，所以改用把自相關放進模型的 GLS。",
            dg$alpha, dg$acf_lag1)
  } else {
    sprintf("自相關檢定的 p 值都 ≥ %.2f，沒有足夠證據顯示殘差有自相關（lag 1 的 ACF = %.2f）。這不代表證明獨立，但目前的模型足以處理序列相關。",
            dg$alpha, dg$acf_lag1)
  }
  paste(c(paste0("- ", lines), "", concl), collapse = "\n")
}

its_explain_sens <- function(sens) {
  s_all <- sens[!is.na(sens$b2) & sens$scale == "原尺度（差值）", ]
  if (!nrow(s_all)) return("")
  # 方向：所有原尺度模型都看；顯著性：排除已知標準誤偏小的對照模型
  # （主要模型為 GLS 時的 OLS、以及所有 Newey–West），見 VALIDATION.md 的模擬。
  main_is_gls <- any(grepl("^主要模型：.*GLS", s_all$model))
  contrast <- grepl("Newey", s_all$model) | (main_is_gls & grepl("^OLS", s_all$model))
  s <- s_all[!contrast, ]
  same_dir <- all(sign(s_all$b2) == sign(s_all$b2[1]))
  n_sig <- sum(s$b2_upper < 0 | s$b2_lower > 0)
  out <- sprintf("在 %d 個原尺度模型中，β2 的估計介於 %.2f 到 %.2f，%s。",
                 nrow(s_all), min(s_all$b2), max(s_all$b2),
                 if (same_dir) "方向全部一致" else "**方向不一致**，結論不穩健")
  if (any(contrast)) {
    out <- c(out, sprintf("%s的標準誤已知偏小（模擬中誤判率約 20–30%%），只當對照，不列入顯著性判斷。",
                          paste(s_all$model[contrast], collapse = "、")))
  }
  if (nrow(s) <= 1) {
    out <- c(out, "撇開只當對照的模型，只剩主要模型本身，沒有其他設定可以比較，所以**無法判斷結論是否穩健**；只能說對照模型的方向與主要模型一致。")
  } else {
  out <- c(out, sprintf("其餘 %d 個模型中，%d 個的 95%% CI 不含 0。", nrow(s), n_sig))
  out <- c(out, if (n_sig == 0) "效果方向一致但都不顯著：資料支持「可能有效果」，但不確定性大。"
                else if (n_sig == nrow(s)) "這些設定都顯著，結論穩健。"
                else "顯著與否取決於模型設定：效果方向可信，但「是否顯著」不穩健，報告時應以主要模型為準並同時呈現這些結果。")
  }
  qp <- sens[grepl("quasi-Poisson", sens$model) & !is.na(sens$b2), ]
  if (nrow(qp)) {
    out <- c(out, sprintf("計數模型（quasi-Poisson）的水準變化率比 = %.3f（95%% CI %.3f 到 %.3f），即介入當下約%s %.1f%%；它的標準誤用 Newey–West，CI 可能偏窄。",
                          qp$b2, qp$b2_lower, qp$b2_upper,
                          if (qp$b2 < 1) "減少" else "增加", 100 * abs(qp$b2 - 1)))
  }
  if (any(grepl("Newey", sens$note))) {
    out <- c(out, "Newey–West 的標準誤比 OLS 還小，這在介入後點數少時常見，不能據此認為結果更精確。")
  }
  paste(out, collapse = "")
}

# ── 11. 保存與讀取分析狀態 ─────────────────────────────────────
its_save_state <- function(state, path) {
  state$saved_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  state$r_version <- R.version.string
  saveRDS(state, path)
  invisible(path)
}
