# validate.R — its-intervention 驗收
# 執行位置：its-intervention/   指令：Rscript tests/validate.R
source("scripts/its.R", encoding = "UTF-8")
its_configure()
out_dir <- "tests/results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log <- data.frame(group = character(), case = character(), expected = character(),
                  observed = character(), pass = logical())
rec <- function(group, case, expected, observed, pass) {
  log <<- rbind(log, data.frame(group = group, case = case, expected = expected,
                                observed = observed, pass = pass))
  cat(sprintf("[%s] %s｜%s：%s\n", if (pass) "PASS" else "FAIL", group, case, observed))
}

# ── A. 模擬：已知真值，比較 OLS / OLS+NW / GLS AR(1) ───────────
sim_one <- function(b2, b3, phi = 0.6, n_pre = 60, n_post = 36, sigma = 5) {
  n <- n_pre + n_post
  T <- seq_len(n)
  X <- as.integer(T > n_pre)
  Tpost <- pmax(T - (n_pre + 1), 0) * X
  month <- ((T - 1) %% 12) + 1
  e <- as.numeric(arima.sim(list(ar = phi), n, sd = sigma * sqrt(1 - phi^2)))
  y <- 100 + 0.2 * T + b2 * X + b3 * Tpost + 8 * sin(2 * pi * month / 12) + e
  start <- as.Date("2015-01-01")
  data.frame(month = format(seq(start, by = "month", length.out = n), "%Y-%m"), y = y)
}
t0_sim <- format(seq(as.Date("2015-01-01"), by = "month", length.out = 61)[61], "%Y-%m")

run_sim <- function(R, b2, b3, label) {
  res <- vector("list", R)
  for (r in seq_len(R)) {
    d <- sim_one(b2, b3)
    p <- its_prepare(d, "month", "y", t0_sim, 12)
    sp <- its_spec(p, "dummy")
    m_o <- its_fit(p, sp, "ols", approved = TRUE)
    m_g <- tryCatch(its_fit(p, sp, "gls", 1, approved = TRUE), error = function(e) NULL)
    get <- function(e, nm) {
      t <- e$table
      k <- c(b2 = 2, b3 = 3)[nm]
      c(est = t$estimate[k], lo = t$lower[k], hi = t$upper[k])
    }
    effs <- list(OLS = its_effects(m_o, p),
                 `OLS + Newey–West` = its_effects(m_o, p, vcov_type = "newey_west"))
    if (!is.null(m_g)) effs$`GLS AR(1)` <- its_effects(m_g, p)
    res[[r]] <- do.call(rbind, lapply(names(effs), function(nm) {
      rbind(data.frame(method = nm, par = "β2", t(get(effs[[nm]], "b2")), truth = b2),
            data.frame(method = nm, par = "β3", t(get(effs[[nm]], "b3")), truth = b3))
    }))
  }
  z <- do.call(rbind, res)
  names(z)[3:5] <- c("est", "lo", "hi")
  z$cover <- z$lo <= z$truth & z$hi >= z$truth
  z$reject <- z$lo > 0 | z$hi < 0
  agg <- aggregate(cbind(est, cover, reject) ~ method + par, data = z, FUN = mean)
  agg$truth <- ifelse(agg$par == "β2", b2, b3)
  agg$bias <- agg$est - agg$truth
  agg$scenario <- label
  agg$R <- R
  agg
}

set.seed(20261118)
R <- 1000
cat("模擬中（R =", R, "，每個情境）……\n")
sim_eff <- run_sim(R, b2 = -10, b3 = -0.5, "有效果（β2 = −10，β3 = −0.5）")
sim_null <- run_sim(R, b2 = 0, b3 = 0, "沒有效果（β2 = β3 = 0）")
sim <- rbind(sim_eff, sim_null)
write.csv(sim, file.path(out_dir, "simulation.csv"), row.names = FALSE)
print(sim[, c("scenario", "method", "par", "truth", "est", "bias", "cover", "reject")], digits = 3)

g <- sim_eff[sim_eff$method == "GLS AR(1)", ]
rec("模擬", "GLS 估計偏誤", "β2、β3 平均估計接近真值（|偏誤| < 10% 真值）",
    sprintf("β2 偏誤 %.3f、β3 偏誤 %.4f", g$bias[g$par == "β2"], g$bias[g$par == "β3"]),
    all(abs(g$bias) < 0.1 * abs(g$truth)))
rec("模擬", "GLS 95% CI 涵蓋率", "介於 90%–98%",
    paste(sprintf("%s %.1f%%", g$par, 100 * g$cover), collapse = "、"),
    all(g$cover > 0.90 & g$cover < 0.98))
o <- sim_eff[sim_eff$method == "OLS", ]
rec("模擬", "OLS 涵蓋率偏低（示範問題）", "明顯低於 95%",
    paste(sprintf("%s %.1f%%", o$par, 100 * o$cover), collapse = "、"),
    all(o$cover < 0.90))
gn <- sim_null[sim_null$method == "GLS AR(1)", ]
on <- sim_null[sim_null$method == "OLS", ]
rec("模擬", "虛無情境的型一錯誤", "GLS 接近 5%，OLS 明顯偏高",
    sprintf("GLS %s；OLS %s",
            paste(sprintf("%s %.1f%%", gn$par, 100 * gn$reject), collapse = "、"),
            paste(sprintf("%s %.1f%%", on$par, 100 * on$reject), collapse = "、")),
    all(gn$reject < 0.10) && all(on$reject > 0.10))

# AR 階數選擇
set.seed(1)
picks <- replicate(60, {
  d <- sim_one(-10, -0.5)
  p <- its_prepare(d, "month", "y", t0_sim, 12)
  its_select_ar(p, its_spec(p, "dummy"), p_max = 3)$best_p
})
rec("模擬", "AIC 選 AR 階數（真值 p = 1）", "最常選到 p = 1",
    paste(sprintf("p=%s：%d 次", names(table(picks)), table(picks)), collapse = "、"),
    as.integer(names(which.max(table(picks)))) == 1)

# 涵蓋率圖
library(ggplot2)
sim$lab <- paste(sim$par)
gp <- ggplot(sim, aes(100 * cover, method, colour = par)) +
  geom_vline(xintercept = 95, linetype = "dashed") +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  facet_wrap(~ scenario) +
  labs(x = "95% CI 的實際涵蓋率（%）", y = NULL, colour = NULL,
       title = sprintf("模擬驗證：AR(1) 誤差（φ = 0.6），每情境 %d 次", R)) +
  theme_minimal() + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "coverage.png"), gp, width = 9, height = 3.6, dpi = 120)

# ── B. 錯誤案例：應該停止 ─────────────────────────────────────
base <- read.csv("assets/seatbelts.csv")
expect_stop <- function(case, data, t0, freq = 12, pattern) {
  p <- its_prepare(data, "month", "DriversKilled", t0, freq)
  msg <- paste(p$issues$message[p$issues$level == "stop"], collapse = "；")
  rec("錯誤案例", case, paste("停止：", pattern), if (nzchar(msg)) msg else "（沒有停止）",
      p$status == "stop" && grepl(pattern, msg))
}
expect_stop("介入點早於序列", base, "1960-01", pattern = "範圍")
expect_stop("介入點晚於序列", base, "1990-01", pattern = "範圍")
d <- base; d$month[10] <- d$month[9]
expect_stop("時間重複", d, "1983-02", pattern = "重複")
expect_stop("間隔不等（少一個月）", base[-50, ], "1983-02", pattern = "間隔")
d <- base; d$DriversKilled[30] <- NA
expect_stop("結果變數有缺值", d, "1983-02", pattern = "缺值")
expect_stop("介入後只有 2 點", base, "1984-11", pattern = "少於 3 點")
p <- its_prepare(base, "month", "DriversKilled", "1983-02", 12)
e <- tryCatch({its_fit(p, its_spec(p, "dummy"), "ols"); "（沒有報錯）"},
              error = function(e) conditionMessage(e))
rec("錯誤案例", "未確認就配適（approved = FALSE）", "報錯拒絕配適", e, grepl("尚未確認", e))
e <- tryCatch({its_spec(its_prepare(read.csv("assets/nile.csv"), "year", "flow", 1899, 1), "dummy"); "（沒有報錯）"},
              error = function(e) conditionMessage(e))
rec("錯誤案例", "年資料要求季節項", "報錯", e, grepl("年資料", e))

# ── C. 邊界案例：警告但可以繼續 ───────────────────────────────
p <- its_prepare(base, "month", "DriversKilled", "1984-05", 12)
rec("邊界案例", "介入後只有 8 點", "警告 < 12 點，狀態 ok",
    paste(p$status, paste(p$issues$message, collapse = "；")),
    p$status == "ok" && any(grepl("少於 12", p$issues$message)))
p <- its_prepare(base[sample(nrow(base)), ], "month", "DriversKilled", "1983-02", 12)
rec("邊界案例", "資料順序打亂", "自動排序並說明",
    paste(p$status, paste(p$issues$message, collapse = "；")),
    p$status == "ok" && any(grepl("重新排序", p$issues$message)) && p$info$T0 == 170)
p <- its_prepare(base, "month", "DriversKilled", "1983-01-15", 12)
rec("邊界案例", "介入時間不在月初", "以下一個觀測當起點並說明",
    paste(format(p$info$t0_time), paste(p$issues$message, collapse = "；")),
    p$status == "ok" && p$info$T0 == 170)
nile <- read.csv("assets/nile.csv")
p <- its_prepare(nile, "year", "flow", 1899, 1)
m <- its_fit(p, its_spec(p, "none"), "ols", approved = TRUE)
ef <- its_effects(m, p)
rec("邊界案例", "年資料（Nile，無季節）", "流程可完成",
    sprintf("β2 = %.1f（%.1f ~ %.1f）", ef$table$estimate[2], ef$table$lower[2], ef$table$upper[2]),
    is.finite(ef$table$estimate[2]))
wk <- data.frame(week = format(seq(as.Date("2020-01-06"), by = "week", length.out = 156)),
                 y = rpois(156, 50) + rep(c(0, 10), c(100, 56)))
p <- its_prepare(wk, "week", "y", format(as.Date("2020-01-06") + 7 * 100), 52)
m <- its_fit(p, its_spec(p, "fourier", K = 2), "ols", approved = TRUE)
rec("邊界案例", "週資料 + Fourier 季節項", "流程可完成，β2 接近 10",
    sprintf("β2 = %.2f", coef(m)["X"]), abs(coef(m)["X"] - 10) < 5)
p <- its_prepare(base, "month", "DriversKilled", "1983-02", 12, transition = 3)
rec("邊界案例", "過渡期 3 期", "排除 3 點、介入後剩 20 點",
    sprintf("排除 %d 點、介入後 %d 點", p$info$n_excluded, p$info$n_post),
    p$info$n_excluded == 3 && p$info$n_post == 20)
lc <- data.frame(month = base$month, y = rpois(192, 4))
p <- its_prepare(lc, "month", "y", "1983-02", 12)
rec("邊界案例", "低計數（平均約 4）", "警告並建議 quasi-Poisson",
    paste(p$issues$message, collapse = "；"), any(grepl("quasi-Poisson", p$issues$message)))

# ── D. 一致性：skill 結果與獨立手算一致 ─────────────────────
p <- its_prepare(base, "month", "DriversKilled", "1983-02", 12)
m <- its_fit(p, its_spec(p, "dummy"), "gls", 1, approved = TRUE)
dd <- data.frame(y = base$DriversKilled, T = 1:192, X = as.integer(1:192 >= 170),
                 Tpost = pmax(1:192 - 170, 0), month = factor(rep(1:12, 16)))
ref <- gls(y ~ T + X + Tpost + month, data = dd, correlation = corAR1(form = ~ T))
rec("一致性", "GLS AR(1) 係數與手寫 nlme 模型一致", "β2、β3 差 < 1e-6",
    sprintf("差 = %.2e、%.2e", abs(coef(m)["X"] - coef(ref)["X"]), abs(coef(m)["Tpost"] - coef(ref)["Tpost"])),
    abs(coef(m)["X"] - coef(ref)["X"]) < 1e-6 && abs(coef(m)["Tpost"] - coef(ref)["Tpost"]) < 1e-6)

write.csv(log, file.path(out_dir, "validation_log.csv"), row.names = FALSE)
cat(sprintf("\n共 %d 項，通過 %d 項\n", nrow(log), sum(log$pass)))
