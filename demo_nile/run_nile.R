# Nile 示範：年資料、無季節性；1898 年起 Aswan 水壩工程後流量下降
# 執行位置：its-intervention/
source("scripts/its.R", encoding = "UTF-8"); its_configure()
source("scripts/finalize_report.R", encoding = "UTF-8")
d <- read.csv("assets/nile.csv")
prep <- its_prepare(d, "year", "flow", 1899, 1)
print(prep)
ex <- its_explore(prep, "demo_nile/out", "年流量（10^8 m³）"); print(ex)
sp <- its_spec(prep, "none")
m_ols <- its_fit(prep, sp, "ols", approved = TRUE)
dg <- its_diagnose(m_ols, prep, out_dir = "demo_nile/out"); print(dg)
ar <- its_select_ar(prep, sp); print(ar)
method <- if (dg$autocorrelation) "gls" else "ols"
p <- if (dg$autocorrelation) ar$best_p else 0
main <- its_fit(prep, sp, method, p, approved = TRUE)
if (method == "gls") print(its_diagnose(main, prep, out_dir = "demo_nile/out"))
e <- its_effects(main, prep); print(e)
its_plot(main, prep, "demo_nile/out/its_main.png", "年流量（10^8 m³）")
cat(its_explain(e, prep, "年流量", " 億立方公尺"), "\n")

state <- list(
  title = "尼羅河年流量在 1898 年後的變化",
  research_question = "1898 年後（Aswan 水壩工程與灌溉擴張），尼羅河 Aswan 的年流量水準是否下降？",
  data_file = "assets/nile.csv", time = "year", y = "flow", t0 = 1899, freq = 1,
  transition = 0, covariates = NULL, y_label = "年流量（10^8 m³）", unit = " 億立方公尺",
  context_events = "文獻認為 1898 年附近的流量下降與水壩工程及氣候變化都有關，無法完全歸因於單一事件",
  season = "none", K = 2, model_covariates = NULL, p_max = 4,
  method = method, p = p, p_reason = "",
  sens = list(which = c("newey_west"), covariates = NULL, transition = NULL, pre_window = NULL),
  completed_stages = c("explore", "ols", if (method == "gls") "ar_select", "final", "sensitivity"),
  decisions = data.frame(item = c("結果變數", "介入起點", "資料頻率", "模型設定"),
                         answer = c("flow", "1899（1898 年工程開始後第一年）", "年（1）",
                                    sprintf("分段迴歸、無季節項、%s", .its_label(main))))
)
html <- its_write_report(state, "demo_nile/report", ".", analysis_finished = TRUE)
cat("報告：", html, "\n")
