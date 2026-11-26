# Seatbelts 示範：用分析狀態產生最終報告
# 執行位置：its-intervention/
source("scripts/finalize_report.R", encoding = "UTF-8")

state <- list(
  title = "英國安全帶法對駕駛死亡人數的影響",
  research_question = "1983 年 1 月 31 日英國強制前座繫安全帶的法律生效後，每月駕駛死亡人數的水準或趨勢是否改變？",
  data_file = "assets/seatbelts.csv",
  time = "month", y = "DriversKilled", t0 = "1983-02", freq = 12,
  transition = 0, covariates = c("kms", "PetrolPrice"),
  y_label = "每月駕駛死亡數", unit = " 人",
  context_events = "1973–74 年石油危機使介入前趨勢先升後降；同期交通量（kms）與油價（PetrolPrice）也在變動，以敏感度分析控制",
  season = "dummy", K = 2, model_covariates = NULL,
  p_max = 4, method = "gls", p = 1, p_reason = "",
  sens = list(which = c("newey_west", "covariates", "transition", "pre_window", "quasipoisson"),
              covariates = c("kms", "PetrolPrice"), transition = 2, pre_window = 96),
  completed_stages = c("explore", "ols", "ar_select", "final", "sensitivity"),
  decisions = data.frame(
    item = c("結果變數", "介入起點", "資料頻率", "結果型態", "過渡期", "其他事件",
             "模型設定", "自相關結構", "敏感度分析"),
    answer = c("DriversKilled", "1983-02（法律 1983-01-31 生效，第一個完整月）", "月（12）",
               "計數，但平均約 120，用線性模型＋quasi-Poisson 敏感度分析", "無（立即生效）",
               "石油危機；交通量與油價", "分段迴歸 + 月份虛擬變數", "GLS AR(1)（AIC 最小）",
               "Newey–West、共變數、過渡期 2 個月、只用 1975 年後資料、quasi-Poisson")
  )
)

html <- its_write_report(state, out_dir = "demo_seatbelts/report", skill_dir = ".",
                         analysis_finished = TRUE)
cat("報告：", html, "\n")
