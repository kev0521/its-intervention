---
name: its-intervention
description: >-
  介入分析（中斷時間序列 ITS／分段迴歸）的互動分析流程：判斷一個已知事件（政策、法規、新療程、建設）
  發生後，時間序列的水準或趨勢有沒有改變。用 R 做資料檢查與探索、分段迴歸、殘差自相關診斷、
  依規則選 OLS 或 GLS AR(p)、估計水準變化與斜率變化及反事實、敏感度分析，
  最後產出繁體中文、可重現的 Quarto HTML 報告。使用者提到介入分析、中斷時間序列、
  interrupted time series、ITS、segmented regression、分段迴歸、政策評估、
  「某事件前後有沒有改變」「政策實施後有沒有效」並提供時間序列資料時，使用本 skill。
---

# 介入分析：中斷時間序列（ITS）

本 skill 帶使用者從一條時間序列和一個已知介入時間點，走到效果估計與報告。
**每個確認點都要停下來等使用者明確回答。**「幫我分析」不代表同意所有步驟。

**所有分析數字都必須來自 `scripts/its.R` 的函數。** 不心算、不捏造、不改用 Python、不臨時補算函數沒提供的統計量。
自己寫的呼叫程式只負責讀檔、呼叫函數、保存結果。缺功能時先說明限制，不用未審查的程式繞過。

## 執行方式

- skill 路徑記為 `SKILL_DIR`（本檔所在資料夾）。在使用者的工作資料夾建立分析腳本（例如 `its_work/run.R`），開頭 `source(file.path(SKILL_DIR, "scripts/its.R"), encoding = "UTF-8"); its_configure()`。
- **用 `Rscript 檔案.R` 執行，不要用多行的 `Rscript -e`**（Windows 上可能崩潰）。
- 圖輸出到 `its_work/figures/`，畫完要實際開圖檢查再描述。
- 需要的套件：`nlme`、`sandwich`、`lmtest`、`ggplot2`；報告需要 Quarto（或 RStudio 內建的 Quarto）。
- 每階段的物件用 `saveRDS()` 存到 `its_work/`，下一階段的腳本再讀回來。
- 工作資料夾裡已經有 `its_work/` 時，**另建 `its_work_<YYYYMMDD-HHMM>/`**，不要跟舊分析的檔案混在一起；最後告訴使用者這次用的是哪個資料夾。
- **函數介面看下表就好，不需要去讀 `its.R` 原始碼。**

### 函數速查（`scripts/its.R`）

| 函數 | 參數 | 回傳／重點 |
|---|---|---|
| `its_prepare` | `(data, time, y, t0, freq, transition = 0, covariates = NULL)`；`t0` 用 `"1983-02"`、日期或數字；`freq` 只能是 1、4、7、12、52 | `$status`（`"ok"`／`"stop"`）、`$issues`、`$info`（`n_pre`、`n_post`…）、`$data` |
| `its_explore` | `(prep, out_dir, y_label)` | 時序圖、季節圖；`$suggestion$season`、`$suggestion$family`；說明文字用 `its_explain_explore(ex, prep)` |
| `its_spec` | `(prep, season = "dummy"/"fourier"/"none", K = 2, covariates = NULL)` | 模型設定 |
| `its_fit` | `(prep, spec, method = "ols"/"gls", p = 0, approved = TRUE)` | 沒有 `approved = TRUE` 會報錯 |
| `its_diagnose` | `(fit, prep, alpha = 0.05, out_dir, tag)` | 檢定表、圖 `diag_<tag>.png`；說明用 `its_explain_diag(dg)` |
| `its_select_ar` | `(prep, spec, p_max = 4)` | AIC／BIC 表、`$best_p` |
| `its_effects` | `(fit, prep, k = NULL, level = 0.95)` | `$table`（β1、β2、β3、β1+β3、各期效果）、`$relative`；說明用 `its_explain(eff, prep, y_label, unit)` |
| `its_plot` | `(fit, prep, file, y_label, title)` | 觀測、配適、反事實線 |
| `its_sensitivity` | `(main, prep, which, covariates, transition, pre_window, data)`；`which` 可選 `"newey_west"`、`"covariates"`、`"transition"`、`"pre_window"`、`"quasipoisson"`；`transition` = 排除介入後前幾期；`pre_window` = 介入前保留最近幾**期**（月資料保留 8 年就是 96）；這兩項都要搭配 `data`（原始資料框） | 每個模型的 β2、β3；說明用 `its_explain_sens(sens)` |

### 報告的 `state` 欄位（`its_write_report(state, out_dir, SKILL_DIR, analysis_finished = TRUE)`）

`title`、`research_question`、`data_file`、`time`、`y`、`t0`、`freq`、`transition`、`covariates`、`y_label`、`unit`、`context_events`、`season`、`K`、`model_covariates`、`p_max`、`method`、`p`、`p_reason`、`sens = list(which, covariates, transition, pre_window)`、`completed_stages`、`decisions = data.frame(item, answer)`。完整範例見 `demo_seatbelts/make_report.R`。

## 範圍

| 範圍內 | 範圍外（說明理由後停止或建議其他方法） |
|---|---|
| 一條間隔相等的序列（年、季、月、週、日） | 介入時間未知、想找斷點 → `strucchange` |
| 一個已知介入時間點（可設過渡期） | 多個介入點、有控制序列 → 第一版不支援 |
| 連續或計數的結果變數 | 多個單位的面板資料 → 混合模型或 DiD |
| | 想預測未來 → 預測模型 |

方法細節見 `references/methods.md`，決策流程見 `references/decision-tree.md`。

## 提問方式：一律用選擇題

所有要使用者回答的問題（第 1 階段的設定、每個停止點的確認）**都用選擇題問**，不要用開放式提問。

- **有 `AskUserQuestion` 工具時一定要用它**。限制：每次最多 4 題、每題 2–4 個選項，工具會自動附「其他」讓使用者自由輸入，所以**不要自己加「其他」選項**。
- **選項要從資料產生**：欄位名稱、偵測到的頻率、AIC 排名等都直接放進選項，不要給空泛的選項。
- **推薦的選項放第一個**，標籤結尾加「（建議）」，`description` 寫一句理由（例如「偵測到間隔 1 個月」）。
- 選項之間必須互斥；可以複選的題目（例如共變數）設 `multiSelect: true`。
- 敏感度分析**不問**，依第 5 階段的條件自動執行。
- 每次呼叫完就停下來等回答；拿到答案後，先用一句話複述使用者的選擇，再進下一步。
- **沒有 `AskUserQuestion` 的環境**（例如網頁版）：用文字列出同樣的選擇題，選項編成 A、B、C…，最後一個固定是「其他（請說明）」，請使用者回覆代號。
- 使用者選「其他」並輸入內容時，照他的輸入做；內容不清楚時，再用一題選擇題追問，不要自己猜。

---

## 第 1 階段：讀資料與確認設定（一定要停）

1. 讀檔，列出欄位、筆數、時間範圍，推測時間欄位、頻率、結果變數型態。
2. **分兩次 `AskUserQuestion` 問完下面 8 題**（第一次 Q1–Q4，等回答後第二次 Q5–Q8），然後停止等待：

   | # | header | 問題 | 選項怎麼產生 | 複選 |
   |---|---|---|---|---|
   | Q1 | 結果變數 | 要分析哪一欄的變化？ | 數值欄位最多 4 個，名稱最像研究問題的放第一個（建議） | 否 |
   | Q2 | 時間欄位 | 哪一欄是時間？ | 可解析成日期或年份的欄位；只有一個候選時另一個選項為「資料已依時間排序，沒有時間欄」 | 否 |
   | Q3 | 介入時間 | 介入從哪一個時間點開始（第一個受影響的觀測）？ | 使用者提過日期時：「該日期所在期」「下一期（第一個完整期）」，各附理由。資料裡有 0/1 介入指標欄（例如 `law`）時，可以把「指標第一次等於 1 的那一期（依 X 欄）」列為建議選項，description 請使用者確認跟實際生效日一致。**兩者都沒有時不可猜**，選項為「我會在『其他』輸入」「還不確定，先只看探索圖」 | 否 |
   | Q4 | 資料頻率 | 資料的頻率是？ | 偵測到的頻率放第一個（建議），其餘從 年／季／月／週／日 補到 4 個 | 否 |
   | Q5 | 變數型態 | 結果變數是哪一種？ | 「計數（人數、件數）」「連續數值」；依資料是否全為非負整數決定哪個建議；平均計數 < 10 時在 description 註明會建議 quasi-Poisson | 否 |
   | Q6 | 生效方式 | 介入是立即生效，還是有過渡期？ | 「立即生效（建議）」「有 1 期過渡期」「有 2–3 期過渡期」 | 否 |
   | Q7 | 其他事件 | 介入前後有沒有其他事件，或定義、測量方式改變？ | 「沒有」「有（請在『其他』說明）」「不知道」 | 否 |
   | Q8 | 共變數 | 要不要把同期的變數當共變數（只用在敏感度分析）？ | 其餘數值欄位最多 3 個 ＋「都不加（建議）」；**排除**介入指標欄（例如 0/1 的 `law`）和其他結果變數（例如同一事件的其他死傷數） | 是 |

   - Q3 選「先只看探索圖」：只做第 2 階段的資料檢查與時序圖，畫完再用選擇題問介入時間。**不要自己替使用者決定介入時間。**
   - Q7 回答「有」或「不知道」：記進 `decisions`，結論與限制一定要提。
3. 研究問題不清楚時（例如使用者只附了檔案），在 Q5–Q8 之後**另外呼叫一次**，用一題選擇題問研究目的（例如「水準有沒有跳動」「趨勢有沒有改變」「兩者都看（建議）」）。

## 第 2 階段：資料檢查與探索（一定要停）

1. `prep <- its_prepare(data, time, y, t0, freq, transition, covariates)`，`print(prep)`。
   - `status == "stop"`：說明原因（缺值、重複、間隔不等、介入點不在範圍內、點數 < 3），**停止**，請使用者處理。不自動插補、不自動刪資料。
   - 有 `warn`：照實告知（例如任一段 < 12 點、檢定力有限），可以繼續。
2. `ex <- its_explore(prep, "its_work/figures", y_label)`，打開時序圖與季節圖，描述趨勢、季節性、介入前後的變化外觀。
3. 介入前趨勢的二次項檢定 p < .10 時，提醒反事實可能不可信，並預告第 5 階段會自動做「前段視窗」敏感度分析。
4. 提出模型設定建議：季節項（`ex$suggestion$season`）、模型家族（`ex$suggestion$family`），寫出模型公式並解釋 β2、β3。
5. **用選擇題問「是否依此設定配適模型？」並停止等待。** 選項：「依建議設定配適（建議）」「改季節項設定」「改模型家族」「先停在探索，不配適」。選「改季節項」時再追問一題：「月份／季別虛擬變數」「Fourier 項（K = 2）」「不加季節項」。

## 第 3 階段：OLS 與殘差診斷（自相關時要停）

1. `spec <- its_spec(prep, season, K)`；`m_ols <- its_fit(prep, spec, "ols", approved = TRUE)`。
2. `dg <- its_diagnose(m_ols, prep, out_dir = ...)`，打開診斷圖，用 `its_explain_diag(dg)` 的內容解釋。
3. 判定：
   - **沒有自相關** → 主要模型 = OLS，直接進第 4 階段。
   - **有自相關** → 說明為什麼不能只看 OLS（SE 被低估）。執行 `ar <- its_select_ar(prep, spec)`，呈現 AIC 表，建議 AIC 最小的 p。**用選擇題問「主要模型要用哪一個？」並停止等待。** 選項：「GLS AR(最佳 p)（建議）」附 AIC；「GLS AR(次佳 p)」附 AIC；「仍用 OLS（不建議：SE 會被低估）」。
4. 使用者確認後 `main <- its_fit(prep, spec, "gls", p, approved = TRUE)`，再對 `main` 做 `its_diagnose()`。標準化殘差仍有自相關時，照實說明並建議 ARIMA（範圍外），不要假裝已經處理好。

## 第 4 階段：效果估計與解釋

1. `eff <- its_effects(main, prep)`；`its_plot(main, prep, file, y_label)`，打開圖檢查。
2. 依 `its_explain(eff, prep, y_label, unit)` 的內容解釋，再用研究問題的語言補充：
   - 一定要報 β2、β3、介入後斜率、指定時間點的效果、整體相對變化，**全部附 95% CI**。
   - p ≥ .05 寫「沒有足夠證據顯示…」，不寫「沒有效果」；p < .05 也不寫「證明介入造成…」。
   - 整體相對變化跟 β2 的結論不同時，要解釋原因（一個是單點、一個是整段加總）。
   - 主要模型是 GLS 且 p 在 .03–.10 之間時，提醒 GLS 在小樣本略為寬鬆（`VALIDATION.md`：涵蓋率約 92–93%），結論要保守。
3. 解釋完**直接進第 5 階段，不問使用者要不要做敏感度分析**。

## 第 5 階段：敏感度分析（自動執行，不問使用者）

只要下面任一條件成立，就**一律執行所有符合條件的項目**，不提供「不做」或「自己挑」的選項。一項都不符合時，說明「沒有需要檢查的疑慮」，跳過這個階段。

| 條件 | 項目 | 預設參數 |
|---|---|---|
| OLS 殘差有自相關 | `newey_west`（只當對照） | 自動落後期 |
| 使用者在 Q8 選了共變數 | `covariates` | 使用者選的欄位 |
| 使用者在 Q6 選了過渡期，或介入後頭幾期在圖上明顯異常 | `transition` | 使用者選的期數；沒選時用 2 期 |
| 介入前趨勢的二次項 p < .10 | `pre_window` | 最近 8 個季節週期（月資料 96 期）；年資料或介入前不足 16 個週期時，用介入前的後一半 |
| 結果是計數資料 | `quasipoisson`（報率比） | — |

1. 先用一句話告訴使用者要做哪幾項、為什麼（對應上表的條件），然後執行 `sens <- its_sensitivity(main, prep, which, covariates, transition, pre_window, data)`。
2. 用 `its_explain_sens(sens)` 的內容說明：方向是否一致、有幾個模型顯著、結論穩不穩健。
3. Newey–West 的 SE 比 OLS 還小時，一定要提醒它在小樣本不可靠。
4. **主要模型是 GLS 時，「OLS（未處理自相關）」和所有 Newey–West 結果只當對照**，不能拿來說「結論穩健」（模擬中誤判率約 20–30%）。`its_explain_sens()` 已經自動排除，自己補充說明時也要遵守。
5. 使用者看完結果後，主動要求加做其他項目或改參數時，照做並更新報告；skill 不主動詢問。

## 結束與報告

使用者說結束，或流程走完時，先說「此次分析結束」，然後直接產生報告，不再多問一次。`decisions` 表記錄每題的**問題、選到的選項**（選「其他」時記使用者輸入的原文）。

1. 組 `state`（欄位見 `demo_seatbelts/make_report.R`）：設定、使用者的回答（`decisions` 表，用使用者原話）、`completed_stages`（只放**實際完成**的階段：`explore`、`ols`、`ar_select`、`final`、`sensitivity`）。
2. `source(file.path(SKILL_DIR, "scripts/finalize_report.R"))`；`its_write_report(state, out_dir, SKILL_DIR, analysis_finished = TRUE)`。
3. 報告是 HTML（單一檔案、圖嵌入），每段「程式碼 → 結果 → 繁中解釋」。未完成的階段整段移除。
4. 編譯失敗時保留 qmd 並照實回報錯誤，不宣稱完成。

## 重現性

- `its_configure()` 固定選項與種子；資料 MD5、R 與套件版本寫進報告。
- 資料、介入時間或模型設定改變後，後面的確認全部失效，回到第 2 階段。

## 寫作規則

- 繁體中文、台灣用語；保留必要的英文名詞與 R 函數名稱。
- 數字直接引用 R 輸出；表格用 Markdown。
- 預設 α = .05、雙尾、95% CI；使用者指定時照指定。
- ITS 是準實驗設計：結論都要搭配「趨勢延續」和「沒有同時事件」兩個假設一起說。

## 示範

`assets/seatbelts.csv`（英國安全帶法，1983-02 介入，月資料）是完整示範；`demo_seatbelts/` 有分析腳本和報告。
