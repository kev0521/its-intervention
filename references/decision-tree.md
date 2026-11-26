# 決策樹

```mermaid
flowchart TD
  A(["使用者附上時間序列資料"]) --> Q1{{"第 1 階段 選擇題 ①：結果變數、時間欄、介入時間、頻率"}}
  Q1 --> Q2{{"選擇題 ②：變數型態、生效方式、其他事件、共變數"}}
  Q2 --> C["第 2 階段：資料檢查 its_prepare()"]
  C -->|"缺值、重複、間隔不等、點數 < 3"| X1["說明原因並停止"]
  C -->|"通過"| D["探索：時序圖、季節圖、趨勢彎曲檢定"]
  D --> Q3{{"選擇題 ③：依建議設定配適模型？"}}
  Q3 -->|"先停在探索"| X2["保留探索結果，結束"]
  Q3 -->|"配適"| F["第 3 階段：OLS 分段迴歸 ＋ 殘差診斷"]
  F --> G{"有自相關？"}
  G -->|"否（BG、Ljung–Box p ≥ .05）"| H["主要模型 = OLS"]
  G -->|"是"| I["AIC 比較 GLS AR(0)–AR(4)"]
  I --> Q4{{"選擇題 ④：主要模型用哪一個？"}}
  Q4 --> K["GLS AR(p) ＋ 標準化殘差檢查"]
  H --> L["第 4 階段：β2、β3、相對變化、反事實圖、繁中解釋"]
  K --> L
  L --> S{"有需要檢查的疑慮？"}
  S -->|"有（自相關、趨勢彎曲、計數、共變數、過渡期）"| N["第 5 階段：自動做所有符合條件的敏感度分析"]
  S -->|"都沒有"| O(["產生 HTML 報告"])
  N --> O
  classDef ask fill:#fff4e5,stroke:#d97706,color:#7c2d12
  classDef stop fill:#fde8e8,stroke:#d1242f,color:#7f1d1d
  classDef calc fill:#e8f1fe,stroke:#1f6feb,color:#0b3d91
  classDef io fill:#e9f7ef,stroke:#1a7f37,color:#0f5323
  class Q1,Q2,Q3,Q4 ask
  class X1,X2 stop
  class C,D,F,I,K,H,L,N,G,S calc
  class A,O io
```

## 判斷規則（文字版）

| 節點 | 規則 |
|---|---|
| 資料檢查 stop | 時間缺值或重複、間隔不等、Y 有缺值、介入點不在範圍內、任一段 < 3 點 |
| 資料檢查 warn | 任一段 < 12 點；介入前不到一個季節週期；平均計數 < 10 |
| 季節項建議 | 頻率 1 → none；4、12 → dummy；7、52 → fourier（K = 2） |
| 模型家族建議 | 計數且平均 < 10 → quasi-Poisson 為主；其他 → gaussian |
| 有自相關 | OLS 殘差的 Breusch–Godfrey 或 Ljung–Box 任一 p < α（DW 只當參考） |
| AR 階數 | ML 估計比較 AIC，取最小；差距 < 2 時可以選較簡單的模型，但要說明 |
| 趨勢彎曲 | 介入前二次項 p < .10 → 自動做「前段視窗」敏感度分析 |
| Newey–West | 只當敏感度分析；SE 比 OLS 小時加警告 |

## 範圍外：說明理由後停止

- 介入時間未知 → 斷點偵測（`strucchange`）
- 多個介入點、有控制序列、多個單位的面板資料 → 第一版不支援
- 想預測未來 → 預測模型
