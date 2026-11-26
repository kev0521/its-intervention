#!/usr/bin/env sh
# 由 README.md 產生單一檔案的 README.html（在 repo 根目錄執行：sh tools/build_readme_html.sh）
# 需要 Quarto；mermaid 流程圖會轉成 Quarto 的 {mermaid} 區塊一起嵌入。
set -e
cd "$(dirname "$0")/.."
tmp=README_build.qmd
{
  printf -- '---\ntitle: "介入分析（ITS）Skill 使用說明"\nlang: zh-TW\nformat:\n  html:\n    embed-resources: true\n    toc: true\n    theme: cosmo\n---\n\n'
  sed -e '1d' -e 's/^```mermaid$/```{mermaid}/' -e 's/^> \[!NOTE\]$/> **注意**/' README.md
} > "$tmp"
quarto render "$tmp" --output README.html
rm -f "$tmp"
echo "已產生 README.html"
