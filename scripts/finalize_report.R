# finalize_report.R — 分析結束後產生 Quarto 報告（HTML）
# 只放實際完成的階段；未完成的階段整段移除，不留空章節。
#
# 用法：
#   source("<skill>/scripts/finalize_report.R")
#   its_write_report(state, out_dir = "its_report", skill_dir = "<skill 路徑>")

.find_quarto <- function() {
  q <- Sys.which("quarto")
  if (nzchar(q)) return(unname(q))
  cand <- c("C:/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe",
            "C:/Program Files/Quarto/bin/quarto.exe",
            "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto",
            "/usr/local/bin/quarto")
  hit <- cand[file.exists(cand)]
  if (!length(hit)) stop("找不到 Quarto；請安裝 Quarto 或 RStudio")
  hit[1]
}

its_write_report <- function(state, out_dir, skill_dir, analysis_finished = FALSE,
                             render = TRUE, name = "its_report") {
  if (!isTRUE(analysis_finished)) stop("分析尚未結束，不能產生報告")
  valid <- c("explore", "ols", "ar_select", "final", "sensitivity")
  st <- state$completed_stages
  if (!length(st) || any(!st %in% valid)) stop("completed_stages 必須是：", paste(valid, collapse = ", "))
  needs <- list(ols = "explore", ar_select = "ols", final = "ols", sensitivity = "final")
  for (s in st) if (!all(needs[[s]] %in% st)) stop("階段 ", s, " 缺少前置階段 ", needs[[s]])
  if ("final" %in% st && state$method == "gls" && !"ar_select" %in% st) {
    stop("主要模型是 GLS，但沒有完成 ar_select 階段")
  }

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  file.copy(file.path(skill_dir, "scripts", "its.R"), file.path(out_dir, "its.R"), overwrite = TRUE)
  data_src <- normalizePath(state$data_file, mustWork = TRUE)
  file.copy(data_src, file.path(out_dir, basename(data_src)), overwrite = TRUE)
  state$data_file <- basename(data_src)
  source(file.path(skill_dir, "scripts", "its.R"), encoding = "UTF-8", local = TRUE)
  its_save_state(state, file.path(out_dir, "state.rds"))

  lines <- readLines(file.path(skill_dir, "assets", "report.qmd"), encoding = "UTF-8", warn = FALSE)
  lines <- gsub("{{TITLE}}", state$title, lines, fixed = TRUE)
  lines <- gsub("{{DATE}}", format(Sys.Date()), lines, fixed = TRUE)
  keep <- rep(TRUE, length(lines))
  cur <- NULL
  for (k in seq_along(lines)) {
    if (grepl("^<!-- stage: [a-z_]+ -->$", lines[k])) {
      cur <- sub("^<!-- stage: ([a-z_]+) -->$", "\\1", lines[k])
      keep[k] <- FALSE
    } else if (lines[k] == "<!-- /stage -->") {
      cur <- NULL
      keep[k] <- FALSE
    } else if (!is.null(cur) && !cur %in% st) {
      keep[k] <- FALSE
    }
  }
  qmd <- file.path(out_dir, paste0(name, ".qmd"))
  con <- file(qmd, open = "w", encoding = "UTF-8")
  writeLines(lines[keep], con)
  close(con)
  if (!render) return(invisible(qmd))

  old <- setwd(out_dir)
  on.exit(setwd(old))
  status <- system2(.find_quarto(), c("render", shQuote(basename(qmd))),
                    stdout = TRUE, stderr = TRUE)
  html <- file.path(out_dir, paste0(name, ".html"))
  if (!file.exists(file.path(basename(html)))) {
    stop("Quarto 編譯失敗：\n", paste(tail(status, 30), collapse = "\n"))
  }
  invisible(normalizePath(file.path(basename(html))))
}
