"""重新打包 its-intervention.zip（在任何位置執行：python tools/build_skill.py）。

只包含 skill 執行需要的檔案；README、期中報告、tools 不打包。
zip 內的根資料夾為 its-intervention/，解壓到 ~/.claude/skills/ 即可使用。
"""
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parent.parent
INCLUDE = ["SKILL.md", "VALIDATION.md", "scripts", "assets", "references",
           "tests", "demo_nile", "demo_seatbelts"]
SKIP_DIRS = {"out", "report", "__pycache__"}

out = ROOT / "its-intervention.zip"
with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as zf:
    for item in INCLUDE:
        p = ROOT / item
        files = [p] if p.is_file() else sorted(x for x in p.rglob("*") if x.is_file())
        for f in files:
            rel = f.relative_to(ROOT)
            if SKIP_DIRS & set(rel.parts[:-1]):
                continue
            zf.write(f, Path("its-intervention") / rel)
print(f"已產生 {out.name}")
