#!/usr/bin/env bash
# build_pdf.sh -- regenerate PDF + DOCX from markdown reports.
#
# STIX Two Text + Math (macOS native) handles Greek and most math
# symbols. A small set of rare chars (right-arrow, sqrt, checkmark,
# angle brackets) are substituted with LaTeX math-mode equivalents in
# a temp pdfsafe copy of the markdown before xelatex runs.
#
# Usage:
#   bash docs/build_pdf.sh                  # build both calibration reports
#   bash docs/build_pdf.sh HDD              # build only HDD report
#   bash docs/build_pdf.sh PDD              # build only PDD report
#   bash docs/build_pdf.sh KYLE             # build Kyle's onboarding doc
#
# Requirements:
#   - MacTeX (xelatex) at /Library/TeX/texbin
#   - pandoc (any recent version)
#   - macOS with STIX Two Text + STIX Two Math fonts (system default)
#   - python3

set -euo pipefail

export PATH="/Library/TeX/texbin:$PATH"
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$SCRIPT_DIR"

substitute_unicode() {
  local src="$1"
  local dst="$2"
  python3 - "$src" "$dst" <<'PY'
import sys
src, dst = sys.argv[1], sys.argv[2]
with open(src) as f:
    text = f.read()
replacements = {
    '→': r'$\to$',           # right arrow
    '√': r'$\sqrt{\,}$',     # sqrt
    '✓': r'$\checkmark$',    # checkmark
    '⟨': r'$\langle$',       # left angle bracket
    '⟩': r'$\rangle$',       # right angle bracket
}
for old, new in replacements.items():
    text = text.replace(old, new)
with open(dst, 'w') as f:
    f.write(text)
PY
}

build_one() {
  local base="$1"
  local title="$2"
  local date="$3"

  echo "=== $base ==="

  # DOCX -- no preprocessing needed (Word renders Unicode fine)
  pandoc "${base}.md" \
    -o "${base}.docx" \
    --toc --toc-depth=2 \
    --metadata title="${title}" \
    --metadata author="T. Mehlhorn" \
    --metadata date="${date}"
  echo "  -> ${base}.docx"

  # PDF -- preprocess for xelatex, then build
  substitute_unicode "${base}.md" "${base}_pdfsafe.md"
  pandoc "${base}_pdfsafe.md" \
    -o "${base}.pdf" \
    --pdf-engine=xelatex \
    --toc --toc-depth=2 \
    --metadata title="${title}" \
    --metadata author="T. Mehlhorn" \
    --metadata date="${date}" \
    -V geometry:margin=1in \
    -V mainfont="STIX Two Text" \
    -V mathfont="STIX Two Math" \
    -V monofont="Menlo" \
    2>&1 | (grep -v "Missing character" || true)
  rm -f "${base}_pdfsafe.md"
  echo "  -> ${base}.pdf"
}

TARGET="${1:-BOTH}"

if [ "$TARGET" = "BOTH" ] || [ "$TARGET" = "HDD" ]; then
  build_one Xcimer_HDD_calibration_report \
    "Helios HDD Calibration: WT_cthomas Production Anchor and RHINO Cross-Validation" \
    "2026-06-02"
fi

if [ "$TARGET" = "BOTH" ] || [ "$TARGET" = "PDD" ]; then
  build_one Xcimer_PDD_calibration_report \
    "Helios PDD Calibration: Foam-Burn Deficit and Root-Cause Analysis" \
    "2026-05-28"
fi

if [ "$TARGET" = "KYLE" ]; then
  build_one Kyle_summer_onboarding \
    "Coordinating with the helios_postprocess Project -- Summer 2026" \
    "2026-06-04"
fi

echo "Done."
