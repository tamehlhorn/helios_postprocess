#!/usr/bin/env bash
# ============================================================================
# run_pdd_scan.sh
# ---------------
# Drive a PDD geometric-parameter scan on the Studio:
#   1. For each .rhw under $SCAN_DIR matching a pattern, run Helios in batch
#      mode (skip if .exo already exists)
#   2. After all runs complete, invoke dump_burn_region_density.py per run
#      to extract no-burn implosion metrics into a single CSV
#
# Designed to complement examples/make_pdd_scan_rhw.py for the foam-burn
# recovery design study (FL_prism=0.012 with geometric defocus tuning).
#
# Usage
# -----
#   bash examples/run_pdd_scan.sh
#       Run the default pattern (wf_fl012_*.rhw under $SCAN_DIR).
#
#   bash examples/run_pdd_scan.sh --pattern "myscan_*.rhw"
#       Run a custom pattern.
#
#   bash examples/run_pdd_scan.sh --skip-helios
#       Re-extract metrics from already-completed runs, without re-running.
#
#   bash examples/run_pdd_scan.sh --commit
#       Same as default, plus git-commit and push the scan CSV at the end.
# ============================================================================
set -u
set -o pipefail

HELIOS_DIR="${HELIOS_DIR:-$HOME/Codes/Prism/Helios_11.0.0}"
HELIOS_BIN="${HELIOS_BIN:-$HELIOS_DIR/Helios.app/Contents/MacOS/Helios}"
HELIOS_ARCH_PREFIX="${HELIOS_ARCH_PREFIX:-}"
SCAN_DIR="${SCAN_DIR:-$HOME/Sims/Xcimer/Olson_PDD/PDD_scan}"
REPO_DIR="${REPO_DIR:-$HOME/helios_postprocess}"
PATTERN="${PATTERN:-wf_fl012_*.rhw}"
CSV_OUT="${CSV_OUT:-$REPO_DIR/notebooks/pdd_scan_results.csv}"

DO_COMMIT=0
SKIP_HELIOS=0
SKIP_METRICS=0
while [ $# -gt 0 ]; do
    case "$1" in
        --commit)          DO_COMMIT=1; shift ;;
        --skip-helios)     SKIP_HELIOS=1; shift ;;
        --skip-metrics)    SKIP_METRICS=1; shift ;;
        --pattern)         PATTERN="$2"; shift 2 ;;
        --scan-dir)        SCAN_DIR="$2"; shift 2 ;;
        -h|--help)         sed -n '2,30p' "$0"; exit 0 ;;
        *)                 echo "Unknown arg: $1"; exit 2 ;;
    esac
done

if [ ! -d "$SCAN_DIR" ]; then
    echo "ERROR: SCAN_DIR not found: $SCAN_DIR" >&2; exit 1
fi
if [ $SKIP_HELIOS -eq 0 ] && [ ! -x "$HELIOS_BIN" ]; then
    echo "ERROR: Helios binary not found: $HELIOS_BIN" >&2; exit 1
fi

shopt -s nullglob
RHWS=("$SCAN_DIR"/$PATTERN)
shopt -u nullglob
if [ ${#RHWS[@]} -eq 0 ]; then
    echo "No .rhw files matching '$PATTERN' under $SCAN_DIR"
    exit 1
fi

echo "PDD scan driver"
echo "  Helios:   $HELIOS_BIN"
echo "  Scan dir: $SCAN_DIR"
echo "  Pattern:  $PATTERN  (${#RHWS[@]} files)"
echo "  CSV:      $CSV_OUT"
echo

run_one_helios () {
    local rhw="$1"
    local stem
    stem=$(basename "$rhw" .rhw)
    local outdir="$SCAN_DIR/$stem"
    local exo_nested="$outdir/$stem/$stem.exo"
    local exo_flat="$outdir/$stem.exo"

    if [ -f "$exo_nested" ] || [ -f "$exo_flat" ]; then
        echo "  $stem: .exo already exists -- skipping Helios"
        return
    fi
    mkdir -p "$outdir"
    echo "  $stem: running Helios"
    ( cd "$HELIOS_DIR" && $HELIOS_ARCH_PREFIX "$HELIOS_BIN" -b -i "$rhw" -d "$outdir" -o "$stem" -x )
    local rc=$?
    if [ $rc -ne 0 ]; then
        echo "  $stem: Helios FAILED (rc=$rc)"
    fi
}

if [ $SKIP_HELIOS -eq 0 ]; then
    echo "===== Step 1: Helios runs ====="
    for rhw in "${RHWS[@]}"; do
        run_one_helios "$rhw"
    done
    echo
fi

if [ $SKIP_METRICS -eq 0 ]; then
    echo "===== Step 2: extract foam-region metrics ====="
    # Build the list of completed run base paths (handle nested vs flat layouts)
    BASES=()
    for rhw in "${RHWS[@]}"; do
        stem=$(basename "$rhw" .rhw)
        outdir="$SCAN_DIR/$stem"
        if [ -f "$outdir/$stem/$stem.exo" ]; then
            BASES+=("$outdir/$stem/$stem")
        elif [ -f "$outdir/$stem.exo" ]; then
            BASES+=("$outdir/$stem")
        fi
    done
    if [ ${#BASES[@]} -eq 0 ]; then
        echo "  No completed runs found."
    else
        ( cd "$REPO_DIR" && PYTHONPATH=. python3 examples/dump_burn_region_density.py \
            --csv "$CSV_OUT" "${BASES[@]}" )
    fi
    echo
fi

if [ $DO_COMMIT -eq 1 ]; then
    echo "===== Step 3: commit + push CSV ====="
    cd "$REPO_DIR"
    git add notebooks/pdd_scan_results.csv
    if git diff --cached --quiet; then
        echo "  No CSV changes to commit."
    else
        git commit -m "pdd_scan: results $(date +%Y-%m-%d)"
        git push
    fi
fi

echo "Done."
