#!/usr/bin/env bash
# ============================================================================
# run_zero_d_matrix.sh
# --------------------
# End-to-end driver for the zero-D verification matrix on the Studio.
# For each .rhw input file:
#   1. Run Helios in batch mode (skip if .exo already exists)
#   2. Run examples/verify_zero_d.py to extract rates and append to CSV
#
# Usage
# -----
#   bash examples/run_zero_d_matrix.sh
#       Runs the default matrix below.
#
#   bash examples/run_zero_d_matrix.sh --commit
#       Same, plus git-commit and push the updated CSV at the end.
#
#   bash examples/run_zero_d_matrix.sh --rhw FILE --label LBL --comp PRESET
#       Run a single (filename, label, composition) triple.  Useful for
#       resuming a single failed test.
#
# Requirements
# ------------
#   - Helios v11.0.0 installed at ~/Codes/Prism/Helios_11.0.0  (path adjustable
#     via $HELIOS_DIR env var)
#   - .rhw input files already generated under $SIM_DIR (use
#     examples/make_static_rhw.py to create them).
#   - python3 with helios_postprocess on path (PYTHONPATH=.)
# ============================================================================
set -u
set -o pipefail

# ── Paths -- override via env if needed ──
HELIOS_DIR="${HELIOS_DIR:-$HOME/Codes/Prism/Helios_11.0.0}"
HELIOS_BIN="${HELIOS_BIN:-$HELIOS_DIR/Helios}"
SIM_DIR="${SIM_DIR:-$HOME/Sims/Xcimer/Olson_PDD/Fraley_burn}"
REPO_DIR="${REPO_DIR:-$HOME/helios_postprocess}"

# ── Default matrix (label, .rhw stem, composition preset) ──
# Each row: stem_in_SIM_DIR | composition_preset | label_for_CSV
MATRIX=(
    "DT_static_5keV          | pure-DT | DT_static_5keV"
    "DT_static_7keV          | pure-DT | DT_static_7keV"
    "DT_static_10keV         | pure-DT | DT_static_10keV"
    "DT_static_15keV         | pure-DT | DT_static_15keV"
    "DT_static_WF_5keV       | 1.7C-CH | foam_static_5keV"
    "DT_static_WF_7keV       | 1.7C-CH | foam_static_7keV"
    "DT_static_WF_10keV      | 1.7C-CH | foam_static_10keV"
    "DT_static_WF_15keV      | 1.7C-CH | foam_static_15keV"
    "DT_static_10keV_n1e23   | pure-DT | DT_static_10keV_n1e23"
    "DT_static_10keV_n1e24   | pure-DT | DT_static_10keV_n1e24"
    "DT_static_10keV_n1e25   | pure-DT | DT_static_10keV_n1e25"
)

# ── Arg parsing ──
SINGLE_RHW=""
SINGLE_LBL=""
SINGLE_CMP=""
DO_COMMIT=0
SKIP_HELIOS=0
SKIP_VERIFY=0
while [ $# -gt 0 ]; do
    case "$1" in
        --commit)          DO_COMMIT=1; shift ;;
        --skip-helios)     SKIP_HELIOS=1; shift ;;
        --skip-verify)     SKIP_VERIFY=1; shift ;;
        --rhw)             SINGLE_RHW="$2"; shift 2 ;;
        --label)           SINGLE_LBL="$2"; shift 2 ;;
        --comp)            SINGLE_CMP="$2"; shift 2 ;;
        -h|--help)
            sed -n '2,28p' "$0"; exit 0 ;;
        *) echo "Unknown arg: $1"; exit 2 ;;
    esac
done

# ── Sanity checks ──
if [ ! -x "$HELIOS_BIN" ] && [ $SKIP_HELIOS -eq 0 ]; then
    echo "ERROR: Helios binary not found / not executable: $HELIOS_BIN" >&2
    echo "  Override path with: HELIOS_BIN=/path/to/Helios bash $0" >&2
    exit 1
fi
if [ ! -d "$SIM_DIR" ]; then
    echo "ERROR: SIM_DIR not found: $SIM_DIR" >&2
    exit 1
fi
if [ ! -f "$REPO_DIR/examples/verify_zero_d.py" ]; then
    echo "ERROR: verify_zero_d.py not found under $REPO_DIR/examples/" >&2
    exit 1
fi

# Single-run mode overrides the matrix
if [ -n "$SINGLE_RHW" ]; then
    if [ -z "$SINGLE_LBL" ] || [ -z "$SINGLE_CMP" ]; then
        echo "ERROR: --rhw requires both --label and --comp" >&2; exit 2
    fi
    MATRIX=( "$SINGLE_RHW | $SINGLE_CMP | $SINGLE_LBL" )
fi

# ── Helper: run one .rhw ──
run_one () {
    local stem="$1"           # e.g. DT_static_5keV
    local comp="$2"           # pure-DT or 1.7C-CH
    local label="$3"          # CSV label
    local rhw="$SIM_DIR/${stem}.rhw"
    local outdir="$SIM_DIR/${stem}"
    local exo="$outdir/${stem}.exo"
    local base="$outdir/${stem}"        # without extension, for verify_zero_d

    echo
    echo "=============================================================="
    echo "  $stem    (label=$label, comp=$comp)"
    echo "=============================================================="

    if [ ! -f "$rhw" ]; then
        echo "  SKIP: RHW input not found: $rhw"
        return
    fi

    # --- Helios ---
    if [ $SKIP_HELIOS -eq 0 ]; then
        if [ -f "$exo" ]; then
            echo "  Helios: .exo already exists -- skipping run"
        else
            echo "  Helios: running ($HELIOS_BIN -b -i $rhw -d $outdir -o $stem -x)"
            mkdir -p "$outdir"
            ( cd "$HELIOS_DIR" && "./$(basename "$HELIOS_BIN")" -b -i "$rhw" -d "$outdir" -o "$stem" -x )
            local rc=$?
            if [ $rc -ne 0 ]; then
                echo "  Helios FAILED (rc=$rc) on $stem -- continuing to next"
                return
            fi
        fi
    fi

    # --- Verification ---
    if [ $SKIP_VERIFY -eq 0 ]; then
        if [ ! -f "$exo" ]; then
            echo "  verify: skipped (no .exo present)"
            return
        fi
        echo "  verify: extracting rates"
        ( cd "$REPO_DIR" && PYTHONPATH=. python3 examples/verify_zero_d.py \
            "$base" --comp "$comp" --label "$label" )
    fi
}

# ── Main loop ──
echo "Studio zero-D verification matrix"
echo "  Helios: $HELIOS_BIN"
echo "  Sims:   $SIM_DIR"
echo "  Repo:   $REPO_DIR"
echo "  Cases:  ${#MATRIX[@]}"

for row in "${MATRIX[@]}"; do
    # Strip whitespace, split on '|'
    stem=$(echo "$row" | cut -d'|' -f1 | xargs)
    comp=$(echo "$row" | cut -d'|' -f2 | xargs)
    lbl=$(echo  "$row" | cut -d'|' -f3 | xargs)
    run_one "$stem" "$comp" "$lbl"
done

# ── Optional commit + push ──
if [ $DO_COMMIT -eq 1 ]; then
    echo
    echo "=============================================================="
    echo "  Committing + pushing updated CSV"
    echo "=============================================================="
    cd "$REPO_DIR"
    git add notebooks/verification_results.csv
    if git diff --cached --quiet; then
        echo "  No CSV changes to commit."
    else
        git commit -m "verify_zero_d: matrix run $(date +%Y-%m-%d)"
        git push
    fi
fi

echo
echo "Done."
