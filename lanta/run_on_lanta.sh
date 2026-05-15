#!/usr/bin/env bash
###############################################################################
# run_on_lanta.sh — example invocation against the 226-sample hatchery cohort
#
# Run on LANTA login node (no SLURM job needed; ngsTracts is CPU-light).
###############################################################################
set -euo pipefail

BASE=/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04
STAGE3_DIR=${BASE}/ngsPedigree/stage3_out             # produced by ngsPedigree
FAI=${BASE}/00-samples/fClaHyb_Gar_LG.fa.fai          # 28-LG ref index
OUT=${BASE}/ngsTracts/out

# Optional inputs
INV_ATLAS=${BASE}/MS_Inversions/inversion_atlas.tsv   # if present, used by Stage 3
KARYOTYPES=${BASE}/MS_Inversions/parent_karyotypes.tsv # for inversion-aware NCO rule

# Sanity: required inputs
[[ -d "${STAGE3_DIR}" ]] || { echo "[FAIL] Stage 3 outputs not found: ${STAGE3_DIR}"; exit 1; }
[[ -s "${STAGE3_DIR}/departure_intervals.tsv" ]] || { echo "[FAIL] departure_intervals.tsv missing"; exit 1; }
[[ -s "${STAGE3_DIR}/stage3.args.tsv" ]] || { echo "[FAIL] stage3.args.tsv missing"; exit 1; }
[[ -s "${FAI}" ]] || { echo "[FAIL] FAI missing: ${FAI}"; exit 1; }

mkdir -p "${OUT}"

# Activate the same env ngsPedigree uses (has numpy + pandas)
source ~/.bashrc
mamba activate assembly

# Main classifier
python3 "$(dirname "$0")/../scripts/STEP_TRC_01_classify_intervals.py" \
    --stage3-dir "${STAGE3_DIR}" \
    --fai        "${FAI}" \
    --outdir     "${OUT}" \
    ${KARYOTYPES:+--parent-karyotypes "${KARYOTYPES}"}

# Optional: refine CO breakpoints if per-site file was emitted
if [[ -s "${STAGE3_DIR}/per_site_haplotype_calls.tsv" ]]; then
  echo ""
  echo "[$(date '+%F %T')] per_site_haplotype_calls.tsv present — running TRC_02"
  python3 "$(dirname "$0")/../scripts/STEP_TRC_02_traversal_scan.py" \
      --tract-classifications "${OUT}/tract_classifications.tsv" \
      --per-site-file         "${STAGE3_DIR}/per_site_haplotype_calls.tsv" \
      --outdir                "${OUT}"
else
  echo ""
  echo "[$(date '+%F %T')] per_site_haplotype_calls.tsv not present — skipping TRC_02"
  echo "  (run ngsPedigree Stage 3 with --emit-per-site to enable breakpoint refinement)"
fi

echo ""
echo "[$(date '+%F %T')] [DONE] ngsTracts outputs in ${OUT}/"
ls -lh "${OUT}/"
