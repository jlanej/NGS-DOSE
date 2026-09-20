#!/usr/bin/env bash
#SBATCH --job-name=ngsdose_cohort
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/cohort_%j.out
# Counts -> per-sample estimates -> cohort calibration -> covariate adjustment -> transmission.
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}/config.sh"
OUT="$WORK_DIR/cohort_$MODE"; mkdir -p "$OUT"

ngsdose_py ngsdose estimate "$COUNTS_DIR"/*.json.gz -r "$BUNDLE" -o "$EST_DIR" -t "$OUT/single_sample.tsv" -j "${SLURM_CPUS_PER_TASK:-8}"
ngsdose_py ngsdose cohort "$EST_DIR"/*.estimate.json.gz -r "$BUNDLE" -t "$OUT/cohort.tsv" --save-efficiencies "$OUT/efficiencies.json"

# candidate estimators of 45S dosage, from the one used in the literature to the calibrated one
# ... and three traits whose answer is known in advance, as controls of the transmission analysis
# itself: held-out autosomal sequence (no variance but error), and the EBV and mitochondrial
# loads of the cell line (large variance, none of it transmitted through the nuclear genome)
COLS="rDNA45S.18S.flat rDNA45S.cn_all_flat rDNA45S.18S rDNA45S.cn_single rDNA45S.cn rDNA5S.cn DJ.cn truth.auto chrEBV.copies chrM.copies"
# Two adjustments: NGS-PCA's coverage PCs, and the internal control-region PCs from `cohort`. How many
# PCs: by default those above the Marchenko-Pastur edge of the noise bulk (N_PC=mp; a number
# overrides it), chosen from the SVD's own spectrum...
N_PC="${N_PC:-mp}"
ngsdose_py ngsdose adjust "$OUT/cohort.tsv" --pcs "$NGSPCA_DIR/svd.pcs.txt" --n-pc "$N_PC" -c $COLS -o "$OUT/cohort.adj_ngspca.tsv"
ngsdose_py ngsdose adjust "$OUT/cohort.tsv" --n-pc "$N_PC" -c $COLS -o "$OUT/cohort.adj_ctrlpc.tsv"
# ... and the evidence to keep or overrule that default: for every number of PCs, the cross-validated
# error of the known truths (held-out autosomal, chrX, chrY, DJ) and the transmission reliability
# of the classes with its paired difference from no adjustment
SWEEP="rDNA45S.18S.flat rDNA45S.cn_single rDNA45S.cn rDNA5S.cn"
ngsdose_py ngsdose pcsweep "$OUT/cohort.tsv" --pcs "$NGSPCA_DIR/svd.pcs.txt" -c $SWEEP -p "$PEDIGREE" -o "$OUT/pcsweep.ngspca.tsv" 2> >(tee "$OUT/pcsweep.ngspca.summary.txt" >&2)
ngsdose_py ngsdose pcsweep "$OUT/cohort.tsv" -c $SWEEP -p "$PEDIGREE" -o "$OUT/pcsweep.ctrlpc.tsv" 2> >(tee "$OUT/pcsweep.ctrlpc.summary.txt" >&2)

# the ripeness table: reliability of every estimator with bootstrap CIs, and the paired
# bootstrap of each against the literature's 18S depth ratio
for tab in cohort cohort.adj_ngspca cohort.adj_ctrlpc; do
  cols="$COLS"; [ "$tab" = cohort ] || cols="$(for c in $COLS; do printf '%s.adj ' "$c"; done)"
  base="rDNA45S.18S.flat"; [ "$tab" = cohort ] || base="rDNA45S.18S.flat.adj"
  ngsdose_py ngsdose trios "$OUT/$tab.tsv" -p "$PEDIGREE" -c $cols --compare-to "$base" --json "$OUT/transmission.$tab.json" | tee "$OUT/transmission.$tab.tsv"
done
