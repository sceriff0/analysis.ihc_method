#!/usr/bin/env bash
# =============================================================================
# cibersortx_fractions.sh — CIBERSORTx Fractions on the study's bulk RNA
#
# The molecular page (analysis/_children/molecular_body.Rmd, chunk
# deconv-cibersortx) writes the mixture this job reads and reads back what this
# job writes, so the two directories below are the contract:
#   data/cibersortx/input/mixture.txt      gene x sample, GeneSymbol first  (page -> job)
#   data/cibersortx/output/*_Results.txt   fractions; *_Adjusted.txt with  (job -> page)
#                                          B-mode batch correction
#
# The token never touches the repo: it lives in ~/.cibersortx_token (one line, the
# token alone). The username is the e-mail that token was issued to — the two are
# bound together — and is set below. Tokens expire; request a new one at
# https://cibersortx.stanford.edu/getoken.php when the run says so.
#
# Settings, and why (README_CIBERSORTxFractions.txt, v1.0):
#   --sigmatrix LM22.txt    the 22-type leukocyte signature (microarray-derived)
#   --QN FALSE              quantile normalisation is for microarray mixtures;
#                           this is RNA-seq
#   --rmbatchBmode TRUE     B-mode corrects the RNA-seq mixture against the
#                           microarray signature — the case LM22 + RNA-seq is
#   --perm 100              permutations for the per-sample p-value (0 = none)
#
# Run from the project root on the cluster:
#   sbatch code/cibersortx_fractions.sh        # or: bash code/cibersortx_fractions.sh
# then re-knit analysis/molecular_massimo2.Rmd.
# =============================================================================
#SBATCH --job-name=cibersortx
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
set -euo pipefail

# The project root. Under sbatch the script runs from Slurm's spool copy, so its
# own path says nothing; the submit directory does (submit from the project root,
# as the header says). Outside Slurm, walk up from the script's own location.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  ROOT="$SLURM_SUBMIT_DIR"
else
  ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
ROOT="${IHC_ROOT:-$ROOT}"                             # or name it outright
[[ -f "$ROOT/code/cibersortx_fractions.sh" ]] || { echo "$ROOT is not the ihc_method root — sbatch from the project root, or export IHC_ROOT" >&2; exit 1; }
INPUT_DIR="$ROOT/data/cibersortx/input"
OUTPUT_DIR="$ROOT/data/cibersortx/output"
MIXTURE="mixture.txt"
SIGMATRIX="${CIBERSORTX_SIGMATRIX:-LM22.txt}"       # must sit in INPUT_DIR
SIF="${CIBERSORTX_SIF:-$HOME/containers/cibersortx_fractions.sif}"
# Build the image once with:  singularity pull "$SIF" docker://cibersortx/fractions

# The account the token was issued to. Set in the script so the job is one
# command; an exported CIBERSORTX_USERNAME still overrides it.
CIBERSORTX_USERNAME="${CIBERSORTX_USERNAME:-mohammadreza.javadinamin@ieo.it}"
[[ -s "$HOME/.cibersortx_token" ]] || { echo "no token in ~/.cibersortx_token — paste the token from https://cibersortx.stanford.edu/getoken.php there" >&2; exit 1; }
CIBERSORTX_TOKEN="$(cat "$HOME/.cibersortx_token")"

[[ -f "$INPUT_DIR/$MIXTURE"   ]] || { echo "no $INPUT_DIR/$MIXTURE — knit the molecular page first" >&2; exit 1; }
[[ -f "$INPUT_DIR/$SIGMATRIX" ]] || { echo "no $INPUT_DIR/$SIGMATRIX — copy LM22.txt there" >&2; exit 1; }
[[ -f "$SIF" ]]                  || { echo "no image at $SIF — singularity pull docker://cibersortx/fractions" >&2; exit 1; }
mkdir -p "$OUTPUT_DIR"

echo "[INFO] CIBERSORTx Fractions: $MIXTURE vs $SIGMATRIX (QN off, B-mode on)"
singularity exec \
  --bind "$INPUT_DIR:/src/data" \
  --bind "$OUTPUT_DIR:/src/outdir" \
  "$SIF" \
  /src/CIBERSORTxFractions \
    --username "$CIBERSORTX_USERNAME" \
    --token "$CIBERSORTX_TOKEN" \
    --mixture "$MIXTURE" \
    --sigmatrix "$SIGMATRIX" \
    --QN FALSE \
    --rmbatchBmode TRUE \
    --perm 100 \
    --verbose TRUE
echo "[DONE] output in $OUTPUT_DIR:"; ls -1 "$OUTPUT_DIR"
