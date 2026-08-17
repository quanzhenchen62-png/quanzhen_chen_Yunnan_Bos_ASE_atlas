#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash 00_run_rfmix.sh \
    --reference-bcf data/rfmix_inputs/chr1.reference.clean.bcf \
    --query-bcf data/rfmix_inputs/chr1.query.clean.bcf \
    --sample-map data/rfmix_inputs/reference_sample_map.tsv \
    --genetic-map data/rfmix_inputs/chr1.uniform_fullchr_1Mb_1cM.map.tsv \
    --chromosome 1 \
    --output-prefix results/example/rfmix/chr1.rfmix.fullchr

Required arguments:
  --reference-bcf   phased reference-panel BCF for one chromosome
  --query-bcf       phased query-panel BCF for the same chromosome
  --sample-map      two-column RFMix reference sample map
  --genetic-map     three-column genetic map used by RFMix
  --chromosome      chromosome identifier passed to RFMix
  --output-prefix   output prefix for RFMix result files

Optional arguments:
  --rfmix-bin       RFMix executable name or path [rfmix]
  --bcftools-bin    bcftools executable name or path [bcftools]
  --threads         number of threads [4]
  --generations     generations since admixture [8]
  --random-seed     RFMix random seed [12345]
EOF
}

RFMIX_BIN="rfmix"
BCFTOOLS_BIN="bcftools"
THREADS="4"
GENERATIONS="8"
RANDOM_SEED="12345"
WEIGHT="3"

REFERENCE_BCF=""
QUERY_BCF=""
SAMPLE_MAP=""
GENETIC_MAP=""
CHROMOSOME=""
OUTPUT_PREFIX=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --reference-bcf) REFERENCE_BCF="$2"; shift 2 ;;
    --query-bcf) QUERY_BCF="$2"; shift 2 ;;
    --sample-map) SAMPLE_MAP="$2"; shift 2 ;;
    --genetic-map) GENETIC_MAP="$2"; shift 2 ;;
    --chromosome) CHROMOSOME="$2"; shift 2 ;;
    --output-prefix) OUTPUT_PREFIX="$2"; shift 2 ;;
    --rfmix-bin) RFMIX_BIN="$2"; shift 2 ;;
    --bcftools-bin) BCFTOOLS_BIN="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --generations) GENERATIONS="$2"; shift 2 ;;
    --random-seed) RANDOM_SEED="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

for value_name in REFERENCE_BCF QUERY_BCF SAMPLE_MAP GENETIC_MAP CHROMOSOME OUTPUT_PREFIX; do
  if [[ -z "${!value_name}" ]]; then
    echo "Missing required argument: ${value_name}" >&2
    usage >&2
    exit 1
  fi
done

for file_path in "$REFERENCE_BCF" "$QUERY_BCF" "$SAMPLE_MAP" "$GENETIC_MAP"; do
  if [[ ! -s "$file_path" ]]; then
    echo "Missing or empty input: $file_path" >&2
    exit 1
  fi
done

if ! command -v "$RFMIX_BIN" >/dev/null 2>&1; then
  echo "RFMix executable not found: $RFMIX_BIN" >&2
  exit 1
fi

if ! command -v "$BCFTOOLS_BIN" >/dev/null 2>&1; then
  echo "bcftools executable not found: $BCFTOOLS_BIN" >&2
  exit 1
fi

mkdir -p "$(dirname "$OUTPUT_PREFIX")"

if [[ ! -s "${REFERENCE_BCF}.csi" && ! -s "${REFERENCE_BCF}.tbi" ]]; then
  "$BCFTOOLS_BIN" index -f "$REFERENCE_BCF"
fi

if [[ ! -s "${QUERY_BCF}.csi" && ! -s "${QUERY_BCF}.tbi" ]]; then
  "$BCFTOOLS_BIN" index -f "$QUERY_BCF"
fi

ref_nsites="$("$BCFTOOLS_BIN" index -n "$REFERENCE_BCF")"
query_nsites="$("$BCFTOOLS_BIN" index -n "$QUERY_BCF")"

if [[ "$ref_nsites" != "$query_nsites" ]]; then
  echo "Reference/query site numbers differ: $ref_nsites vs $query_nsites" >&2
  exit 1
fi

stdout_log="${OUTPUT_PREFIX}.rfmix.stdout.log"
stderr_log="${OUTPUT_PREFIX}.rfmix.stderr.log"

rm -f \
  "${OUTPUT_PREFIX}.fb.tsv" \
  "${OUTPUT_PREFIX}.msp.tsv" \
  "${OUTPUT_PREFIX}.sis.tsv" \
  "${OUTPUT_PREFIX}.rfmix.Q" \
  "$stdout_log" \
  "$stderr_log"

"$RFMIX_BIN" \
  -f "$QUERY_BCF" \
  -r "$REFERENCE_BCF" \
  -m "$SAMPLE_MAP" \
  -g "$GENETIC_MAP" \
  -o "$OUTPUT_PREFIX" \
  -w "$WEIGHT" \
  -G "$GENERATIONS" \
  --chromosome="$CHROMOSOME" \
  --n-threads="$THREADS" \
  --random-seed="$RANDOM_SEED" \
  > "$stdout_log" \
  2> "$stderr_log"

for expected in \
  "${OUTPUT_PREFIX}.fb.tsv" \
  "${OUTPUT_PREFIX}.msp.tsv" \
  "${OUTPUT_PREFIX}.sis.tsv" \
  "${OUTPUT_PREFIX}.rfmix.Q"
do
  if [[ ! -s "$expected" ]]; then
    echo "Expected RFMix output missing or empty: $expected" >&2
    exit 1
  fi
done

if grep -Eiq 'Initial analysis[[:space:]]*-[[:space:]]*logl[[:space:]]+[-+]?(inf|nan)' "$stderr_log"; then
  echo "RFMix produced a non-finite log-likelihood" >&2
  exit 1
fi

if ! grep -Eq 'Initial analysis[[:space:]]*-[[:space:]]*logl[[:space:]]+[-+]?[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?' "$stderr_log"; then
  echo "Could not confirm a finite Initial analysis log-likelihood in $stderr_log" >&2
  exit 1
fi

echo "RFMix completed successfully:"
echo "  output prefix: $OUTPUT_PREFIX"
echo "  MSP: ${OUTPUT_PREFIX}.msp.tsv"
echo "  FB:  ${OUTPUT_PREFIX}.fb.tsv"
