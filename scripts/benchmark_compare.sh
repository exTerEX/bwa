#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

REFERENCE=""
READS=""

usage() {
    cat <<'EOF'
Usage:
  ./scripts/benchmark_compare.sh --reference <ref.fasta|fastq> --reads <reads.fasta|fastq> [options]

Options:
  --max-edits <int>        Edit budget for inexact benchmark (default: 2)
  --max-hits <int>         Max hits per read passed to aligner (default: 50)
  --repeat <int>           Number of repeated runs per benchmark case (default: 5)
  --sa-sample-rate <int>   Sampling rate for sampled index (default: 8)
  --seed-candidates <int>  Seed candidate cap for seed mode (default: 8192)
    --csv-out <path>         Write benchmark rows to CSV file
    --csv-append             Append to existing CSV file instead of overwriting
  --help                   Show this message

Environment variable equivalents:
    MAX_EDITS, MAX_HITS, REPEAT, SA_SAMPLE_RATE, SEED_CANDIDATES, CSV_OUT
EOF
}

MAX_EDITS="${MAX_EDITS:-2}"
MAX_HITS="${MAX_HITS:-50}"
REPEAT="${REPEAT:-5}"
SA_SAMPLE_RATE="${SA_SAMPLE_RATE:-8}"
SEED_CANDIDATES="${SEED_CANDIDATES:-8192}"
CSV_OUT="${CSV_OUT:-}"
CSV_APPEND=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --reference)
            REFERENCE="$2"
            shift 2
            ;;
        --reads)
            READS="$2"
            shift 2
            ;;
        --max-edits)
            MAX_EDITS="$2"
            shift 2
            ;;
        --max-hits)
            MAX_HITS="$2"
            shift 2
            ;;
        --repeat)
            REPEAT="$2"
            shift 2
            ;;
        --sa-sample-rate)
            SA_SAMPLE_RATE="$2"
            shift 2
            ;;
        --seed-candidates)
            SEED_CANDIDATES="$2"
            shift 2
            ;;
        --csv-out)
            CSV_OUT="$2"
            shift 2
            ;;
        --csv-append)
            CSV_APPEND=1
            shift
            ;;
        --help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage
            exit 1
            ;;
    esac

done

if [[ -z "$REFERENCE" || -z "$READS" ]]; then
    echo "Both --reference and --reads are required." >&2
    usage
    exit 1
fi

if [[ ! -f "$REFERENCE" ]]; then
    echo "Reference file not found: $REFERENCE" >&2
    exit 1
fi
if [[ ! -f "$READS" ]]; then
    echo "Reads file not found: $READS" >&2
    exit 1
fi

if [[ "$REPEAT" -lt 1 ]]; then
    echo "REPEAT must be >= 1" >&2
    exit 1
fi
if [[ "$SA_SAMPLE_RATE" -lt 1 ]]; then
    echo "SA_SAMPLE_RATE must be >= 1" >&2
    exit 1
fi

if [[ -n "$CSV_OUT" ]]; then
    mkdir -p "$(dirname "$CSV_OUT")"
fi

make >/dev/null

TMP_DIR="$(mktemp -d)"
cleanup() {
    rm -rf "$TMP_DIR"
}
trap cleanup EXIT

FULL_INDEX="$TMP_DIR/full.bwaidx"
SAMPLED_INDEX="$TMP_DIR/sampled.bwaidx"

./build/bwa index --reference "$REFERENCE" --output "$FULL_INDEX" --sa-sample-rate 1 >/dev/null
./build/bwa index --reference "$REFERENCE" --output "$SAMPLED_INDEX" --sa-sample-rate "$SA_SAMPLE_RATE" >/dev/null

csv_quote() {
    local value="$1"
    value="${value//\"/\"\"}"
    printf '"%s"' "$value"
}

csv_write_header_if_needed() {
    if [[ -z "$CSV_OUT" ]]; then
        return
    fi

    if [[ "$CSV_APPEND" -eq 0 ]]; then
        : >"$CSV_OUT"
    fi

    if [[ ! -s "$CSV_OUT" ]]; then
        echo 'section,label,index_kind,search_mode,max_edits,repeat,max_hits,sa_sample_rate,seed_candidates,total_ms,avg_ms,total_hits,reference,reads' >>"$CSV_OUT"
    fi
}

csv_append_row() {
    if [[ -z "$CSV_OUT" ]]; then
        return
    fi

    local section="$1"
    local label="$2"
    local index_kind="$3"
    local search_mode="$4"
    local max_edits="$5"
    local total_ms="$6"
    local avg_ms="$7"
    local total_hits="$8"

    {
        csv_quote "$section"
        printf ','
        csv_quote "$label"
        printf ','
        csv_quote "$index_kind"
        printf ','
        csv_quote "$search_mode"
        printf ',%s,%s,%s,%s,%s,%s,%s,%s,' \
            "$max_edits" "$REPEAT" "$MAX_HITS" "$SA_SAMPLE_RATE" "$SEED_CANDIDATES" "$total_ms" "$avg_ms" "$total_hits"
        csv_quote "$REFERENCE"
        printf ','
        csv_quote "$READS"
        printf '\n'
    } >>"$CSV_OUT"
}

run_case() {
    local section="$1"
    local label="$2"
    local index_kind="$3"
    local mode="$4"
    local max_edits="$5"
    local index_path="$6"

    local start_ns end_ns elapsed_ns
    local out=""

    start_ns="$(date +%s%N)"
    for ((i = 0; i < REPEAT; ++i)); do
        out="$(./build/bwa align \
            --index "$index_path" \
            --reads "$READS" \
            --max-edits "$max_edits" \
            --max-hits "$MAX_HITS" \
            --search-mode "$mode" \
            --seed-candidates "$SEED_CANDIDATES")"
    done
    end_ns="$(date +%s%N)"

    elapsed_ns="$((end_ns - start_ns))"

    local total_hits avg_ms total_ms
    total_hits="$(awk '{
        for (i = 1; i <= NF; ++i) {
            split($i, kv, "=")
            if (kv[1] == "hits") {
                sum += kv[2]
            }
        }
    } END { print sum + 0 }' <<<"$out")"

    total_ms="$(awk -v ns="$elapsed_ns" 'BEGIN { printf "%.3f", ns / 1000000.0 }')"
    avg_ms="$(awk -v ns="$elapsed_ns" -v r="$REPEAT" 'BEGIN { printf "%.3f", (ns / 1000000.0) / r }')"

    printf "%-26s total_ms=%8s avg_ms=%8s total_hits=%s\n" "$label" "$total_ms" "$avg_ms" "$total_hits"

    csv_append_row "$section" "$label" "$index_kind" "$mode" "$max_edits" "$total_ms" "$avg_ms" "$total_hits"
}

csv_write_header_if_needed

echo "Benchmark settings: repeat=$REPEAT max_edits=$MAX_EDITS max_hits=$MAX_HITS sa_sample_rate=$SA_SAMPLE_RATE seed_candidates=$SEED_CANDIDATES"
echo "Reference=$REFERENCE"
echo "Reads=$READS"
echo

echo "Inexact search comparison (backtracking vs seed)"
run_case "inexact" "full_sa + backtracking" "full_sa" "backtracking" "$MAX_EDITS" "$FULL_INDEX"
run_case "inexact" "full_sa + seed" "full_sa" "seed" "$MAX_EDITS" "$FULL_INDEX"
run_case "inexact" "sampled_sa + backtracking" "sampled_sa" "backtracking" "$MAX_EDITS" "$SAMPLED_INDEX"
run_case "inexact" "sampled_sa + seed" "sampled_sa" "seed" "$MAX_EDITS" "$SAMPLED_INDEX"
echo

echo "Exact locate latency comparison (full SA vs sampled SA)"
run_case "exact_locate" "full_sa + exact locate" "full_sa" "auto" 0 "$FULL_INDEX"
run_case "exact_locate" "sampled_sa + exact locate" "sampled_sa" "auto" 0 "$SAMPLED_INDEX"

if [[ -n "$CSV_OUT" ]]; then
    echo
    echo "CSV written to $CSV_OUT"
fi
