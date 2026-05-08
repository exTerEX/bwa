#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

make >/dev/null

./build/bwa index --reference data/example_reference.fasta --output build/example.bwaidx >/dev/null

exact_output="$(./build/bwa align --index build/example.bwaidx --reads data/example_reads.fasta --max-edits 0 --max-hits 10)"
if ! grep -Eq "read=read_exact .*hits=[1-9][0-9]*" <<<"$exact_output"; then
    echo "Expected read_exact to have at least one exact hit"
    exit 1
fi

inexact_output="$(./build/bwa align --index build/example.bwaidx --reads data/example_reads.fasta --max-edits 1 --max-hits 10)"
if ! grep -Eq "read=read_mismatch .*hits=[1-9][0-9]*" <<<"$inexact_output"; then
    echo "Expected read_mismatch to align with one edit"
    exit 1
fi
if ! grep -Eq "read=read_indel .*hits=[1-9][0-9]*" <<<"$inexact_output"; then
    echo "Expected read_indel to align with one edit"
    exit 1
fi

echo "Example test passed"
echo "$inexact_output"
