#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

# Configuration: allow override of interpreter commands
PYTHON=${PYTHON:-python3}
QUARTO=${QUARTO:-quarto}

# Usage help
usage() {
  cat <<EOF
Usage: ${0##*/} [-h|--help]

Run the Sanger sequencing processing pipeline:
  1. Convert AB1 files to FASTA (raw & filtered) and generate quality plots
  2. Align forward/reverse reads and build consensus sequences
  3. Merge HVSI and HVSII consensus sequences
  4. Generate an HTML QC report

Environment variables:
  PYTHON: path to Python interpreter (default: python3)
  QUARTO: path to quarto binary (default: quarto)
EOF
}

if [[ "${1-}" =~ ^(-h|--help)$ ]]; then
  usage
  exit 0
fi

# Ensure dependencies are installed
command -v mafft >/dev/null 2>&1 || { echo "Error: mafft is not installed. Exiting." >&2; exit 1; }
command -v "$PYTHON" >/dev/null 2>&1 || { echo "Error: $PYTHON is not available. Exiting." >&2; exit 1; }
"$PYTHON" - <<'PYTHONTST' >/dev/null 2>&1
import Bio
PYTHONTST
if [[ $? -ne 0 ]]; then
  echo "Error: Biopython is not installed in $PYTHON environment. Exiting." >&2
  exit 1
fi
command -v "$QUARTO" >/dev/null 2>&1 || { echo "Error: quarto is not installed. Exiting." >&2; exit 1; }

# Create required directories
dirs=(fasta filtered consensus plots aligned final logs)
for d in "${dirs[@]}"; do
  mkdir -p "$d"
done

#### Step 1: Convert AB1 to FASTA with Phred filtering ####
echo "Step 1: Converting AB1 to FASTA with quality filtering..."
: > logs/conversion.log
ab1_files=( *.ab1 )
if [[ ${#ab1_files[@]} -eq 0 ]]; then
  echo "Warning: No .ab1 files found in current directory." | tee -a logs/conversion.log
else
  for ab1 in "${ab1_files[@]}"; do
    base="${ab1%.ab1}"
    "$PYTHON" 01_convert_ab1_quality.py \
      "$ab1" \
      "fasta/${base}.fasta" \
      "filtered/${base}_filtered.fasta" \
      "plots/${base}_quality.png" \
      >> logs/conversion.log 2>&1
    echo "Processed $ab1" >> logs/conversion.log
  done
fi

#### Step 2: Align forward/reverse reads and build consensus ####
echo "Step 2: Aligning forward/reverse reads and building consensus..."
: > logs/align.log
shopt -s nullglob
for f_fasta in filtered/*-F_filtered.fasta; do
  prefix="${f_fasta##*/}"
  prefix="${prefix%-F_filtered.fasta}"
  r_fasta="filtered/${prefix}-R_filtered.fasta"
  if [[ ! -f "$r_fasta" ]]; then
    echo "Warning: Missing reverse read for $prefix. Skipping." >> logs/align.log
    continue
  fi
  r_fasta_rc="filtered/${prefix}-R_rc_filtered.fasta"
  aln_file="aligned/${prefix}_aligned.fasta"
  cons_file="consensus/${prefix}_consensus.fasta"

  # Reverse-complement reverse read
  "$PYTHON" - <<PYRC >> logs/align.log 2>&1
from Bio import SeqIO
record = SeqIO.read("$r_fasta", "fasta")
record.seq = record.seq.reverse_complement()
record.id += "_RC"
record.description = "Reverse-complemented"
SeqIO.write(record, "$r_fasta_rc", "fasta")
PYRC

  # Align reads using temporary file
  tmp=$(mktemp)
  cat "$f_fasta" "$r_fasta_rc" > "$tmp"
  mafft --auto "$tmp" > "$aln_file" 2>> logs/align.log
  rm -f "$tmp" "$r_fasta_rc"

  # Build consensus
  "$PYTHON" 02_make_consensus.py "$aln_file" "$cons_file" >> logs/align.log 2>&1
  echo "Consensus built for $prefix" >> logs/align.log
done

#### Step 3: Merge HVSI and HVSII consensus ####
echo "Step 3: Merging HVSI and HVSII consensus sequences..."
: > logs/merge.log
shopt -s nullglob
for hvsi in consensus/*-HSV1_consensus.fasta; do
  sample="${hvsi##*/}"
  sample="${sample%-HSV1_consensus.fasta}"
  hvsii="consensus/${sample}-HSV2_consensus.fasta"
  final_out="final/${sample}_merged.fasta"

  if [[ ! -f "$hvsii" ]]; then
    echo "Warning: Missing HVSI or HVSII for $sample. Skipping." >> logs/merge.log
    continue
  fi

  seq1=$(grep -v '^>' "$hvsi" | tr -d '\n')
  seq2=$(grep -v '^>' "$hvsii" | tr -d '\n')
  {
    printf ">%s_HVSI_HVSII\n" "$sample"
    printf "%s%s\n" "$seq1" "$seq2"
  } > "$final_out"
  echo "Merged $sample" >> logs/merge.log
done

#### Step 4: Generate HTML QC report ####
echo "Step 4: Generating HTML QC report..."

: > logs/qc.log
# Ensure ipykernel is installed for this environment
if ! "$PYTHON" -m pip show ipykernel >/dev/null 2>&1; then
  echo "Installing ipykernel..." | tee -a logs/qc.log
  "$PYTHON" -m pip install --quiet ipykernel >> logs/qc.log 2>&1
fi
# Ensure the 'ab1-qc' Jupyter kernel is available for Quarto execution
echo "Registering Jupyter kernel 'ab1-qc'..." | tee -a logs/qc.log
"$PYTHON" -m ipykernel install --user --name ab1-qc --display-name "AB1 QC" >> logs/qc.log 2>&1 || true
# Render the QC report, saving both output and errors
echo "Rendering QC report..." | tee -a logs/qc.log
"$QUARTO" render 03_qc_report.qmd --to html --output qc_report.html 2>&1 | tee -a logs/qc.log
echo "Report generated: qc_report.html" | tee -a logs/qc.log

echo "Pipeline complete!"