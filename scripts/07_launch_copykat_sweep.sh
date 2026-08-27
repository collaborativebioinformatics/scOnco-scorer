#!/usr/bin/env bash

set -u
set -o pipefail

WORK="$HOME/scOnco"
IMAGE="sconco-copykat:1.1"
REMOTE="Group9_2026:/results/expression_callers/copykat"

mkdir -p "$WORK/logs"
mkdir -p "$WORK/results/copykat"

echo "=============================================="
echo "CopyKAT LLU sweep"
echo "Started: $(date)"
echo "=============================================="

fractions=(01 05 10 20 40 80)

for p in "${fractions[@]}"; do
  for r in $(seq -w 1 10); do

    dataset="LLU_p${p}_r${r}"
    outdir="$WORK/results/copykat/$dataset"
    logfile="$WORK/logs/${dataset}_copykat.log"
    remote_dataset="$REMOTE/$dataset"

    echo
    echo "----------------------------------------------"
    echo "Dataset: $dataset"
    echo "Time: $(date)"
    echo "----------------------------------------------"

    # Already safely stored remotely?
    if dx ls "$remote_dataset/" 2>/dev/null | grep -qx "DONE"; then
      echo "Remote DONE found."
      echo "COMPLETED + ALREADY UPLOADED: $dataset"
      continue
    fi

    # Run only when not completed locally
    if [[ ! -f "$outdir/DONE" ]]; then

      docker run --rm \
        -v "$WORK:/work" \
        "$IMAGE" \
        Rscript /work/code/06_run_copykat_dataset.R "$dataset" \
        > "$logfile" 2>&1

      run_status=$?

      tail -20 "$logfile"

      if [[ $run_status -ne 0 ]]; then
        echo "ERROR: CopyKAT failed for $dataset"
        echo "See $logfile"
        continue
      fi
    else
      echo "Local DONE found; computation skipped."
    fi

    if [[ ! -f "$outdir/DONE" ]]; then
      echo "ERROR: no local DONE marker for $dataset"
      continue
    fi

    dx mkdir -p "$remote_dataset"

    echo "Uploading $dataset..."

    upload_failed=0

    # Upload results except DONE
    for f in "$outdir"/*; do
      [[ "$(basename "$f")" == "DONE" ]] && continue

      dx upload "$f" \
        --path "$remote_dataset/" \
        --brief || {
          echo "ERROR uploading $f"
          upload_failed=1
          break
        }
    done

    if [[ $upload_failed -ne 0 ]]; then
      echo "UPLOAD FAILED: $dataset"
      continue
    fi

    # Upload log
    if [[ -f "$logfile" ]]; then
      dx upload "$logfile" \
        --path "$remote_dataset/" \
        --brief || {
          echo "ERROR uploading log for $dataset"
          continue
        }
    fi

    # DONE uploaded last
    dx upload "$outdir/DONE" \
      --path "$remote_dataset/" \
      --brief || {
        echo "ERROR uploading DONE for $dataset"
        continue
      }

    echo "COMPLETED + UPLOADED: $dataset"

  done
done

echo
echo "=============================================="
echo "CopyKAT sweep launcher finished"
echo "Finished: $(date)"
echo "=============================================="
