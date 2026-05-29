#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <input.fits>"
  exit 1
fi

input_fits="$1"
base="${input_fits%.fits}"

python peel.py "$input_fits" \
  --detection-mode hybrid \
  --sigma-threshold 4.5 \
  --global-sigma-threshold 5.0 \
  --local-bg-sigma 20 \
  --local-noise-sigma 8 \
  --detect-smooth-sigma 1.0 \
  --seed-with-bright \
  --psf-radius 5 \
  --cutoff 0.10 \
  --lowpass-mode gaussian \
  --iterations 20 \
  --mask-output "${base}_hybrid_mask.fits" \
  -o "${base}_hybrid_inpainted.fits"

echo "Wrote: ${base}_hybrid_mask.fits"
echo "Wrote: ${base}_hybrid_inpainted.fits"
