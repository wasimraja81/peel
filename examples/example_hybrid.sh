#!/bin/bash
# Example: Inpaint compact sources from a radio observation
#
# This example uses the hybrid detection mode with Gaussian low-pass
# filtering and 20 iterations to get clean inpainting.

FITS_FILE="${1:?Usage: $0 <input.fits>}"

python ../peel.py "$FITS_FILE" \
  --detection-mode hybrid \
  --sigma-threshold 4.5 \
  --global-sigma-threshold 5.0 \
  --psf-radius 5 \
  --cutoff 0.10 \
  --lowpass-mode gaussian \
  --iterations 20 \
  --detect-smooth-sigma 1.0 \
  --seed-with-bright \
  --mask-output "${FITS_FILE%.fits}_mask.fits" \
  -o "${FITS_FILE%.fits}_inpainted.fits"

echo "Done! Outputs:"
echo "  Inpainted image: ${FITS_FILE%.fits}_inpainted.fits"
echo "  Detection mask: ${FITS_FILE%.fits}_mask.fits"
