#!/bin/bash
# Demo: Full workflow including visualization
# Usage: ./demo/complete_workflow.sh input.fits

set -e

if [ $# -eq 0 ]; then
    echo "Usage: $0 input.fits"
    exit 1
fi

INPUT=$1
BASENAME=$(basename "$INPUT" .fits)

echo "=== PEEL Inpainting + Visualization Workflow ==="
echo "Input: $INPUT"
echo ""

# Step 1: Run hybrid inpainting
echo "Step 1: Running hybrid inpainting..."
python peel.py "$INPUT" \
    --detection-mode hybrid \
    --sigma-threshold 4.5 \
    --global-sigma-threshold 5.0 \
    --psf-radius 5 \
    --cutoff 0.10 \
    --lowpass-mode gaussian \
    --iterations 20 \
    --seed-with-bright \
    --mask-output "${BASENAME}_mask.fits" \
    -o "${BASENAME}_inpainted.fits"

echo "✓ Inpainting complete"
echo "  - ${BASENAME}_inpainted.fits (inpainted image)"
echo "  - ${BASENAME}_mask.fits (source mask)"
echo ""

# Step 2: Generate comparison PNGs
echo "Step 2: Generating publication-ready PNG comparisons..."
python tools/generate_comparison_pngs.py \
    "$INPUT" \
    "${BASENAME}_inpainted.fits" \
    -o "$BASENAME"

echo "✓ Visualization complete"
echo "  - ${BASENAME}_compare_panel.png (grayscale comparison)"
echo "  - ${BASENAME}_compare_panel_viridis.png (viridis comparison)"
echo "  - ${BASENAME}_original.png"
echo "  - ${BASENAME}_peeled.png"
echo "  - ${BASENAME}_residual.png"
echo ""

echo "=== Workflow complete ==="
