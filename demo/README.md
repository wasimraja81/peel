# Demo Quickstart

This directory contains runnable demos for peel.

## Prerequisites

Install dependencies in your active environment:

pip install -r requirements.txt

## Demo 1: Hybrid inpainting (recommended)

Run from repository root:

./demo/run_hybrid_demo.sh /path/to/input.fits

Outputs:
- <input>_hybrid_mask.fits
- <input>_hybrid_inpainted.fits

## Demo 2: Compare 10 vs 20 iterations

Run from repository root:

./demo/compare_iterations.sh /path/to/input.fits

Outputs:
- <input>_iter10_inpainted.fits
- <input>_iter20_inpainted.fits

## Optional: Visual comparison in DS9

ds9 /path/to/input.fits /path/to/input_hybrid_inpainted.fits
