# Visualization Tools

## `generate_comparison_pngs.py`

Generate publication-ready comparison PNG visualizations for FITS inpainting results.

### Usage

```bash
python tools/generate_comparison_pngs.py original.fits inpainted.fits [output_prefix]
```

### Arguments

- `original` (required): Path to original FITS file
- `inpainted` (required): Path to inpainted/peeled FITS file
- `-o, --output` (optional): Output prefix for PNG files. If not specified, uses the basename of the original FITS file.
- `--dpi` (optional): DPI for comparison panel images (default: 220)

### Output Files

The script generates up to 5 PNG files:

1. **`{prefix}_original.png`**: Single-panel grayscale view of original image
2. **`{prefix}_peeled.png`**: Single-panel grayscale view of inpainted image
3. **`{prefix}_residual.png`**: Single-panel viridis view of residual (Original - Peeled)
4. **`{prefix}_compare_panel.png`**: Three-panel comparison (Original, Peeled, Residual) in grayscale
5. **`{prefix}_compare_panel_viridis.png`**: Three-panel comparison (Original, Peeled, Residual) in viridis colormap

### Features

- **Smart NaN handling**: Renders NaN/invalid pixels as white in all panels
- **Matched color scales**: Original and peeled images use the same 1-99 percentile normalization
- **Data-driven residual scaling**: Residual uses its own 1-99 percentile range for optimal contrast
- **Consistent masking**: NaN regions from the original FITS are masked as white in all panels, including the peeled image
- **asinh stretch**: Applies nonlinear asinh stretch for enhanced dynamic range visibility

### Example

After running the peel inpainting tool:

```bash
python peel.py myimage.fits --hybrid

# Generate comparison PNGs
python tools/generate_comparison_pngs.py myimage.fits myimage_inpainted_hybrid.fits myimage_comparison
```

This creates files like:
- `myimage_comparison_original.png`
- `myimage_comparison_peeled.png`
- `myimage_comparison_residual.png`
- `myimage_comparison_compare_panel.png`
- `myimage_comparison_compare_panel_viridis.png`
