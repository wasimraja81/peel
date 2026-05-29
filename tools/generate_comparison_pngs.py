#!/usr/bin/env python3
"""
Generate publication-ready comparison PNG visualizations for FITS inpainting results.

This script creates comparison plots showing original, peeled (inpainted), and residual
images with matched color scales and white NaN masking.

Usage:
    python generate_comparison_pngs.py original.fits inpainted.fits [output_prefix]

Example:
    python generate_comparison_pngs.py myimage.fits myimage_inpainted.fits myimage
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.visualization import AsinhStretch, ImageNormalize


def render_with_white_nans(ax, data: np.ndarray, norm, cmap):
    """Render data with NaN regions as white, using the given norm and colormap."""
    nan_mask = ~np.isfinite(data)
    
    # Replace NaN with a safe value for normalization
    safe_data = np.where(nan_mask, np.nanmin(data[~nan_mask]) if (~nan_mask).any() else 0, data)
    
    # Normalize the safe data
    normalized = norm(safe_data)
    
    # Get RGBA from colormap
    rgba = cmap(normalized)
    
    # Set NaN pixels to white
    rgba[nan_mask] = [1.0, 1.0, 1.0, 1.0]
    
    # Display the RGBA array
    im = ax.imshow(rgba, origin="lower", interpolation="nearest")
    
    # Create ScalarMappable for colorbar
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    return sm


def render_with_external_mask(ax, data: np.ndarray, norm, cmap, external_mask: np.ndarray):
    """Render data with an external NaN mask applied, rendering masked regions as white."""
    # Replace masked regions with a safe value for normalization
    safe_data = np.where(external_mask, np.nanmin(data[~external_mask]) if (~external_mask).any() else 0, data)
    
    # Normalize the safe data
    normalized = norm(safe_data)
    
    # Get RGBA from colormap
    rgba = cmap(normalized)
    
    # Set masked pixels to white
    rgba[external_mask] = [1.0, 1.0, 1.0, 1.0]
    
    # Display the RGBA array
    im = ax.imshow(rgba, origin="lower", interpolation="nearest")
    
    # Create ScalarMappable for colorbar
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    return sm


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate comparison PNG visualizations for FITS inpainting results."
    )
    parser.add_argument(
        "original",
        type=str,
        help="Path to original FITS file",
    )
    parser.add_argument(
        "inpainted",
        type=str,
        help="Path to inpainted FITS file",
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        default=None,
        help="Output prefix for PNG files (default: basename of original FITS)",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=220,
        help="DPI for comparison panel images (default: 220)",
    )

    args = parser.parse_args()

    orig_path = Path(args.original)
    inpainted_path = Path(args.inpainted)

    if not orig_path.exists():
        print(f"Error: Original FITS file not found: {orig_path}")
        return
    if not inpainted_path.exists():
        print(f"Error: Inpainted FITS file not found: {inpainted_path}")
        return

    # Determine output prefix
    if args.output:
        out_prefix = args.output
    else:
        out_prefix = orig_path.stem

    # Load data
    orig = np.squeeze(fits.getdata(orig_path)).astype(float)
    inpainted = np.squeeze(fits.getdata(inpainted_path)).astype(float)
    res = orig - inpainted

    # Compute statistics
    finite_orig = np.isfinite(orig)
    finite_res = np.isfinite(res)
    orig_nan_mask = ~finite_orig

    if finite_orig.any():
        disp_lo, disp_hi = np.nanpercentile(orig[finite_orig], [1.0, 99.0])
    else:
        disp_lo, disp_hi = -1.0, 1.0

    if finite_res.any():
        res_lo, res_hi = np.nanpercentile(res[finite_res], [1.0, 99.0])
    else:
        res_lo, res_hi = -1.0, 1.0

    display_norm = ImageNormalize(
        vmin=disp_lo,
        vmax=disp_hi,
        stretch=AsinhStretch(a=0.03),
        clip=True,
    )

    residual_norm = ImageNormalize(
        vmin=res_lo,
        vmax=res_hi,
        stretch=AsinhStretch(a=0.03),
        clip=True,
    )

    # Setup colormaps
    gray_cmap = plt.colormaps["gray"].copy()
    gray_cmap.set_bad(color="white", alpha=1.0)

    viridis_cmap = plt.colormaps["viridis"].copy()
    viridis_cmap.set_bad(color="white", alpha=1.0)

    residual_cmap = plt.colormaps["viridis"].copy()
    residual_cmap.set_bad(color="white", alpha=1.0)

    # Single-panel images (grayscale)
    plt.figure(figsize=(8, 8))
    ax = plt.gca()
    sm = render_with_white_nans(ax, orig, display_norm, gray_cmap)
    plt.title("Original")
    plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_original.png", dpi=200)
    plt.close()

    plt.figure(figsize=(8, 8))
    ax = plt.gca()
    sm = render_with_external_mask(ax, inpainted, display_norm, gray_cmap, orig_nan_mask)
    plt.title("Peeled (hybrid)")
    plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_peeled.png", dpi=200)
    plt.close()

    plt.figure(figsize=(8, 8))
    ax = plt.gca()
    sm = render_with_external_mask(ax, res, residual_norm, residual_cmap, orig_nan_mask)
    plt.title("Residual (Orig - Hybrid)")
    plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_residual.png", dpi=200)
    plt.close()

    # Comparison panel (grayscale)
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), constrained_layout=True)
    render_with_white_nans(axes[0], orig, display_norm, gray_cmap)
    axes[0].set_title("Original")

    render_with_external_mask(axes[1], inpainted, display_norm, gray_cmap, orig_nan_mask)
    axes[1].set_title("Peeled (hybrid)")

    sm = render_with_external_mask(axes[2], res, residual_norm, residual_cmap, orig_nan_mask)
    axes[2].set_title("Residual (Orig - Hybrid)")

    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])

    fig.colorbar(sm, ax=axes[2], fraction=0.046, pad=0.04)
    fig.savefig(f"{out_prefix}_compare_panel.png", dpi=args.dpi)
    plt.close(fig)

    # Comparison panel (viridis)
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), constrained_layout=True)
    render_with_white_nans(axes[0], orig, display_norm, viridis_cmap)
    axes[0].set_title("Original")

    render_with_external_mask(axes[1], inpainted, display_norm, viridis_cmap, orig_nan_mask)
    axes[1].set_title("Peeled (hybrid)")

    sm = render_with_external_mask(axes[2], res, residual_norm, residual_cmap, orig_nan_mask)
    axes[2].set_title("Residual (Orig - Hybrid)")

    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])

    fig.colorbar(sm, ax=axes[2], fraction=0.046, pad=0.04)
    fig.savefig(f"{out_prefix}_compare_panel_viridis.png", dpi=args.dpi)
    plt.close(fig)

    print(f"Generated comparison PNGs with prefix: {out_prefix}")
    print(f"  - {out_prefix}_original.png")
    print(f"  - {out_prefix}_peeled.png")
    print(f"  - {out_prefix}_residual.png")
    print(f"  - {out_prefix}_compare_panel.png")
    print(f"  - {out_prefix}_compare_panel_viridis.png")


if __name__ == "__main__":
    main()
