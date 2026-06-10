#!/usr/bin/env python3
"""Peel compact sources from FITS images using iterative Fourier inpainting.

This module detects compact sources, builds a mask, and fills masked pixels
through iterative band-limited interpolation while preserving measured pixels.
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

try:
    from astropy.io import fits
except ImportError as exc:
    raise SystemExit(
        "astropy is required. Install with: pip install astropy"
    ) from exc

try:
    from scipy.ndimage import gaussian_filter
except ImportError:
    gaussian_filter = None


def robust_stats(values: np.ndarray) -> tuple[float, float]:
    """Return robust central value and scale estimate for an array.

    Parameters
    ----------
    values : np.ndarray
        Input image/data array that may contain NaNs or infs.

    Returns
    -------
    tuple[float, float]
        `(median, sigma)` where sigma is MAD-based (`1.4826 * MAD`) and falls
        back to standard deviation if needed.
    """
    finite = np.isfinite(values)
    if not np.any(finite):
        return 0.0, 1.0

    valid = values[finite]
    median = float(np.median(valid))
    mad = float(np.median(np.abs(valid - median)))
    sigma = 1.4826 * mad

    if not np.isfinite(sigma) or sigma <= 0:
        sigma = float(np.std(valid))
    if not np.isfinite(sigma) or sigma <= 0:
        sigma = 1.0

    return median, sigma


def local_maxima_3x3(image: np.ndarray) -> np.ndarray:
    """Find local maxima in a 2D image using a 3x3 neighborhood.

    Parameters
    ----------
    image : np.ndarray
        2D image used for peak detection.

    Returns
    -------
    np.ndarray
        Boolean mask where `True` indicates a pixel greater than or equal to
        all 8 neighbors.
    """
    padded = np.pad(image, 1, mode="edge")
    center = padded[1:-1, 1:-1]

    maxima = np.ones_like(image, dtype=bool)
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            if dy == 0 and dx == 0:
                continue
            neighbor = padded[1 + dy : 1 + dy + image.shape[0], 1 + dx : 1 + dx + image.shape[1]]
            maxima &= center >= neighbor

    return maxima


def build_disk_offsets(radius: int) -> list[tuple[int, int]]:
    """Build `(dy, dx)` offsets inside a filled disk of given radius."""
    offsets: list[tuple[int, int]] = []
    for dy in range(-radius, radius + 1):
        for dx in range(-radius, radius + 1):
            if dy * dy + dx * dx <= radius * radius:
                offsets.append((dy, dx))
    return offsets


def expand_peaks_to_psf_mask(peaks: np.ndarray, radius: int) -> np.ndarray:
    """Expand seed pixels into a circular (disk) mask.

    Parameters
    ----------
    peaks : np.ndarray
        Boolean seed mask.
    radius : int
        Radius in pixels for circular expansion.

    Returns
    -------
    np.ndarray
        Expanded boolean mask.
    """
    if radius <= 0:
        return peaks.copy()

    mask = np.zeros_like(peaks, dtype=bool)
    offsets = build_disk_offsets(radius)

    for dy, dx in offsets:
        src_y0 = max(0, -dy)
        src_y1 = peaks.shape[0] - max(0, dy)
        src_x0 = max(0, -dx)
        src_x1 = peaks.shape[1] - max(0, dx)

        dst_y0 = max(0, dy)
        dst_y1 = dst_y0 + (src_y1 - src_y0)
        dst_x0 = max(0, dx)
        dst_x1 = dst_x0 + (src_x1 - src_x0)

        mask[dst_y0:dst_y1, dst_x0:dst_x1] |= peaks[src_y0:src_y1, src_x0:src_x1]

    return mask


def expand_peaks_to_box_mask(peaks: np.ndarray, box_size: int) -> np.ndarray:
    """Expand seed pixels into a square mask of size `box_size`.

    Parameters
    ----------
    peaks : np.ndarray
        Boolean seed mask.
    box_size : int
        Side length in pixels for square expansion.

    Returns
    -------
    np.ndarray
        Expanded boolean mask.
    """
    if box_size <= 1:
        return peaks.copy()

    low = -(box_size // 2)
    high = low + box_size - 1

    mask = np.zeros_like(peaks, dtype=bool)

    for dy in range(low, high + 1):
        for dx in range(low, high + 1):
            src_y0 = max(0, -dy)
            src_y1 = peaks.shape[0] - max(0, dy)
            src_x0 = max(0, -dx)
            src_x1 = peaks.shape[1] - max(0, dx)

            dst_y0 = max(0, dy)
            dst_y1 = dst_y0 + (src_y1 - src_y0)
            dst_x0 = max(0, dx)
            dst_x1 = dst_x0 + (src_x1 - src_x0)

            mask[dst_y0:dst_y1, dst_x0:dst_x1] |= peaks[src_y0:src_y1, src_x0:src_x1]

    return mask


def detect_point_source_mask(
    image: np.ndarray,
    sigma_threshold: float,
    psf_radius_px: int,
    box_size_px: int,
    detection_mode: str,
    global_sigma_threshold: float,
    local_bg_sigma: float,
    local_noise_sigma: float,
    detect_smooth_sigma: float,
    seed_with_bright: bool,
) -> np.ndarray:
    """Detect compact sources and return mask pixels to inpaint.

    Detection supports `global`, `local`, and `hybrid` modes. The resulting
    seed map is expanded either with circular PSF radius or square box size.

    Parameters
    ----------
    image : np.ndarray
        Input 2D image.
    sigma_threshold : float
        Local/global detection threshold in sigma units.
    psf_radius_px : int
        Circular expansion radius in pixels.
    box_size_px : int
        Square expansion size in pixels; if >1, overrides radius expansion.
    detection_mode : str
        One of `global`, `local`, or `hybrid`.
    global_sigma_threshold : float
        Global-threshold sigma used in hybrid mode.
    local_bg_sigma : float
        Gaussian sigma for local background model.
    local_noise_sigma : float
        Gaussian sigma for local noise model.
    detect_smooth_sigma : float
        Optional pre-detection smoothing sigma.
    seed_with_bright : bool
        If true, use all significant bright pixels as seeds; otherwise peaks.

    Returns
    -------
    np.ndarray
        Boolean source mask restricted to finite input pixels.
    """
    finite = np.isfinite(image)
    med, sigma = robust_stats(image)

    filled = image.copy()
    filled[~finite] = med

    if detection_mode in ("local", "hybrid"):
        if gaussian_filter is None:
            raise SystemExit(
                "scipy is required for --detection-mode local. Install with: pip install scipy"
            )
        if local_bg_sigma <= 0 or local_noise_sigma <= 0:
            raise ValueError("--local-bg-sigma and --local-noise-sigma must be > 0 in local mode")

        background = gaussian_filter(filled, sigma=local_bg_sigma, mode="reflect")
        residual = filled - background

        detection_image = residual
        if detect_smooth_sigma > 0:
            detection_image = gaussian_filter(detection_image, sigma=detect_smooth_sigma, mode="reflect")

        local_scale = 1.4826 * gaussian_filter(np.abs(residual), sigma=local_noise_sigma, mode="reflect")
        finite_scale = local_scale[finite]
        scale_floor = float(np.median(finite_scale)) if finite_scale.size else sigma
        if not np.isfinite(scale_floor) or scale_floor <= 0:
            scale_floor = sigma if sigma > 0 else 1.0
        local_scale = np.maximum(local_scale, scale_floor)

        bright = detection_image > sigma_threshold * local_scale
        local_peaks = bright & local_maxima_3x3(detection_image)
        local_seeds = bright if seed_with_bright else local_peaks

        if detection_mode == "local":
            seeds = local_seeds
        else:
            global_detection = filled
            global_bright = global_detection > med + global_sigma_threshold * sigma
            global_peaks = global_bright & local_maxima_3x3(global_detection)
            seeds = local_seeds | global_peaks
    else:
        detection_image = filled
        if detect_smooth_sigma > 0:
            if gaussian_filter is None:
                raise SystemExit(
                    "scipy is required for --detect-smooth-sigma > 0. Install with: pip install scipy"
                )
            detection_image = gaussian_filter(detection_image, sigma=detect_smooth_sigma, mode="reflect")

        bright = detection_image > med + sigma_threshold * sigma
        peaks = bright & local_maxima_3x3(detection_image)
        seeds = bright if seed_with_bright else peaks

    if box_size_px > 1:
        psf_mask = expand_peaks_to_box_mask(seeds, box_size_px)
    else:
        psf_mask = expand_peaks_to_psf_mask(seeds, psf_radius_px)
    psf_mask &= finite

    return psf_mask


def build_lowpass_filter(
    shape: tuple[int, int],
    cutoff_cyc_per_pix: float,
    mode: str,
) -> np.ndarray:
    """Build a 2D frequency-domain low-pass filter.

    Parameters
    ----------
    shape : tuple[int, int]
        Image shape `(ny, nx)`.
    cutoff_cyc_per_pix : float
        Cutoff scale in cycles per pixel.
    mode : str
        `hard` for binary cutoff, `gaussian` for smooth taper.

    Returns
    -------
    np.ndarray
        2D filter array aligned with `fftshift`-centered spectra.
    """
    h, w = shape
    fy = np.fft.fftshift(np.fft.fftfreq(h))
    fx = np.fft.fftshift(np.fft.fftfreq(w))
    yy, xx = np.meshgrid(fy, fx, indexing="ij")
    radius = np.sqrt(yy * yy + xx * xx)

    if mode == "hard":
        return radius <= cutoff_cyc_per_pix
    if mode == "gaussian":
        return np.exp(-0.5 * (radius / cutoff_cyc_per_pix) ** 2)

    raise ValueError(f"Unknown lowpass mode: {mode}")


def bandlimited_inpaint(
    original: np.ndarray,
    mask_to_fill: np.ndarray,
    cutoff_cyc_per_pix: float,
    lowpass_mode: str,
    iterations: int,
) -> np.ndarray:
    """Fill masked pixels by iterative band-limited Fourier interpolation.

    The method enforces two constraints each iteration:
    1) spectral smoothness via low-pass filtering in Fourier space,
    2) exact preservation of measured (unmasked) finite pixels.

    Parameters
    ----------
    original : np.ndarray
        Input image, potentially with NaNs.
    mask_to_fill : np.ndarray
        Boolean mask of pixels to inpaint.
    cutoff_cyc_per_pix : float
        Low-pass cutoff/scale in cycles per pixel.
    lowpass_mode : str
        `hard` or `gaussian` low-pass behavior.
    iterations : int
        Number of FFT/IFFT refinement iterations.

    Returns
    -------
    np.ndarray
        Inpainted image estimate.
    """
    finite = np.isfinite(original)
    med, _ = robust_stats(original)

    original_filled = original.copy()
    original_filled[~finite] = med

    current = original_filled.copy()
    current[mask_to_fill] = 0.0
    lowpass = build_lowpass_filter(
        current.shape,
        cutoff_cyc_per_pix=cutoff_cyc_per_pix,
        mode=lowpass_mode,
    )

    measured_mask = finite & ~mask_to_fill

    for _ in range(iterations):
        work = current.copy()
        spectrum = np.fft.fftshift(np.fft.fft2(work))
        filtered_spectrum = spectrum * lowpass
        interpolated = np.fft.ifft2(np.fft.ifftshift(filtered_spectrum)).real

        current[mask_to_fill] = interpolated[mask_to_fill]
        current[measured_mask] = original_filled[measured_mask]

    return current


def load_fits_2d(path: Path, hdu_index: int) -> tuple[np.ndarray, fits.Header]:
    """Load and validate a 2D FITS image from the selected HDU.

    Parameters
    ----------
    path : Path
        Path to FITS file.
    hdu_index : int
        HDU index to read.

    Returns
    -------
    tuple[np.ndarray, fits.Header]
        2D float64 image array and original HDU header.
    """
    with fits.open(path) as hdul:
        data = hdul[hdu_index].data
        header = hdul[hdu_index].header

    if data is None:
        raise ValueError(f"HDU {hdu_index} has no data")

    arr = np.asarray(data, dtype=np.float64)
    arr = np.squeeze(arr)
    if arr.ndim != 2:
        raise ValueError(f"Expected 2D image after squeeze; got shape {arr.shape}")

    return arr, header


def save_fits(path: Path, data: np.ndarray, header: fits.Header | None = None) -> None:
    """Write array data to a FITS file, overwriting if it exists."""
    hdu = fits.PrimaryHDU(data=data, header=header)
    hdu.writeto(path, overwrite=True)


def parse_args() -> argparse.Namespace:
    """Parse and validate command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed, validated runtime options for source detection and inpainting.
    """
    parser = argparse.ArgumentParser(
        description="Detect point sources, mask them, and fill via iterative band-limited Fourier interpolation."
    )
    parser.add_argument("input_fits", type=Path, help="Input FITS image path")
    parser.add_argument(
        "-o",
        "--output-fits",
        type=Path,
        default=None,
        help="Output FITS path (default: <input>_inpainted.fits)",
    )
    parser.add_argument("--hdu", type=int, default=0, help="HDU index to read (default: 0)")
    parser.add_argument(
        "--sigma-threshold",
        type=float,
        default=5.0,
        help="Detection threshold in robust sigma above median (default: 5.0)",
    )
    parser.add_argument(
        "--psf-radius",
        type=int,
        default=2,
        help="Radius (pixels) around each peak to mask (default: 2)",
    )
    parser.add_argument(
        "--box-size",
        type=int,
        default=0,
        help=(
            "If >1, use a square box of this size (pixels) around each peak. "
            "Example: 10 gives a 10x10 mask region. Overrides --psf-radius."
        ),
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=0.12,
        help=(
            "Low-pass scale in cycles/pixel (0 < cutoff <= 0.5). "
            "For hard mode it is a radius cutoff; for gaussian mode it is the "
            "frequency sigma. Lower = smoother interpolation (default: 0.12)"
        ),
    )
    parser.add_argument(
        "--lowpass-mode",
        choices=("hard", "gaussian"),
        default="gaussian",
        help="Frequency low-pass shape (default: gaussian)",
    )
    parser.add_argument(
        "--iterations",
        type=int,
        default=8,
        help="Number of Fourier interpolation iterations (default: 8)",
    )
    parser.add_argument(
        "--include-nans",
        action="store_true",
        help="Also fill original NaN pixels in addition to point-source masks",
    )
    parser.add_argument(
        "--detection-mode",
        choices=("global", "local", "hybrid"),
        default="global",
        help="Source detection mode (default: global)",
    )
    parser.add_argument(
        "--global-sigma-threshold",
        type=float,
        default=5.0,
        help=(
            "Global peak threshold (robust sigma) used only in hybrid mode to recover "
            "small compact sources (default: 5.0)"
        ),
    )
    parser.add_argument(
        "--local-bg-sigma",
        type=float,
        default=20.0,
        help="Gaussian sigma (px) for local background estimate in local mode (default: 20.0)",
    )
    parser.add_argument(
        "--local-noise-sigma",
        type=float,
        default=8.0,
        help="Gaussian sigma (px) for local noise map in local mode (default: 8.0)",
    )
    parser.add_argument(
        "--detect-smooth-sigma",
        type=float,
        default=0.0,
        help="Optional Gaussian smoothing sigma (px) before thresholding (default: 0.0)",
    )
    parser.add_argument(
        "--seed-with-bright",
        action="store_true",
        help=(
            "Use all significant bright pixels as seeds (better for larger compact sources). "
            "Default uses only local-maximum peaks."
        ),
    )
    parser.add_argument(
        "--mask-output",
        type=Path,
        default=None,
        help="Optional FITS path to save the binary mask used for filling",
    )

    args = parser.parse_args()

    if args.cutoff <= 0 or args.cutoff > 0.5:
        raise ValueError("--cutoff must be in (0, 0.5]")
    if args.iterations < 1:
        raise ValueError("--iterations must be >= 1")
    if args.sigma_threshold <= 0:
        raise ValueError("--sigma-threshold must be > 0")
    if args.psf_radius < 0:
        raise ValueError("--psf-radius must be >= 0")
    if args.box_size < 0:
        raise ValueError("--box-size must be >= 0")
    if args.local_bg_sigma <= 0:
        raise ValueError("--local-bg-sigma must be > 0")
    if args.local_noise_sigma <= 0:
        raise ValueError("--local-noise-sigma must be > 0")
    if args.detect_smooth_sigma < 0:
        raise ValueError("--detect-smooth-sigma must be >= 0")
    if args.global_sigma_threshold <= 0:
        raise ValueError("--global-sigma-threshold must be > 0")

    return args


def main() -> None:
    """Run CLI workflow: load image, detect mask, inpaint, and write outputs."""
    args = parse_args()

    image, header = load_fits_2d(args.input_fits, args.hdu)
    source_mask = detect_point_source_mask(
        image=image,
        sigma_threshold=args.sigma_threshold,
        psf_radius_px=args.psf_radius,
        box_size_px=args.box_size,
        detection_mode=args.detection_mode,
        global_sigma_threshold=args.global_sigma_threshold,
        local_bg_sigma=args.local_bg_sigma,
        local_noise_sigma=args.local_noise_sigma,
        detect_smooth_sigma=args.detect_smooth_sigma,
        seed_with_bright=args.seed_with_bright,
    )

    mask_to_fill = source_mask.copy()
    if args.include_nans:
        mask_to_fill |= ~np.isfinite(image)

    inpainted = bandlimited_inpaint(
        original=image,
        mask_to_fill=mask_to_fill,
        cutoff_cyc_per_pix=args.cutoff,
        lowpass_mode=args.lowpass_mode,
        iterations=args.iterations,
    )

    output_path = args.output_fits
    if output_path is None:
        output_path = args.input_fits.with_name(args.input_fits.stem + "_inpainted.fits")

    final = image.copy()
    final[mask_to_fill] = inpainted[mask_to_fill]

    # Record provenance in header
    invocation = " ".join(sys.argv)
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    header.add_history(f"peel.py run on {timestamp}")
    # Split long invocation across multiple HISTORY cards (72 chars each)
    chunk = 72
    for i in range(0, len(invocation), chunk):
        header.add_history(invocation[i : i + chunk])

    save_fits(output_path, final, header)

    if args.mask_output is not None:
        save_fits(args.mask_output, mask_to_fill.astype(np.uint8), None)

    print(f"Input image shape: {image.shape}")
    print(f"Detection mode: {args.detection_mode}")
    print(f"Low-pass mode: {args.lowpass_mode}")
    print(f"Detected source pixels in mask: {int(source_mask.sum())}")
    print(f"Total filled pixels: {int(mask_to_fill.sum())}")
    print(f"Output written to: {output_path}")
    if args.mask_output is not None:
        print(f"Mask written to: {args.mask_output}")


if __name__ == "__main__":
    main()
