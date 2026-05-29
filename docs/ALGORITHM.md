# peel — Algorithm Details

## Problem Statement

Given a 2D astronomical image with compact point sources (and possibly NaN pixels), inpaint the regions where sources were detected such that:
- The interpolated values are statistically consistent with the surrounding noise
- High-frequency Fourier components are suppressed (band-limited interpolation)
- The final image appears as if sources were never present

## Key Concepts

### 1. Detection via Robust Statistics

**Global Mode** (`--detection-mode global`):
- Compute image median $m$ and MAD $\sigma = 1.4826 \times \text{MAD}$
- Find pixels $p > m + \sigma_{\text{thresh}} \times \sigma$
- Refine to local maxima (3×3 neighborhood)

**Local Mode** (`--detection-mode local`):
- Estimate spatially-varying background: $b(x,y) = \text{Gauss}(I; \sigma_{bg})$
- Compute residual: $r(x,y) = I(x,y) - b(x,y)$
- Estimate spatially-varying noise: $\sigma_{local}(x,y) = 1.4826 \times \text{MAD}(\text{Gauss}(|r|; \sigma_{noise}))$
- Detect: pixels where $r(x,y) > \sigma_{\text{thresh}} \times \sigma_{local}(x,y)$

**Hybrid Mode** (`--detection-mode hybrid`):
- Union of local and global results
- Recovers both large extended sources (local) and small bright peaks (global)

### 2. Seed Expansion

For each detected pixel, expand a masked region:

**Circular** (`--psf-radius R`):
- Mask all pixels within distance $R$ of the seed: a filled disk

**Square** (`--box-size B`):
- Mask a $B \times B$ region centered on the seed

**Seed Selection** (`--seed-with-bright`):
- Without: only local peaks become seeds (one per source)
- With: all significantly bright pixels become seeds (expands source footprint)

The final mask is the union of all expanded regions.

### 3. Iterative Band-Limited Inpainting

The Gerchberg-Papoulis algorithm solves:
- Constraint 1 (spectral): Image must be band-limited to frequencies $< f_c$
- Constraint 2 (data): Unmasked pixels must equal measured values

**Iteration $n$**:

1. **Zero masked pixels**:
   $$I_{\text{work}}(x,y) = \begin{cases} I_n(x,y) & \text{if unmeasured} \\ 0 & \text{if masked} \end{cases}$$
   
   Or, on later iterations, seed with previous IFFT:
   $$I_{\text{work}}(x,y) = \begin{cases} I_n(x,y) & \text{if unmeasured} \\ I_{n}^{\text{ifft}}(x,y) & \text{if masked} \end{cases}$$

2. **FFT**:
   $$\hat{I} = \text{FFT}(I_{\text{work}})$$

3. **Low-pass filter**:
   $$\hat{I}_{\text{filtered}} = \hat{I} \cdot H(f)$$
   
   where $H(f)$ is either:
   - **Hard cutoff**: $H(f) = \begin{cases} 1 & |f| \le f_c \\ 0 & |f| > f_c \end{cases}$
   - **Gaussian**: $H(f) = \exp\left(-\frac{1}{2}(f/f_c)^2\right)$

4. **IFFT**:
   $$I_{\text{interp}} = \text{IFFT}(\hat{I}_{\text{filtered}})$$

5. **Restore data constraint**:
   $$I_{n+1}(x,y) = \begin{cases} I_{\text{measured}}(x,y) & \text{if unmeasured} \\ I_{\text{interp}}(x,y) & \text{if masked} \end{cases}$$

Repeat until convergence (typically 10–20 iterations).

### 4. Convergence

The RMS change in masked-pixel values between successive iterations decreases exponentially:
$$\Delta_n = \sqrt{\frac{1}{N_{\text{mask}}} \sum_{(x,y) \in \text{mask}} (I_{n+1} - I_n)^2} \propto \lambda^n$$

where $\lambda < 1$ depends on the frequency cutoff and spectral content. By $n \sim 20$, convergence is typically excellent.

## Parameter Choices

### Detection Sensitivity

- **`--sigma-threshold` (local mode)**: Lower = more sensitive. Risk of including diffuse emission.
  - Recommended: 3.5–5.0
  
- **`--local-bg-sigma`**: Scale of background variation. Larger = smoother background, better for extended emission.
  - Recommended: 15–30 px

- **`--local-noise-sigma`**: Scale of noise fluctuations. Smaller = finer noise structure, more sensitivity to faint sources.
  - Recommended: 6–12 px

### Masking Footprint

- **`--psf-radius`**: Should cover ~1–2 times the PSF FWHM
  - For PSF FWHM = 10 px: `--psf-radius 5–6`

- **`--seed-with-bright`**: Trade-off
  - On: Larger footprints per source, better for extended sources, risk of overmasking
  - Off: More conservative, only peaks, better for isolated bright sources

### Interpolation Quality

- **`--cutoff` (cycles/pixel)**: Lower = smoother, less high-frequency noise, fewer artefacts
  - Recommended: 0.08–0.15
  
- **`--lowpass-mode`**: Gaussian recommended to avoid ringing
  
- **`--iterations`**: More = better fidelity but slower
  - Recommended: 8–20
  - Check convergence with per-iteration RMS output (add with `--verbose` flag if implemented)

## NaN Handling

- **Without `--include-nans`**: NaN pixels are filled with median during FFT (to avoid corruption) but remain NaN in output
- **With `--include-nans`**: NaN pixels are also inpainted and replaced with interpolated values

⚠️ Large NaN regions (e.g., > 10% of image) can corrupt the FFT; use without caution.

## Failure Modes

1. **Overmasking**: Too low `--sigma-threshold` with `--seed-with-bright` → extended emission gets included
   - Fix: Raise `--sigma-threshold`, use peaks only, or inspect mask before running

2. **Undermasking**: Faint sources missed
   - Fix: Lower `--sigma-threshold`, reduce `--local-noise-sigma`, use `--detect-smooth-sigma 1.0`

3. **Ringing**: Sharp artefacts at mask boundaries
   - Fix: Use `--lowpass-mode gaussian` and/or lower `--cutoff` further

4. **Slow convergence**: Few iterations chosen, visible discontinuities
   - Fix: Increase `--iterations` to 20–30
