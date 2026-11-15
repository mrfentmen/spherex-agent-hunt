"""FITS HDU1 -> detections CSV. numpy only. No astropy."""
from __future__ import annotations

import numpy as np


def background_mad(img: np.ndarray) -> tuple[float, float]:
    med = float(np.nanmedian(img))
    mad = float(np.nanmedian(np.abs(img - med)))
    return med, 1.4826 * mad


def find_peaks(img: np.ndarray, nsigma: float = 5.0, edge: int = 5) -> list[tuple[int, int, float]]:
    med, sigma = background_mad(img)
    if not np.isfinite(med):
        return []
    if not np.isfinite(sigma) or sigma <= 0:
        thresh = med
    else:
        thresh = med + nsigma * sigma
    work = np.where(np.isfinite(img), img, med)
    ys, xs = np.where(work > thresh)
    h, w = work.shape
    out = []
    for y, x in zip(ys.tolist(), xs.tolist()):
        if x < edge or y < edge or x >= w - edge or y >= h - edge:
            continue
        y0, y1 = max(0, y - 4), min(h, y + 5)
        x0, x1 = max(0, x - 4), min(w, x + 5)
        if np.isfinite(img[y0:y1, x0:x1]).mean() < 0.5:
            continue  # masked/bad region, not a source
        y0, y1 = max(0, y - 1), min(work.shape[0], y + 2)
        x0, x1 = max(0, x - 1), min(work.shape[1], x + 2)
        if work[y, x] == work[y0:y1, x0:x1].max():
            out.append((x, y, float(work[y, x])))
    # NMS: brightest-first, min Chebyshev separation 3 px; pre-cap keeps it cheap.
    out.sort(key=lambda p: p[2], reverse=True)
    out = out[:2000]
    kept: list[tuple[int, int, float]] = []
    for x, y, v in out:
        if all(max(abs(x - kx), abs(y - ky)) >= 3 for kx, ky, _ in kept):
            kept.append((x, y, v))
    return kept


def read_hdu1(path: str):
    """Return (header dict, image float32 2D array). HDU0 skipped, HDU1 read."""
    with open(path, 'rb') as f:
        for hdu in range(2):
            cards = []
            while True:
                blk = f.read(2880)
                if not blk:
                    raise ValueError('EOF in header')
                for i in range(0, 2880, 80):
                    cards.append(blk[i:i+80].decode('ascii', 'replace'))
                if any(c.startswith('END') for c in cards[-36:]):
                    break
            hdr = {c[:8].strip(): c[10:30].strip() for c in cards if '=' in c[:10]}
            naxis = int(hdr.get('NAXIS', '0').split('/')[0])
            naxes = [int(hdr.get(f'NAXIS{n}', '0').split('/')[0]) for n in range(1, naxis+1)]
            import math
            nbytes = math.ceil(abs(int(hdr.get('BITPIX', '8').split('/')[0]))//8 * (math.prod(naxes) if naxes and all(naxes) else 0) / 2880) * 2880 if naxis else 0
            if hdu == 0:
                if nbytes:
                    f.seek(nbytes, 1)
                continue
            raw = f.read(nbytes) if nbytes else b''
            n1, n2 = naxes[0], naxes[1]
            arr = np.frombuffer(raw[:n1*n2*4], dtype='>f4').reshape((n2, n1)).astype(np.float32)
            return hdr, arr
    raise ValueError('no HDU1')


def local_snr(peak: float, med: float, var_val, mad_sigma: float) -> float:
    """Per-pixel SNR: variance plane wins, MAD scatter is the fallback."""
    import math
    try:
        v = float(var_val)
    except (TypeError, ValueError):
        v = float("nan")
    if math.isfinite(v) and v > 0:
        return (peak - med) / math.sqrt(v)
    return (peak - med) / (mad_sigma or 1.0)


def read_image_and_variance(path: str):
    """Return (header, image, variance-or-None) for a SPHEREx L2B frame."""
    with open(path, 'rb') as f:
        hdr0 = img = var = None
        shape = None
        for _ in range(8):
            cards = []
            while True:
                blk = f.read(2880)
                if not blk:
                    break
                for i in range(0, 2880, 80):
                    cards.append(blk[i:i+80].decode('ascii', 'replace'))
                if any(c.startswith('END') for c in cards[-36:]):
                    break
            if not cards:
                break
            hdr = {c[:8].strip(): c[10:30].strip() for c in cards if '=' in c[:10]}
            naxis = int(hdr.get('NAXIS', '0').split('/')[0])
            naxes = [int(hdr.get(f'NAXIS{n}', '0').split('/')[0]) for n in range(1, naxis+1)]
            import math
            bpp = abs(int(hdr.get('BITPIX', '8').split('/')[0])) // 8
            nbytes = math.ceil(bpp * (math.prod(naxes) if naxes and all(naxes) else 0) / 2880) * 2880 if naxis else 0
            raw = f.read(nbytes) if nbytes else b''
            if naxis == 2 and hdr.get('XTENSION', '').strip("' ") == 'IMAGE':
                n1, n2 = naxes[0], naxes[1]
                arr = np.frombuffer(raw[:n1*n2*4], dtype='>f4').reshape((n2, n1)).astype(np.float32)
                if img is None:
                    hdr0, img, shape = hdr, arr, (n1, n2)
                elif hdr.get('EXTNAME', '').strip("' ") == 'VARIANCE' and (n1, n2) == shape:
                    var = arr
                    break
        if img is None:
            raise ValueError('no IMAGE HDU')
        return hdr0, img, var


def band_from_path(path: str) -> str:
    """Detector/band tag from a FITS file name, e.g. '1D2'. '' when absent."""
    import os
    from .sia import detector_key
    return detector_key(os.path.basename(path))


def aperture_flux(img, cx: float, cy: float, r: int = 4) -> float:
    """Fixed-aperture sum minus local background (annulus median)."""
    h, w = img.shape
    xs = slice(max(0, int(cx) - 9), min(w, int(cx) + 10))
    ys = slice(max(0, int(cy) - 9), min(h, int(cy) + 10))
    cut = img[ys, xs].astype(float)
    yy, xx = np.mgrid[ys, xs]
    dist = np.hypot(xx - cx, yy - cy)
    ap = np.where(np.isfinite(cut), cut, np.nan)
    in_ap = (dist <= r) & np.isfinite(cut)
    ring = (dist > r + 2) & (dist <= r + 5) & np.isfinite(cut)
    bg = float(np.median(cut[ring])) if ring.any() else float(np.nanmedian(ap))
    return float(np.nansum(np.where(in_ap, cut - bg, 0.0)))


def _hdr_float(hdr, key: str, default: float = 0.0) -> float:
    v = hdr.get(key)
    if v is None:
        return default
    try:
        return float(str(v).split("/")[0].strip().strip("'"))
    except ValueError:
        return default


def _sip_poly(hdr, prefix: str, order: int, u: float, v: float) -> float:
    # SIP sums run over every stored term, constant included.
    s = 0.0
    for p in range(order + 1):
        for q in range(order + 1 - p):
            c = _hdr_float(hdr, f"{prefix}_{p}_{q}", 0.0)
            if c:
                s += c * (u ** p) * (v ** q)
    return s


def sip_to_radec(hdr, x: float, y: float) -> tuple[float, float]:
    """Pixel to sky via CD/PC matrix, SIP distortion, and TAN deprojection."""
    import math
    ra0 = _hdr_float(hdr, "CRVAL1")
    dec0 = _hdr_float(hdr, "CRVAL2")
    crpix1 = _hdr_float(hdr, "CRPIX1")
    crpix2 = _hdr_float(hdr, "CRPIX2")
    cd11 = _hdr_float(hdr, "CD1_1"); cd12 = _hdr_float(hdr, "CD1_2")
    cd21 = _hdr_float(hdr, "CD2_1"); cd22 = _hdr_float(hdr, "CD2_2")
    if cd11 == 0.0 and cd12 == 0.0 and cd21 == 0.0 and cd22 == 0.0:
        pc11 = _hdr_float(hdr, "PC1_1", 1.0); pc12 = _hdr_float(hdr, "PC1_2", 0.0)
        pc21 = _hdr_float(hdr, "PC2_1", 0.0); pc22 = _hdr_float(hdr, "PC2_2", 1.0)
        cd11, cd12 = pc11 * _hdr_float(hdr, "CDELT1", 1.0), pc12 * _hdr_float(hdr, "CDELT1", 1.0)
        cd21, cd22 = pc21 * _hdr_float(hdr, "CDELT2", 1.0), pc22 * _hdr_float(hdr, "CDELT2", 1.0)
    dx, dy = x - crpix1, y - crpix2
    u = cd11 * dx + cd12 * dy
    v = cd21 * dx + cd22 * dy
    ao = int(_hdr_float(hdr, "A_ORDER", 0)); bo = int(_hdr_float(hdr, "B_ORDER", 0))
    U = u + _sip_poly(hdr, "A", ao, u, v)
    V = v + _sip_poly(hdr, "B", bo, u, v)
    # TAN deprojection (Calabretta & Greisen paper II, inverse).
    ra0r, dec0r = math.radians(ra0), math.radians(dec0)
    xr, yr = math.radians(U), math.radians(V)
    rho = math.hypot(xr, yr)
    if rho == 0.0:
        return ra0, dec0
    c = math.atan(rho)
    dec = math.asin(math.cos(c) * math.sin(dec0r)
                     + yr * math.sin(c) * math.cos(dec0r) / rho)
    ra = ra0r + math.atan2(xr * math.sin(c),
                           rho * math.cos(dec0r) * math.cos(c)
                           - yr * math.sin(dec0r) * math.sin(c))
    return math.degrees(ra) % 360.0, math.degrees(dec)


def moments(img, x, y, r=4):
    y0, y1 = max(0, y-r), min(img.shape[0], y+r+1)
    x0, x1 = max(0, x-r), min(img.shape[1], x+r+1)
    cut = img[y0:y1, x0:x1].astype(float)
    cut = np.where(np.isfinite(cut), cut, np.nanmin(cut) if np.isfinite(cut).any() else 0.0)
    cut -= cut.min()
    s = cut.sum() or 1.0
    gy, gx = np.mgrid[y0:y1, x0:x1]
    cx, cy = (gx*cut).sum()/s, (gy*cut).sum()/s
    muxx = ((gx-cx)**2*cut).sum()/s
    myyy = ((gy-cy)**2*cut).sum()/s
    area = float((cut > cut.max()*0.5).sum())
    rnd = abs(muxx-myyy)/(muxx+myyy+1e-9)
    return cx, cy, rnd, area, float(s)


def pix_scale_deg(hdr, npx=2040) -> float:
    import math
    ctype = (hdr.get('CTYPE1', '') + hdr.get('CTYPE2', '')).upper()
    if 'SIP' not in ctype:
        for k in ('CDELT1', 'CD1_1'):
            v = hdr.get(k)
            if v:
                try:
                    return abs(float(v.split('/')[0]))
                except ValueError:
                    pass
    return 3.485 / npx  # shortcut: linear TAN only, upgrade to SIP when residuals exceed 3 arcsec


def xy_to_radec(hdr, x, y, scale=None) -> tuple[float, float]:
    import math
    ra0 = float(hdr.get('CRVAL1', '0').split('/')[0])
    dec0 = float(hdr.get('CRVAL2', '0').split('/')[0])
    s = scale or pix_scale_deg(hdr)
    # linear TAN approx around field center; cos(dec) correction on RA
    return ra0 + (x - 1020.0) * s / max(math.cos(math.radians(dec0)), 0.2), dec0 + (y - 1020.0) * s


def extract_frame(path: str, epoch_id: str, mjd: float, image_id: str, nsigma=5.0, top_k=400, edge=5,
                band: str | None = None):
    from .detections import Detection
    hdr, img, var = read_image_and_variance(path)
    med, sigma = background_mad(img)
    peaks = find_peaks(img, nsigma=nsigma, edge=edge)
    peaks.sort(key=lambda p: p[2], reverse=True)
    peaks = peaks[:top_k]
    scale = pix_scale_deg(hdr)
    band = band if band else band_from_path(path)
    use_sip = bool(hdr) and "CRPIX1" in hdr
    dets = []
    for i, (x, y, peak) in enumerate(peaks):
        cx, cy, rnd, area, _ = moments(img, x, y)
        if use_sip:
            ra, dec = sip_to_radec(hdr, cx, cy)
        else:
            ra, dec = xy_to_radec(hdr, cx, cy, scale)
        vv = var[int(y), int(x)] if var is not None else float("nan")
        snr = local_snr(peak, med, vv, sigma)
        flux = aperture_flux(img, cx, cy)
        dets.append(Detection(f'FITS-{i}', epoch_id, mjd, ra, dec, snr, flux, rnd, area, cx, cy, image_id,
                              band=band, flux_band=flux if band else None))
    return dets
