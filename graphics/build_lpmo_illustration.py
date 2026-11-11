"""Build the About page's LPMO illustration images from the layered PowerPoint export.

The export (about_page.svg + layers/*.webp, transparent backgrounds) is the package in
about-illustration/about-illustration/ (not in git). Its layout is copied by hand into
templates/catalogue/_lpmo_graphic.html; this script makes the images that layout shows, in
static/images/lpmo_illustration/:

    protein-plain, -ligands, -metals,      the protein at rest, and with the active-site sticks,
    -surface-outline                       the hole-hopping donors or the surface outline marked
    histidine-closeup                      the active-site panel
    hopping-base                           the hole-hopping panel: the donors (green)
    hopping-arrows-blue, -orange           its arrows only (transparent elsewhere): the arrow from
                                           the copper with the blue path, and the orange path
    binding-surface, binding-ribbon        the binding-surface panel: the surface, and its loops

The images are scaled to SCALE: the page shows the illustration about 500 px wide, so half the
export's resolution is still sharp on a high-density screen. The <image> elements keep their
width and height, so the layout does not change.

Usage: python graphics/build_lpmo_illustration.py [layers folder]
"""
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageFilter
from scipy import ndimage

ROOT = Path(__file__).resolve().parent.parent
SOURCE = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / 'about-illustration' / 'about-illustration' / 'layers'
OUT = ROOT / 'lpmocat_web' / 'catalogue' / 'static' / 'images' / 'lpmo_illustration'
SCALE = 0.5
QUALITY = 82
LAYERS = ['protein-plain', 'protein-ligands', 'protein-metals', 'protein-surface-outline',
          'histidine-closeup', 'hopping-base', 'binding-surface', 'binding-ribbon']
# The hole-hopping arrows: in the export's hopping-arrows layer (the copper on the left), the arrow
# from the copper ends before x = STRIPED_END; it goes with the blue path, the orange path is right of it
STRIPED_END = 310
COPPER, COPPER_RADIUS = (74, 168), 42   # the copper sphere, in both layers: never an arrow
SPECK = 80                             # arrow pieces smaller than this (px) are rendering noise


def save(name, image):
    width, height = image.size
    image = image.resize((round(width * SCALE), round(height * SCALE)), Image.LANCZOS)
    target = OUT / f'{name}.webp'
    image.save(target, 'WEBP', quality=QUALITY, alpha_quality=90, method=6)
    print(f'{target.relative_to(ROOT)}  {target.stat().st_size // 1024} KB')


def hopping_arrows():
    """The arrows of hopping-arrows (the strongly blue or orange pixels where it differs from
    hopping-base, whose donors are green) as two transparent layers: blue path, orange path."""
    base = Image.open(SOURCE / 'hopping-base.webp').convert('RGBA')
    arrows = Image.open(SOURCE / 'hopping-arrows.webp').convert('RGBA')
    hsv = np.asarray(arrows.convert('RGB').convert('HSV')).astype(float)
    hue, saturation = hsv[..., 0] * 360 / 255, hsv[..., 1] / 255
    differ = np.abs(np.asarray(base).astype(int) - np.asarray(arrows).astype(int)).sum(-1) > 40
    blue = (hue > 215) & (hue < 285)
    orange = (hue < 35) | (hue > 340)
    arrow = (blue | orange) & (saturation > 0.35) & differ & (np.asarray(arrows)[..., 3] > 0)
    rows, cols = np.indices(arrow.shape)
    arrow &= (cols - COPPER[0]) ** 2 + (rows - COPPER[1]) ** 2 > COPPER_RADIUS ** 2
    pieces, count = ndimage.label(arrow)
    sizes = ndimage.sum(arrow, pieces, range(1, count + 1))
    arrow &= np.isin(pieces, 1 + np.flatnonzero(sizes >= SPECK))
    left = np.zeros_like(arrow)
    left[:, :STRIPED_END] = True
    groups = {'blue': arrow & (blue | left), 'orange': arrow & orange & ~left}
    for name, mask in groups.items():
        alpha = Image.fromarray((mask * 255).astype('uint8'))
        alpha = alpha.filter(ImageFilter.MaxFilter(5)).filter(ImageFilter.GaussianBlur(1.2))   # the arrows' edges
        alpha = Image.fromarray(np.minimum(np.asarray(alpha), np.asarray(arrows)[..., 3]))
        layer = arrows.copy()
        layer.putalpha(alpha)
        yield f'hopping-arrows-{name}', layer


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for name in LAYERS:
        save(name, Image.open(SOURCE / f'{name}.webp').convert('RGBA'))
    for name, image in hopping_arrows():
        save(name, image)


if __name__ == '__main__':
    main()
