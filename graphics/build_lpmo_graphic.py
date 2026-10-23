"""Build the About page's LPMO graphic from the PowerPoint export next to this script.

    python graphics/build_lpmo_graphic.py

LPMO_structural_features.svg (PowerPoint: select the graphic, Save as Picture, SVG) holds four
PyMOL renders and their frames, arrows and labels in one flat group. This script writes
  lpmocat_web/catalogue/static/images/lpmo_graphic/*.webp   the four renders
  lpmocat_web/catalogue/templates/catalogue/_lpmo_graphic.html   the inline SVG, one group per
      feature (data-feat = site | surface | hh), which home.html and info_pages.css link to the
      coloured terms of the "About LPMOs" text
The elements are recognised by their order and colours, and the label anchors are this layout's:
after a changed layout, check the asserts and the anchor() calls below."""
import base64, io, os, re
import xml.etree.ElementTree as ET
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
SRC = f'{HERE}/LPMO_structural_features.svg'
IMG_DIR = f'{REPO}/lpmocat_web/catalogue/static/images/lpmo_graphic'
OUT = f'{REPO}/lpmocat_web/catalogue/templates/catalogue/_lpmo_graphic.html'
SVG, XLINK = 'http://www.w3.org/2000/svg', 'http://www.w3.org/1999/xlink'
ET.register_namespace('', SVG)
NAMES = {'img0': 'structure', 'img2': 'binding_surface', 'img4': 'active_site', 'img6': 'hole_hopping'}

root = ET.parse(SRC).getroot()
q = lambda tag: f'{{{SVG}}}{tag}'
defs = root.find(q('defs'))
os.makedirs(IMG_DIR, exist_ok=True)

# images -> WebP files, referenced through {% static %}
for image in defs.findall(q('image')):
    name = NAMES[image.get('id')]
    data = base64.b64decode(image.get(f'{{{XLINK}}}href').split('base64,')[1])
    Image.open(io.BytesIO(data)).save(f'{IMG_DIR}/{name}.webp', 'WEBP', quality=90, method=6)
    del image.attrib[f'{{{XLINK}}}href']
    image.set('href', f"@@static:images/lpmo_graphic/{name}.webp@@")

# ids prefixed, so they cannot clash with the page's
def prefix(value):
    return re.sub(r'#(img|clip)(\d)', r'#lg-\1\2', value)
for el in root.iter():
    if el.get('id'):
        el.set('id', 'lg-' + el.get('id'))
    for key in list(el.attrib):
        if key == f'{{{XLINK}}}href':
            el.set('href', prefix(el.attrib.pop(key)))
        elif key == 'clip-path':
            el.set(key, prefix(el.get(key)))
    if el.tag == q('text'):
        el.attrib.pop('font-family', None)     # inline SVG: the labels take the page's own font

main = [c for c in root if c.tag == q('g')][0]
items = list(main)
assert len(items) == 21, len(items)
(band, structure, surf_img, surf_frame, surf_arrow, title, site_img, site_frame, site_arrow,
 site_text, surf_text, hh_img, hh_frame, hh_arrow, hh_t1, hh_t2, hh_t3, rot1, rot2, rot3, rot_text) = items
assert band.get('fill') == '#CDE1F5' and site_frame.get('stroke') == '#8F178F' and hh_frame.get('stroke') == '#C00000'
hh_t1.text = 'hole-hopping routes'          # PowerPoint split it into three positioned pieces

def anchor(text, how, x):
    """Anchor a label at x (the group's coordinates), keeping its baseline."""
    y = text.get('transform').rstrip(')').split()[-1]
    text.set('transform', f'matrix(1 0 0 1 {x} {y})')
    text.set('text-anchor', how)

# labels centred under their frames; the right-most one ends at the frame's right edge, so no
# font can push it out of the picture
anchor(site_text, 'middle', 4170)      # frame 3655-4685
anchor(surf_text, 'middle', 4137)      # frame 3558-4716
anchor(hh_t1, 'end', 7330)             # the picture ends at 7367

for child in items:
    main.remove(child)

def group(feat, label, *children):
    g = ET.SubElement(main, q('g'), {'class': 'lg-feat', 'data-feat': feat, 'tabindex': '0',
                                     'role': 'img', 'aria-label': label})
    for child in children:
        g.append(child)
    return g

band.set('class', 'lg-band')
main.append(band)                 # the binding face on the main structure, under the render
main.append(structure)
# the title "LPMO structural features" is left out (user, 2026-10-01): the card says it
group('surface', 'Substrate-binding surface: the face around the copper, seen from above',
      surf_img, surf_frame, surf_arrow, surf_text, rot1, rot2, rot3, rot_text)
group('site', 'Active site: the histidine brace, the copper and the second-sphere residues',
      site_img, site_frame, site_arrow, site_text)
group('hh', 'Hole-hopping routes: aromatic residues between the copper and the surface',
      hh_img, hh_frame, hh_arrow, hh_t1)
for frame in (surf_frame, site_frame, hh_frame):
    frame.set('class', 'lg-frame')

for key in ('width', 'height', 'overflow'):
    root.attrib.pop(key, None)
root.set('viewBox', '0 0 3847 2809')
root.set('class', 'lpmo-graphic')
root.set('role', 'img')
root.set('aria-label', 'LPMO structural features: an LPMO structure with its active site, '
                       'substrate-binding surface and hole-hopping routes shown in three insets')

svg = ET.tostring(root, encoding='unicode')
svg = re.sub(r'@@static:([^@]+)@@', r"{% static '\1' %}", svg)
svg = svg.replace('><', '>\n<')
header = ("{% load static %}\n{% comment %}\nThe About page's LPMO graphic, built from the PowerPoint export "
          "graphics/LPMO_structural_features.svg by\ngraphics/build_lpmo_graphic.py (the four renders are static/images/lpmo_graphic/*.webp). One group per "
          "feature (data-feat = site |\nsurface | hh), which info_pages.css and home.html link to the coloured "
          "terms of the \"About LPMOs\" text.\n{% endcomment %}\n")
open(OUT, 'w', encoding='utf-8').write(header + svg + '\n')
print('partial', os.path.getsize(OUT), 'bytes;', {f: os.path.getsize(f'{IMG_DIR}/{f}') // 1024 for f in sorted(os.listdir(IMG_DIR))})
