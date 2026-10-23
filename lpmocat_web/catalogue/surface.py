"""Per-residue solvent accessibility as SurfaceAccesability.sasa_file stores it (pipeline step 5a,
loaded by populate_db): residue_num,chain,residue_name,sasa,rsa -- the SASA in A^2 and the
relative SASA in % (SASA / the residue type's maximum, Tien et al. 2013)."""
import math

from Bio.Data.PDBData import residue_sasa_scales

# Relative SASA (%) from which a residue counts as surface exposed (the compare page's RSA plot
# and 3D colouring, the overview's surface charge and hydrophobicity)
SURFACE_RSA = 25

# The maxima of pipeline step 5a, for files loaded before it stored the relative SASA
_MAX_ASA = residue_sasa_scales['Wilke']


def read_sasa_file(text):
    """[(residue number, residue name, SASA, relative SASA)] in file order; the relative SASA is
    NaN where unknown. Files without an rsa column (loaded before step 5a stored it) get it from
    the SASA and the Tien et al. maxima."""
    lines = [line for line in (text or '').splitlines() if line.strip() and not line.startswith('#')]
    if not lines:
        return []
    header = [h.strip() for h in lines[0].split(',')]
    rsa_col = header.index('rsa') if 'rsa' in header else None
    rows = []
    for line in lines[1:]:
        parts = [p.strip() for p in line.split(',')]
        try:
            number, name, sasa = int(parts[0]), parts[2].upper(), float(parts[3])
        except (ValueError, IndexError):
            continue
        if rsa_col is not None:
            try:
                rsa = float(parts[rsa_col])
            except (ValueError, IndexError):
                rsa = math.nan
        else:
            rsa = 100 * sasa / _MAX_ASA[name] if name in _MAX_ASA else math.nan
        rows.append((number, name, sasa, rsa))
    return rows
