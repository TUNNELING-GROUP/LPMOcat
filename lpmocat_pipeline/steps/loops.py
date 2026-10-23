"""
NO LONGER USED BY THE PIPELINE (2026-09-23): step 7c now fits a plane to the copper face of
each structure and takes the residues it touches (binding_surface.py). This module is kept as
the curated reference the plane method is checked against -- 80-95 % of the plane residues
fall inside these template regions, which are 3-4 times larger because they are whole loop
stretches rather than the residues that can touch a flat substrate.

Former step 7c: binding surface of every aligned structure, mapped from one
template structure per family through the structural alignment. (The loop regions themselves
come from the beta-sandwich core, see beta_core.py.)

Template ranges are ungapped residue positions of the template (1 = His-1, the first residue
of the aligned structure). Each region's first and last template positions are converted to
alignment columns, and every structure of the family gets the residues at, or nearest to,
those columns. Regions are capped at the structure's catalytic domain end (step 6). Output
positions are ungapped sequence positions, as stored in Protein.sequence and used by the
compare-page plots.

Templates: AA9/AA10/AA14/AA17 the crystal structures of the earlier catalogue; AA11 7P3U
(AfAA11; 4MAI has 29 unresolved residues in LS and L8/LC); AA13 4OPB (AoAA13); AA15 5MSZ
(TdAA15A); AA16 7ZE9 (TtAA16 = MtAA16).
"""

BINDING_SURFACE_TEMPLATES = {
    "AA9": {
        "template": "7OVA",
        "regions": {"L2a": "19-27", "L2": "37-50", "L3": "80-87",
                    "LS": "118-139", "L8": "163-172", "LC": "208-217"},
    },
    "AA10": {
        "template": "4YN1",
        "regions": {"L2a": "17-28", "L2": "49-79", "L2b": "95-108",
                    "L3": "142-146", "LS": "177-194", "L8": "218-222"},
    },
    "AA11": {
        # 7P3U, not 4MAI (unresolved residues) -- see the module docstring.
        "template": "7P3U",
        "regions": {"L2": "13-43", "L3": "67-73", "LS": "96-102",
                    "L8": "125-132", "LC": "162-172"},
    },
    "AA13": {
        # 4OPB, not 5LSV: step 2a picks 4OPB since its brace holds the copper and 5LSV's a
        # zinc. Both are AoAA13 and their aligned chains are identical (233 residues, no
        # mismatch, no gap), so the ranges below did not change with the template.
        "template": "4OPB",
        "regions": {"L2": "18-34", "L2b": "50-56", "L3": "87-92",
                    "LS": "171-191", "L8": "218-221"},
    },
    "AA14": {
        "template": "5NO7",
        "regions": {"L2": "13-33", "L3": "78-102", "L8": "167-174",
                    "LC": "229-242"},
    },
    "AA15": {
        "template": "5MSZ",
        "regions": {"L2": "21-40", "L3": "89-93", "LS": "130-137",
                    "L8": "163-181"},
    },
    "AA16": {
        "template": "7ZE9",
        "regions": {"L2": "22-45", "L3": "82-86", "LS": "111-119",
                    "L8": "152-155"},
    },
    "AA17": {
        "template": "6Z5Y",
        "regions": {"L2": "14-24", "L2b": "30-36", "L3": "92-96",
                    "LS": "122-126", "L8": "153-160"},
    },
}


def _ungapped_to_msa_col(seq, ungapped_pos):
    """Return 1-based MSA column of the Nth non-gap character."""
    count = 0
    for i, c in enumerate(seq):
        if c != "-":
            count += 1
            if count == ungapped_pos:
                return i + 1
    return None


def _msa_col_to_ungapped(seq, msa_col):
    """Return ungapped position at 1-based MSA column; 0 if gap."""
    if msa_col > len(seq):
        return 0
    char = seq[msa_col - 1]
    if char == "-":
        return 0
    return sum(1 for c in seq[:msa_col] if c != "-")


def _nearest_ungapped(seq, msa_col, window=50):
    """Ungapped position of nearest non-gap within ±window of msa_col."""
    for delta in range(0, window + 1):
        for sign in (-1, 1):
            col = msa_col + sign * delta
            if 1 <= col <= len(seq) and seq[col - 1] != "-":
                return sum(1 for c in seq[:col] if c != "-")
    return None


def map_regions(aligned, template_id, template_regions, cat_end_positions=None):
    """Map template residue ranges through the alignment to every aligned sequence.

    aligned: {struct_id: aligned sequence}; cat_end_positions: {struct_id: catalytic domain
    end as an ungapped position}. Returns {struct_id: {region: ["start-end"] or []}}, or {}
    when the template is not in the alignment."""
    template_seq = aligned.get(template_id)
    if template_seq is None:
        return {}

    # Step 1: map template ungapped positions → MSA columns
    template_last_col = _ungapped_to_msa_col(template_seq, len(template_seq.replace("-", "")))
    msa_template = {}
    for region, range_str in template_regions.items():
        start, end = map(int, range_str.split("-"))
        msa_start = _ungapped_to_msa_col(template_seq, start)
        msa_end = _ungapped_to_msa_col(template_seq, end)
        if msa_start is None:
            continue
        if msa_end is None:
            # end position is beyond the template sequence — clamp to last residue
            msa_end = template_last_col
        msa_template[region] = (msa_start, msa_end)

    # Step 2: for each sequence, map MSA columns → ungapped positions
    result = {}
    for struct_id, seq in aligned.items():
        cat_end = (cat_end_positions or {}).get(struct_id)
        record_selection = {}

        for region, (msa_start, msa_end) in msa_template.items():
            # the entire alignment segment is gaps for this sequence
            segment = seq[msa_start - 1:msa_end]
            if not segment.replace("-", ""):
                record_selection[region] = []
                continue

            ug_start = _msa_col_to_ungapped(seq, msa_start)
            if ug_start == 0:
                ug_start = _nearest_ungapped(seq, msa_start)

            ug_end = _msa_col_to_ungapped(seq, msa_end)
            if ug_end == 0:
                ug_end = _nearest_ungapped(seq, msa_end)

            # Keep every region within the catalytic domain
            if cat_end is not None:
                if ug_start and ug_start > cat_end:
                    record_selection[region] = []      # region lies entirely past the domain
                    continue
                if ug_end:
                    ug_end = min(ug_end, cat_end)

            if ug_start and ug_end and ug_start <= ug_end:
                record_selection[region] = [f"{ug_start}-{ug_end}"]
            else:
                record_selection[region] = []

        result[struct_id] = record_selection

    return result


def assign_binding_surface(aligned, families, cat_end_positions):
    """Binding surface of every structure from its family's template.
    Returns ({struct_id: {region: ["start-end"] or []} or None},
    [families whose template is missing from the alignment])."""
    result = {struct_id: None for struct_id in aligned}
    missing = []
    for family, template in sorted(BINDING_SURFACE_TEMPLATES.items()):
        members = [s for s in aligned if families.get(s) == family]
        if not members:
            continue
        mapping = map_regions(aligned, template["template"], template["regions"],
                              cat_end_positions=cat_end_positions)
        if not mapping:
            missing.append(f"{family} binding surface template {template['template']}")
            continue
        for struct_id in members:
            result[struct_id] = mapping[struct_id]
    return result, missing
