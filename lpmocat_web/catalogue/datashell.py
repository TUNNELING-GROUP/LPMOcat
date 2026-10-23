"""Small helpers the views use. BigDataHandler, the Django-side calculations that read the old
85-structure data in raw_data_new/, was removed on 2026-09-25: the pipeline computes everything
(lpmocat_pipeline/pipeline.py) and populate_db only loads it."""


def calculate_flexibility_stats(rmsf_text):
    """
    Calculate min and max RMSF values from RMSF text data.
    Returns (flex_min, flex_max) tuple or (None, None) if parsing fails.
    """
    try:
        rmsf_values = []
        for line in rmsf_text.split('\n'):
            line = line.strip()
            if line and not line.startswith('#'):
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        rmsf_val = float(parts[1])
                        rmsf_values.append(rmsf_val)
                    except (ValueError, IndexError):
                        continue
        
        if rmsf_values:
            return (min(rmsf_values), max(rmsf_values))
        return (None, None)
    except Exception:
        return (None, None)

def calculate_sequence_length(sequence):
    """
    Calculate sequence length excluding gaps, whitespace, and non-amino acid characters.
    Returns integer length or 0 if calculation fails.
    """
    try:
        if not sequence:
            return 0
        # Remove gaps (-), whitespace, and newlines
        sequence_clean = sequence.replace('-', '').replace(' ', '').replace('\n', '').replace('\r', '').replace('\t', '')
        # Count only valid amino acid letters (A-Z, excluding gaps and special chars)
        sequence_clean = ''.join([c for c in sequence_clean if c.isalpha()])
        return len(sequence_clean)
    except Exception:
        return 0
