**`alignment.fasta`**: one record per structure, all of the same length. The header is the
file name, the chain and US-align's statistics, tab-separated; the sequence uses `-` for gaps.

```text
>AF-B6EQB6.pdb:A	L=467	d0=7.72	seqID=0.046	TM-score=0.25379
-HGYVSE--TND-GIAGS-RA------ALCK--Y--PTSD--------T------...
```

The *n*-th letter of a row is the *n*-th residue of `aligned/<id>.pdb`, in file order. This is
how every later step maps residues to [alignment columns](../conventions.md#numbering).

**`matrices.json`**: list with one item per structure, in input order.

| Field | Meaning |
| --- | --- |
| `label` | `<family>_<organism>_<protein name>`, cut to 60 characters |
| `source_file` | `<struct_id>.pdb` |
| `usalign_key` | US-align's name for the chain, `<file>:<chain>` |
| `matrix_4x4` | 16 numbers, row-major: the rotation and translation that place the prepared structure in the common frame. The reference structure has the identity |

**`aligned/<id>.pdb`**: the prepared structure in the common coordinate frame, otherwise
unchanged (same residues, numbering, B-factor column and metal).

Also kept: `matrices.txt` and `usalign_stdout.txt` (US-align's raw output),
`structure_list.txt` (the input order) and `result.<n>.pdb` (US-align's own output files,
copied to `aligned/`).
