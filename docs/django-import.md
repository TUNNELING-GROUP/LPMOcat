# Loading into the catalogue

The web app does not calculate anything about the proteins (hydrophobicity is looked up per
residue type, see [Hydrophobicity scales](reference/hydrophobicity-scales.md)). After a pipeline run, a chain of
Django management commands loads the files into the database and prepares what the pages
need. Run them from `lpmocat_web/`, in this order:

```bash
python manage.py populate_db --force
python manage.py update_protein_codes
python manage.py rename_lpmo_to_aa --apply
python manage.py clean_other_names --apply
python manage.py backfill_references
python manage.py backfill_cazy_references
python manage.py map_domain_ranges
python manage.py pregenerate_plots --force
python manage.py build_downloads
```

Back up `db.sqlite3` first; `--force` replaces existing proteins.

## The commands

| Command | What it does | Options |
| --- | --- | --- |
| `populate_db` | Creates or replaces one database record per structure of `pdb_selection.csv` from the pipeline files. With `--force` it also removes proteins whose structure is no longer selected | `--force`, `--proteins ID …`, `--dry-run`, `--data-dir DIR` |
| `update_protein_codes` | Gives each protein its short code (`SmAA10A`): a curated code from `curation/protein_codes.csv`, else the LPMO-style name found in the protein name, else one built from organism and family | `--dry-run`, `--proteins` |
| `rename_lpmo_to_aa` | Turns `…LPMO9A` codes into `…AA9A` and keeps the old name as an alias. Without `--apply` it only previews | `--apply` |
| `clean_other_names` | Removes activity notes and the protein's own code from its other names. Without `--apply` it only previews | `--apply` |
| `backfill_references` | Fetches the primary citation of every PDB entry of a protein from RCSB | `--force`, `--proteins`, `--dry-run`, `--limit`, `--delay` |
| `backfill_cazy_references` | Resolves the PubMed links of `cazy_lpmos.csv` and `other_lpmos.csv` to citations (NCBI E-utilities) | `--force`, `--proteins`, `--dry-run`, `--csv`, `--delay` |
| `map_domain_ranges` | Places the catalytic domain and CBMs on each structure's own residue numbers, for the domain bar | |
| `pregenerate_plots` | Renders the per-protein plots of the compare page once, so pages load quickly | `--force`, `--protein CODE`, `--type NAME` |
| `build_downloads` | Writes the files of the Downloads page (tables, per-residue data, alignments, structures, one archive of everything) and `release.json`, the data release every page names, into `data/django_outputs/downloads/` | `--release NAME` (default: today's date) |

The order matters: the codes must be final before references and plots are attached to them,
and the plots embed the code.

## What `populate_db` reads

Every file below is described under [Output formats](reference/outputs.md).

| Step | Files |
| --- | --- |
| 1 | `cazy_lpmos.csv`, `other_lpmos.csv` |
| 2 | `pdb_selection.csv`, `uniprot_ids.json`, `uniprot_references.json`, `domain_annotations.json`, `taxonomy.json`, `activities.json` |
| 4 | `alignment.fasta`, `aligned/<id>.pdb` |
| 5 | `sasa.json`, `secondary_structure.json` |
| 6 | `disulfides.json`, `histidine_brace.json` |
| 7 | `beta_core.json`, `loop_regions.json`, `binding_surface.json` (the plane and the copper, for the 3D viewer's views) |
| 8 | `brace_interactions.json` |
| 10 | `flexibility_raw.json` |
| 11 | `<id>.pqr` |
| 12 | `hole_hopping.json` |

A missing file is reported and its fields stay empty.

## What is derived on loading

A few values are decided at import, because they combine sources:

- **Substrates**: those a paper reports activity on (step 2j), together with those of the
  EC numbers and the first substrate the protein name mentions. Only when none of these
  gives one, the family's usual substrate.
- **Regioselectivity**: the oxidised positions papers report (step 2j); without any, the
  positions of the EC numbers, using the IUBMB definitions (ExplorEnz):

    | EC number | Enzyme | Oxidised position |
    | --- | --- | --- |
    | 1.14.99.53 | lytic chitin monooxygenase | C1 and C4 |
    | 1.14.99.54 | lytic cellulose monooxygenase (C1-hydroxylating) | C1 |
    | 1.14.99.55 | lytic starch monooxygenase | C1 |
    | 1.14.99.56 | lytic cellulose monooxygenase (C4-dehydrogenating) | C4 |

    Several EC numbers give the union; none gives "Unknown".

- **Organisms**: one per species (the first two words of the organism name, three for
  "Genus sp. X"), with the NCBI lineage of step 2i.
- **References**: structure papers (RCSB), characterization papers (step 1 lists) and the
  papers UniProt entries cite (step 2f), merged and deduplicated by DOI and PubMed ID.
- **CBMs**: the CBM families any of the protein's UniProt entries annotate (step 2g).

## Files written by the import

`map_domain_ranges`, `pregenerate_plots` and `build_downloads` write to `data/django_outputs/`, which
has no step number because it is not a pipeline output:

| File | Content |
| --- | --- |
| `domain_ranges_auth.json` | Per structure, the catalytic domain and CBM ranges in residue numbers |
| `sifts_cache.json`, `uniprot_sequences.json` | Cached SIFTS responses and UniProt sequences |
| `pregenerated_plots/<code>/<type>.json.gz` | Pre-rendered plots |
| `downloads/` | The Downloads page's files and `release.json` (the data release: its name, date, protein counts per family, BRENDA release and code version) |

## Files the app reads directly

Besides the database, the running app reads three things from the data folder:
`06_catalytic_domain/catalytic_domain_boundaries.json`,
`09_distance_matrix/distance_matrix.csv` and `django_outputs/`. A deployment therefore needs
the database and these, not the whole data folder.
