import json
import math
from types import SimpleNamespace
from unittest import mock

from django.test import SimpleTestCase, TestCase

from .management.commands.populate_db import _sasa_csv, organism_for, taxonomy_fields
from .models import Organism
from .overview import (catalogue_overview, catalogue_summary, gravy, hydrophobicity_scales, lineage_text,
                       loop_lengths, loop_profile, mean_hydropathy, regioselectivity_label, residue_profile,
                       residue_set_profile, structural_features,
                       species_name, structure_source, taxon_group)
from .protein_mass import residues_kda
from .surface import read_sasa_file

NEUROSPORA = [['cellular root', 'cellular organisms'], ['domain', 'Eukaryota'], ['clade', 'Opisthokonta'],
              ['kingdom', 'Fungi'], ['phylum', 'Ascomycota'], ['class', 'Sordariomycetes'],
              ['species', 'Neurospora crassa']]


def _protein(code, family, struct_id, sequence, cat_length=1, cbms=None, kingdom='Eukaryota', lineage=None,
             loops=None, regio='C1-hydroxylating', activities=None):
    substrates = [SimpleNamespace(name='Cellulose')]
    regioselectivity = [SimpleNamespace(name=regio)] if regio else []
    return SimpleNamespace(
        protein_code=code, clade=SimpleNamespace(name=family),
        organism=SimpleNamespace(name='Neurospora crassa OR74A', kingdom=kingdom,
                                 lineage=json.dumps(lineage) if lineage else None),
        substrate=SimpleNamespace(all=lambda: substrates), associated_CBMs=cbms,
        regioselectivity=SimpleNamespace(all=lambda: regioselectivity),
        activities=json.dumps(activities) if activities else None,
        selected_structure=struct_id, cat_domain_length=cat_length, masked=sequence,
        sequence=sequence.replace('-', '') + 'W',
        loop_regions=str(loops) if loops else None)


# a step-2j record whose papers report C1 and C4 oxidation of chitin
REPORTED_C1_C4 = {
    'ec': {}, 'notes': [], 'references': {},
    'substrates': {'chitin': {'result': 'active', 'forms': {}, 'active': ['PMID:1'], 'inactive': [],
                              'oxidation': {'C1': ['PMID:1'], 'C4': ['PMID:1']}, 'sources': ['BRENDA'],
                              'check': False}},
    'observations': [
        {'source': 'BRENDA', 'substrate': 'chitin', 'form': 'beta-chitin', 'result': 'active',
         'oxidation': ['C1', 'C4'], 'only': False, 'refs': ['PMID:1'], 'text': ''},
        # an EC number is not a paper's statement about the substrate
        {'source': 'CAZy', 'substrate': 'cellulose', 'form': None, 'result': 'active',
         'oxidation': [], 'only': False, 'refs': ['PMID:2'], 'text': 'EC 1.14.99.54'},
    ],
}


class GravyTests(SimpleTestCase):
    def test_kyte_doolittle_mean(self):
        self.assertAlmostEqual(gravy('IR'), (4.5 - 4.5) / 2)
        self.assertAlmostEqual(gravy('AG'), (1.8 - 0.4) / 2)

    def test_gaps_case_and_unknown_letters_are_ignored(self):
        self.assertEqual(gravy('a-G\nXB'), gravy('AG'))

    def test_no_residues(self):
        self.assertIsNone(gravy('--'))
        self.assertIsNone(gravy(None))


class RegioselectivityLabelTests(SimpleTestCase):
    def test_short_labels(self):
        self.assertEqual(regioselectivity_label(['C1-hydroxylating']), 'C1')
        self.assertEqual(regioselectivity_label(['C1-hydroxylating/C4-dehydrogenating']), 'C1/C4')
        # several names: the positions of all of them
        self.assertEqual(regioselectivity_label(['C4-dehydrogenating', 'C1-hydroxylating']), 'C1/C4')
        self.assertEqual(regioselectivity_label(['Unknown']), 'Unknown')
        self.assertEqual(regioselectivity_label([]), 'Unknown')


class StructureSourceTests(SimpleTestCase):
    def test_sources(self):
        self.assertEqual(structure_source('6TWE'), 'PDB')
        self.assertEqual(structure_source('AF-Q9RFX5'), 'AlphaFold DB')
        self.assertEqual(structure_source('af-q9rfx5'), 'AlphaFold DB')
        self.assertEqual(structure_source('C3D-BDB32680_1'), 'CAZyme3D')
        self.assertEqual(structure_source('AF3-BaAA9A'), 'own model')
        self.assertEqual(structure_source(None), 'PDB')


def _pqr_line(serial, name, resname, number, charge, icode=' '):
    return (f"ATOM  {serial:5d} {name:<4} {resname:>3} A{number:4d}{icode}   "
            f"{0.0:8.3f}{0.0:8.3f}{0.0:8.3f} {charge:7.4f} 1.8240")


class ResidueProfileTests(SimpleTestCase):
    # residues 1 (exposed Lys, +1), 2 (buried Asp, -1), 2A (exposed Glu, -1), 3 (exposed Arg, +1)
    PQR = '\n'.join([
        'REMARK   6 Total charge on this protein: 0.0000 e',
        _pqr_line(1, 'N', 'LYS', 1, 0.5), _pqr_line(2, 'NZ', 'LYS', 1, 0.5),
        _pqr_line(3, 'OD1', 'ASP', 2, -1.0),
        _pqr_line(4, 'OE1', 'GLU', 2, -1.0, icode='A'),
        _pqr_line(5, 'NH1', 'ARG', 3, 1.0),
    ])
    # as pipeline step 5a stores it: SASA (A^2) and relative SASA (%)
    SASA = '\n'.join(['residue_num,chain,residue_name,sasa,rsa',
                      '1,A,LYS,120.0,50.85', '2,A,ASP,5.0,2.59', '2,A,GLU,90.0,40.36', '3,A,ARG,200.0,72.99'])

    def test_all_and_surface_residues(self):
        profile = residue_profile(self.PQR, self.SASA, catalytic_end=2)
        self.assertEqual(profile['charge'], {'all': 0.0, 'cd': -1.0})
        # the insertion-code residue 2A pairs with its own SASA row (by position): exposed
        self.assertEqual(profile['surface_charge'], {'all': 1.0, 'cd': 0.0})
        kd = hydrophobicity_scales()['Kyte-Doolittle']
        self.assertEqual(profile['surface_hydrophobicity']['Kyte-Doolittle'],
                         {'all': mean_hydropathy('KER', kd), 'cd': mean_hydropathy('KE', kd)})
        # Lys 1, Glu 2A, Arg 3 exposed; Asp 2 buried
        self.assertEqual(profile['surface_length'], {'all': 3, 'cd': 2})
        self.assertEqual(profile['surface_mass'],
                         {'all': round(residues_kda('KER'), 2), 'cd': round(residues_kda('KE'), 2)})

    def test_without_catalytic_end_or_files(self):
        profile = residue_profile(self.PQR, '', catalytic_end=None)
        self.assertEqual(profile['charge'], {'all': 0.0, 'cd': None})
        self.assertEqual(profile['surface_charge'], {'all': None, 'cd': None})
        self.assertEqual(profile['surface_length'], {'all': None, 'cd': None})
        self.assertIsNone(profile['surface_hydrophobicity'])
        profile = residue_profile('', self.SASA, 2)
        self.assertEqual(profile['charge'], {'all': None, 'cd': None})
        self.assertEqual(profile['surface_hydrophobicity']['Kyte-Doolittle']['cd'],
                         mean_hydropathy('KE', hydrophobicity_scales()['Kyte-Doolittle']))

    def test_rows_matched_by_number_when_the_files_differ(self):
        sasa = '\n'.join(['residue_num,chain,residue_name,sasa,rsa', '1,A,LYS,120.0,50.85', '3,A,ARG,1.0,0.36'])
        profile = residue_profile(self.PQR, sasa, catalytic_end=None)
        self.assertEqual(profile['surface_charge']['all'], 1.0)   # only Lys 1 is exposed


class SasaFileTests(SimpleTestCase):
    def test_relative_sasa_column(self):
        rows = read_sasa_file('residue_num,chain,residue_name,sasa,rsa\n1,A,HIS,64.4809,28.79\n2,A,THR,6.4,')
        self.assertEqual(rows[0], (1, 'HIS', 64.4809, 28.79))
        self.assertTrue(math.isnan(rows[1][3]))   # empty = unknown

    def test_files_without_relative_sasa(self):
        # loaded before step 5a stored it: from the Tien et al. maximum (His 224 A^2)
        rows = read_sasa_file('residue_num,chain,residue_name,sasa\n1,A,HIS,56.0\n2,A,UNK,5.0')
        self.assertAlmostEqual(rows[0][3], 25.0)
        self.assertTrue(math.isnan(rows[1][3]))
        self.assertEqual(read_sasa_file('# no SASA: run pipeline step 5 (dssp_sasa)'), [])

    def test_loader_format(self):
        entry = {'residues': [1, 2], 'chains': ['A', 'A'], 'names': ['HIS', 'UNK'],
                 'sasa': [64.48087, 5.0], 'rsa': [28.79, None]}
        text = _sasa_csv(entry)
        self.assertEqual(text.splitlines(), ['residue_num,chain,residue_name,sasa,rsa',
                                             '1,A,HIS,64.4809,28.79', '2,A,UNK,5.0000,'])
        self.assertEqual(read_sasa_file(text)[0], (1, 'HIS', 64.4809, 28.79))


class HydrophobicityScaleTests(SimpleTestCase):
    def test_scales_from_the_compare_page_file(self):
        scales = hydrophobicity_scales()
        self.assertEqual(list(scales)[0], 'Kyte-Doolittle')
        self.assertIn('Hopp-Woods', scales)
        self.assertEqual(scales['Kyte-Doolittle']['I'], 4.5)

    def test_mean_skips_unknown_letters(self):
        scale = {'A': 1.0, 'G': -1.0}
        self.assertEqual(mean_hydropathy('A-AG X', scale), 0.333)
        self.assertIsNone(mean_hydropathy('--', scale))


class LoopLengthTests(SimpleTestCase):
    def test_lengths_from_the_stored_ranges(self):
        self.assertEqual(loop_lengths("{'L2': ['6-80'], 'L3': ['95-99']}"),
                         {'L2': [75, '6-80'], 'L3': [5, '95-99']})

    def test_pieces_are_summed(self):
        self.assertEqual(loop_lengths("{'LC': ['1-3', '10-12']}"), {'LC': [6, '1-3, 10-12']})

    def test_unset_or_unreadable(self):
        self.assertEqual(loop_lengths(None), {})
        self.assertEqual(loop_lengths('not a dict ('), {})
        self.assertEqual(loop_lengths("{'L2': ['x-y']}"), {})


class LoopProfileTests(SimpleTestCase):
    SEQUENCE = 'HGDKLE'          # residues 1-6
    CHARGES = [0.0, 0.0, -1.0, 1.0, 0.0, -1.0]
    RSA = [10.0, 50.0, 30.0, 5.0, math.nan, 80.0]
    RMSF = [1.0, 2.0, 3.0, None, 5.0, 6.0]

    def test_values_of_the_loop_residues(self):
        # the binding surface: D3 of L2 and E6 of LC
        profile = loop_profile("{'L2': ['2-4'], 'LC': ['5-5', '6-6']}", self.SEQUENCE,
                               self.CHARGES, self.RSA, self.RMSF,
                               "{'L2': ['3-3'], 'LC': ['6-6']}")
        kd = hydrophobicity_scales()['Kyte-Doolittle']
        l2 = profile['L2']
        # residue 4 has no RMSF value; exposed: G2, D3; substrate-binding: D3
        self.assertEqual(l2['rmsf'], [2.5, 2.5, 3.0])
        self.assertEqual(l2['residues'], [3, 2, 1])
        # all; exposed: G2 (0), D3 (-1); substrate-binding: D3 (-1)
        self.assertEqual(l2['charge'], [0.0, -1.0, -1.0])
        self.assertEqual(l2['hydrophobicity'][0],
                         [mean_hydropathy('GDK', kd), mean_hydropathy('GD', kd), kd['D']])
        self.assertEqual(len(l2['hydrophobicity']), len(hydrophobicity_scales()))
        # pieces together; a NaN relative SASA is never exposed
        self.assertEqual(profile['LC']['rmsf'], [5.5, 6.0, 6.0])
        self.assertEqual(profile['LC']['charge'], [-1.0, -1.0, -1.0])
        self.assertEqual(profile['LC']['hydrophobicity'][0],
                         [mean_hydropathy('LE', kd), kd['E'], kd['E']])

    def test_without_a_binding_surface(self):
        profile = loop_profile("{'L2': ['2-4']}", self.SEQUENCE, self.CHARGES, self.RSA,
                               self.RMSF)
        self.assertEqual(profile['L2']['charge'][2], None)
        self.assertEqual(profile['L2']['hydrophobicity'][0][2], None)
        self.assertEqual(profile['L2']['residues'][2], None)
        self.assertIsNone(profile['L2']['rmsf'][2])

    def test_a_loop_the_binding_surface_does_not_touch(self):
        # the surface is in LC, so L2 has none of it: a count of 0, not "unknown"
        profile = loop_profile("{'L2': ['2-4'], 'LC': ['6-6']}", self.SEQUENCE, self.CHARGES,
                               self.RSA, self.RMSF, "{'LC': ['6-6']}")
        self.assertEqual(profile['L2']['residues'], [3, 2, 0])
        self.assertIsNone(profile['L2']['rmsf'][2])
        self.assertEqual(profile['LC']['residues'][2], 1)

    def test_missing_or_mismatched_data(self):
        profile = loop_profile("{'L2': ['2-4']}", self.SEQUENCE, None, self.RSA[:3], self.RMSF)
        l2 = profile['L2']
        self.assertEqual(l2['charge'], [None, None, None])
        self.assertIsNone(l2['hydrophobicity'][0][1])            # no surface without SASA
        self.assertIsNotNone(l2['hydrophobicity'][0][0])
        self.assertEqual(l2['rmsf'][0], 2.5)
        self.assertIsNone(l2['residues'][2])                     # no binding surface given
        # no surface-exposed residue: no surface value
        buried = loop_profile("{'L3': ['4-4']}", self.SEQUENCE, self.CHARGES, self.RSA, None)['L3']
        self.assertEqual(buried['charge'], [1.0, None, None])
        self.assertEqual(buried['rmsf'], [None, None, None])
        self.assertEqual(loop_profile(None, self.SEQUENCE, None, None, None), {})


class TaxonGroupTests(SimpleTestCase):
    def test_groups(self):
        self.assertEqual(taxon_group(NEUROSPORA), ('Fungi', 'Ascomycota'))
        insect = [['domain', 'Eukaryota'], ['kingdom', 'Metazoa'], ['phylum', 'Arthropoda'],
                  ['class', 'Insecta'], ['species', 'Thermobia domestica']]
        self.assertEqual(taxon_group(insect), ('Animals', 'Insecta'))
        snail = [['domain', 'Eukaryota'], ['kingdom', 'Metazoa'], ['phylum', 'Mollusca']]
        self.assertEqual(taxon_group(snail), ('Animals', 'Mollusca'))
        oomycete = [['domain', 'Eukaryota'], ['clade', 'Sar'], ['phylum', 'Oomycota']]
        self.assertEqual(taxon_group(oomycete), ('Oomycetes', None))
        bacterium = [['domain', 'Bacteria'], ['kingdom', 'Bacillati'], ['phylum', 'Actinomycetota']]
        self.assertEqual(taxon_group(bacterium), ('Bacteria', 'Actinomycetota'))
        virus = [['acellular root', 'Viruses'], ['kingdom', 'Bamfordvirae'], ['family', 'Poxviridae']]
        self.assertEqual(taxon_group(virus), ('Viruses', 'Poxviridae'))
        metagenome = [['no rank', 'unclassified entries'], ['no rank', 'metagenomes'],
                      ['species', 'compost metagenome']]
        self.assertEqual(taxon_group(metagenome), ('Metagenomes', None))

    def test_without_lineage_the_cazy_kingdom(self):
        self.assertEqual(taxon_group([], 'Bacteria'), ('Bacteria', None))
        self.assertEqual(taxon_group([], 'unclassified'), ('Unclassified', None))
        self.assertEqual(taxon_group([], ''), ('Not recorded', None))

    def test_lineage_text(self):
        self.assertEqual(lineage_text(NEUROSPORA), 'Eukaryota › Fungi › Ascomycota › Sordariomycetes')
        metagenome = [['no rank', 'unclassified entries'], ['no rank', 'metagenomes'],
                      ['species', 'compost metagenome']]
        self.assertEqual(lineage_text(metagenome), 'unclassified entries › metagenomes')

    def test_organism_fields(self):
        entry = {'taxid': 5141, 'name': 'Neurospora crassa', 'rank': 'species', 'matched': 'species',
                 'lineage': NEUROSPORA, 'species': {'taxid': 5141, 'name': 'Neurospora crassa'}}
        self.assertEqual(taxonomy_fields(entry), {'taxid': 5141, 'lineage': json.dumps(NEUROSPORA)})

    def test_organism_fields_of_a_strain_are_the_species(self):
        entry = {'taxid': 367110, 'name': 'Neurospora crassa OR74A', 'rank': 'strain', 'matched': 'name',
                 'lineage': NEUROSPORA + [['strain', 'Neurospora crassa OR74A']],
                 'species': {'taxid': 5141, 'name': 'Neurospora crassa'}}
        self.assertEqual(taxonomy_fields(entry), {'taxid': 5141, 'lineage': json.dumps(NEUROSPORA)})

    def test_organism_fields_above_species_level(self):
        genus = NEUROSPORA[:-1] + [['genus', 'Neurospora']]
        entry = {'taxid': 5140, 'name': 'Neurospora', 'rank': 'genus', 'matched': 'genus',
                 'lineage': genus, 'species': None}
        self.assertEqual(taxonomy_fields(entry), {'taxid': 5140, 'lineage': json.dumps(genus)})
        self.assertEqual(taxonomy_fields(None), {'taxid': None, 'lineage': None})

    def test_species_name(self):
        for name, species in [('Trichoderma reesei QM6A', 'Trichoderma reesei'),
                              ('Bacillus thuringiensis serovar kurstaki ACCC 10066', 'Bacillus thuringiensis'),
                              ('Pestalotiopsis sp. NCi6', 'Pestalotiopsis sp. NCi6'),
                              ('Lentinus similis', 'Lentinus similis'),
                              ('unidentified entomopoxvirus / Melolontha melolontha entomopoxvirus (MMEV)',
                               'unidentified entomopoxvirus'),
                              ('Neurospora', 'Neurospora'), ('', '')]:
            self.assertEqual(species_name(name), species)
            self.assertEqual(species_name(species), species)


class OrganismForTests(TestCase):
    def test_strains_share_their_species(self):
        strain = {'taxid': 367110, 'name': 'Neurospora crassa OR74A', 'rank': 'strain', 'matched': 'name',
                  'lineage': NEUROSPORA + [['strain', 'Neurospora crassa OR74A']],
                  'species': {'taxid': 5141, 'name': 'Neurospora crassa'}}
        first = organism_for('Neurospora crassa OR74A', 'Eukaryota', strain)
        second = organism_for('Neurospora crassa FGSC 2489', '', None)
        self.assertEqual(first.pk, second.pk)
        self.assertEqual(Organism.objects.count(), 1)
        organism = Organism.objects.get()
        self.assertEqual((organism.name, organism.kingdom, organism.taxid),
                         ('Neurospora crassa', 'Eukaryota', 5141))
        self.assertEqual(json.loads(organism.lineage), NEUROSPORA)
        self.assertEqual(taxonomy_fields(None), {'taxid': None, 'lineage': None})


class CatalogueOverviewTests(SimpleTestCase):
    def test_records_and_family_order(self):
        proteins = [
            _protein('SmAA10', 'AA10', 'AF-X', 'HG--A', cbms='CBM2, CBM5/12', lineage=NEUROSPORA,
                     regio='C1-hydroxylating/C4-dehydrogenating', activities=REPORTED_C1_C4),
            _protein('NcAA9C', 'AA9', '4D7U', 'HGA', cbms='nan', loops={'L2': ['20-60'], 'LX': ['70-72']},
                     regio='C4-dehydrogenating'),
        ]
        proteins[0].uniprot = 'Q9X1'                 # NcAA9C has no active UniProt entry
        proteins[1].cat_domain_mass, proteins[1].sequence_mass, proteins[1].sequence_length = 21.3, 30.1, 280
        proteins[1].loop_profile = {'L2': {'rmsf': 1.5}}
        proteins[1].residue_profile = {'charge': {'cd': -6.5, 'all': -8.0},
                                       'surface_charge': {'cd': -4.0, 'all': -5.0},
                                       'surface_hydrophobicity': {'Kyte-Doolittle': {'cd': 0.5, 'all': 0.25}}}
        proteins += [
            _protein('NoEnd', 'AA11', '4MAI', 'HGA', cat_length=None, kingdom='', regio='Unknown'),
        ]
        # the structure holds only the CBM2 (the CBM5/12 lies outside it)
        ranges = {'AF-X': [{'type': 'catalytic', 'start_residue': 1, 'end_residue': 200},
                           {'type': 'CBM', 'accession': 'IPR001919', 'start_residue': 260, 'end_residue': 360,
                            'family': 'CBM2'}]}
        with mock.patch('catalogue.cbms._domain_ranges', return_value=ranges):
            data = catalogue_overview(proteins, lambda p: p.masked)

        self.assertEqual([f['name'] for f in data['families']], ['AA9', 'AA10', 'AA11'])
        self.assertEqual(data['families'][0]['color'], '#4a3aa7')
        # stacked in the tree's family order
        self.assertEqual(data['stack_order'], ['AA9', 'AA11', 'AA10'])
        sm, nc, no_end = data['proteins']
        self.assertEqual(sm['cbms'], {'uniprot': ['CBM2', 'CBM5/12'], 'structure': ['CBM2'], 'uniprot_known': True})
        self.assertEqual(nc['cbms'], {'uniprot': [], 'structure': [], 'uniprot_known': False})
        self.assertEqual(sm['source'], 'AlphaFold DB')
        self.assertEqual(sm['species'], 'Neurospora crassa')
        self.assertEqual(sm['substrates'], ['Cellulose'])
        # regioselectivity, and whether papers report it or only the EC numbers give it
        self.assertEqual(sm['regio'], ['C1/C4', 'papers'])
        self.assertEqual(sm['substrate_papers'], ['Chitin'])
        self.assertEqual(nc['substrate_papers'], [])
        self.assertEqual(nc['regio'], ['C4', 'ec'])
        self.assertEqual(no_end['regio'], ['Unknown', 'none'])
        # GRAVY of the catalytic-domain residues (gaps skipped) and of the whole structure
        self.assertEqual(sm['measures']['hydrophobicity']['Kyte-Doolittle'], {'cd': gravy('HGA'), 'all': gravy('HGAW')})
        self.assertEqual(set(sm['measures']['hydrophobicity']), {s['name'] for s in data['hydrophobicity_scales']})
        self.assertEqual(data['hydrophobicity_scales'][0], {'name': 'Kyte-Doolittle', 'hydrophilic': False})
        # no catalytic domain end: no catalytic length and GRAVY, the whole structure still counts
        self.assertIsNone(no_end['measures']['length']['cd'])
        self.assertIsNone(no_end['measures']['hydrophobicity']['Kyte-Doolittle']['cd'])
        self.assertEqual(no_end['measures']['hydrophobicity']['Kyte-Doolittle']['all'], gravy('HGAW'))
        self.assertEqual(no_end['kingdom'], '')
        self.assertEqual(sm['taxon'], ['Fungi', 'Ascomycota'])
        self.assertEqual(sm['lineage'], 'Eukaryota › Fungi › Ascomycota › Sordariomycetes')
        self.assertEqual(no_end['taxon'], ['Not recorded', None])
        self.assertEqual(nc['loops'], {'L2': [41, '20-60'], 'LX': [3, '70-72']})
        # length, mass and charges as the browse view set them; absent -> None
        self.assertEqual(nc['measures']['length'], {'cd': 1, 'all': 280})
        self.assertEqual(nc['measures']['mass'], {'cd': 21.3, 'all': 30.1})
        self.assertEqual(nc['measures']['charge'], {'cd': -6.5, 'all': -8.0})
        self.assertEqual(nc['measures']['surface_charge'], {'cd': -4.0, 'all': -5.0})
        self.assertEqual(nc['measures']['surface_hydrophobicity']['Kyte-Doolittle'], {'cd': 0.5, 'all': 0.25})
        self.assertEqual(nc['measures']['surface_hydrophobicity']['Hopp-Woods'], {'cd': None, 'all': None})
        self.assertEqual(sm['measures']['mass'], {'cd': None, 'all': None})
        self.assertEqual(sm['measures']['charge'], {'cd': None, 'all': None})
        self.assertEqual(sm['loops'], {})
        self.assertEqual(sm['loop_values'], {})
        self.assertEqual(nc['loop_values'], {'L2': {'rmsf': 1.5}})
        # the known loops in sequence order, then the others
        self.assertEqual(data['loop_names'], ['L2', 'L3', 'LS', 'L8', 'LC', 'LX'])


class CatalogueSummaryTests(SimpleTestCase):
    def test_counts(self):
        proteins = [
            _protein('NcAA9C', 'AA9', '4D7U', 'HGA', lineage=NEUROSPORA),
            _protein('SmAA10', 'AA10', 'AF-X', 'HGA', lineage=NEUROSPORA),
            _protein('BpAA10', 'AA10', 'C3D-Y_1', 'HGA', kingdom=''),
        ]
        proteins[2].organism.name = 'Bacillus sp. X'
        summary = catalogue_summary(proteins)
        self.assertEqual(summary['proteins'], 3)
        self.assertEqual((summary['experimental'], summary['predicted']), (1, 2))
        self.assertEqual((summary['alphafold'], summary['cazyme3d']), (1, 1))
        self.assertEqual(summary['species'], 2)
        self.assertEqual(summary['families'], ['AA9', 'AA10'])
        self.assertEqual(summary['taxon_groups'], [('Fungi', 2), ('Not recorded', 1)])


class ResidueSetProfileTests(SimpleTestCase):
    """The loop regions together and the substrate-binding residues, the two residue sets the
    distribution cards offer besides all and surface-exposed residues."""
    SEQUENCE = 'HGDKLE'          # residues 1-6
    CHARGES = [0.0, 0.0, -1.0, 1.0, 0.0, -1.0]
    RSA = [10.0, 50.0, 30.0, 5.0, math.nan, 80.0]
    LOOPS = "{'L2': ['2-4'], 'LC': ['6-6']}"
    SURFACE = "{'L2': ['3-3'], 'other': ['1-1']}"     # D3 and the N-His

    def test_both_sets(self):
        profile = residue_set_profile(self.LOOPS, self.SURFACE, self.SEQUENCE, self.CHARGES,
                                      self.RSA)
        kd = hydrophobicity_scales()['Kyte-Doolittle']
        self.assertEqual(profile['loops']['residues'], 4)            # G2 D3 K4 E6
        self.assertEqual(profile['loops']['charge'], -1.0)
        self.assertEqual(profile['loops']['hydrophobicity']['Kyte-Doolittle'],
                         mean_hydropathy('GDKE', kd))
        self.assertEqual(profile['binding']['residues'], 2)          # H1 D3, both regions
        self.assertEqual(profile['binding']['charge'], -1.0)
        self.assertEqual(profile['binding']['hydrophobicity']['Kyte-Doolittle'],
                         mean_hydropathy('HD', kd))
        # summed residue masses, no terminal water: the sets are not one chain
        self.assertEqual(profile['loops']['mass'], round(residues_kda('GDKE'), 2))
        self.assertEqual(profile['binding']['mass'], round(residues_kda('HD'), 2))

    def test_empty_and_missing(self):
        profile = residue_set_profile(None, None, self.SEQUENCE, self.CHARGES, self.RSA)
        self.assertEqual(profile['loops'], {'residues': None, 'mass': None, 'charge': None,
                                            'hydrophobicity': {name: None for name
                                                               in hydrophobicity_scales()}})
        # without the PQR file there is no charge, but the sequence still gives hydrophobicity
        no_charges = residue_set_profile(self.LOOPS, self.SURFACE, self.SEQUENCE, None, self.RSA)
        self.assertIsNone(no_charges['loops']['charge'])
        self.assertIsNotNone(no_charges['loops']['hydrophobicity']['Kyte-Doolittle'])

    def test_positions_outside_the_sequence_are_clipped(self):
        profile = residue_set_profile("{'L2': ['5-99']}", "{'L2': ['0-2']}", self.SEQUENCE,
                                      self.CHARGES, self.RSA)
        self.assertEqual(profile['loops']['residues'], 2)            # L5 E6
        self.assertEqual(profile['binding']['residues'], 2)          # H1 G2


class StructuralFeaturesTests(SimpleTestCase):
    """The extra per-protein values the grouping analysis scores. Residue numbers start at 20 here,
    so a test also shows the disulfide regions are looked up by sequence position."""
    SEQUENCE = 'HCNGSKCDEFWYC'                       # positions 1-13
    NUMBERS = list(range(20, 33))                    # residue numbers 20-32
    LOOPS = "{'L2': ['2-4'], 'LC': ['12-13']}"
    CORE = '{"b1": ["5-6"], "b2": ["8-8"]}'
    BONDS = json.dumps([[21, 26, 2.05], [26, 32, 2.1], [21, 40, 2.0]])   # C2-C7, C7-C13, C2-past end
    ALIGNED = 'HC-NGSKCDEF-WYC'                      # residue 30 (W11) in column 13, 31 (Y12) in 14
    PATHS = json.dumps([{'rank': 1, 'path': ['TYR31', 'TRP30', 'CU300']},
                        {'rank': 2, 'path': ['TRP30', 'HIS20', 'NI1231']}])

    def features(self, **kw):
        args = dict(sequence=self.SEQUENCE, domain=self.SEQUENCE, residue_numbers=self.NUMBERS,
                    catalytic_end=32, disulfides=self.BONDS, loop_regions=self.LOOPS, beta_core=self.CORE,
                    binding_surface="{'L2': ['2-3'], 'other': ['9-11']}",
                    aligned_sequence=self.ALIGNED, hole_hopping_paths=self.PATHS)
        args.update(kw)
        return structural_features(**args)

    def test_disulfides_inside_the_domain_by_region(self):
        f = self.features()
        self.assertEqual(f['disulfides'], 2)                       # the bond to residue 40 is past the end
        self.assertEqual(f['disulfide_pattern'], 'L2–other + other–LC')   # C2 (L2) - C7 (other), C7 - C13 (LC)
        self.assertEqual(self.features(disulfides='[]')['disulfide_pattern'], 'none')
        self.assertIsNone(self.features(disulfides=None)['disulfides'])

    def test_core_hairpin_and_sequence_values(self):
        f = self.features()
        self.assertEqual(f['b2'], 'yes')
        self.assertEqual(self.features(beta_core='{"b1": ["5-6"]}')['b2'], 'no')
        self.assertEqual(f['glycosites'], 1)                       # N3-G-S5; not N-P-S/T
        self.assertEqual(self.features(domain='NPSNGT')['glycosites'], 1)
        self.assertIsNotNone(f['pi'])

    def test_binding_surface_make_up(self):
        f = self.features()                                        # C2 N3, E9 F10 W11
        self.assertEqual((f['binding_aromatic'], f['binding_amide'], f['binding_charged']), (2, 1, 1))
        self.assertIsNone(self.features(binding_surface=None)['binding_aromatic'])

    def test_binding_aromatics_include_his_but_not_the_brace(self):
        from .overview import binding_aromatics
        # positions 1-6 = residues 20-25: H20 (brace) Y21 A22 H23 W24 F25, all on the surface
        found = binding_aromatics('HYAHWF', range(20, 26), "{'L2': ['2-3'], 'other': ['1-1', '4-6']}",
                                  "{'L2': ['2-3']}", None, brace=(20, 90))
        self.assertEqual([a['label'] for a in found], ['Tyr21', 'His23', 'Trp24', 'Phe25'])
        self.assertEqual(found[0]['region'], 'L2')
        self.assertIsNone(binding_aromatics('HYA', range(3), None))

    def test_hole_hopping_positions_are_alignment_columns(self):
        self.assertEqual(self.features()['hh_positions'], ['1:H', '13:W', '14:Y'])   # metals left out
        self.assertEqual(self.features(hole_hopping_paths='[]')['hh_positions'], [])
        self.assertIsNone(self.features(hole_hopping_paths='')['hh_positions'])
        self.assertIsNone(self.features(aligned_sequence=None)['hh_positions'])
