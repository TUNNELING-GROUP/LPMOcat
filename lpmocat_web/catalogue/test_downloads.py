from django.test import SimpleTestCase, TestCase

from .downloads import _references_text, residue_csv


class ResidueCsvTests(SimpleTestCase):
    def test_merges_the_plot_exports_by_protein_and_position(self):
        export = {
            'Secondary Structure and Loops': {'content': (
                'protein,residue_number,position,residue,secondary_structure,core_strand,'
                'loop_region,on_binding_surface,disulfide_partner\n'
                'A,30,1,H,coil,,,yes,\n'
                'A,31,2,G,strand,b1,,no,Cys40\n')},
            'Charge': {'content': 'protein,residue_number,position,residue,charge\n'
                                  'A,30,1,H,1.0000\nA,31,2,G,0.0000\n'},
            'Hydrophobicity': {'content': 'protein,residue_number,position,residue,hydrophobicity\n'
                                          'A,30,1,H,-3.2000\nA,31,2,G,-0.4000\n'},
            'Flexibility': {'content': 'protein,residue_number,position,residue,cabsflex,gnm\n'
                                       'A,30,1,H,0.9,\nA,31,2,G,0.5,0.01\n'},
        }
        lines = residue_csv(export, 'Kyte-Doolittle').splitlines()
        self.assertEqual(lines[0], 'protein,residue_number,position,residue,secondary_structure,'
                                   'core_strand,loop_region,on_binding_surface,disulfide_partner,'
                                   'charge,hydrophobicity_kyte_doolittle,flexibility_cabsflex,'
                                   'flexibility_gnm')
        self.assertEqual(lines[2], 'A,31,2,G,strand,b1,,no,Cys40,0.0000,-0.4000,0.5,0.01')
        self.assertEqual(len(lines), 3)


class ReferencesTextTests(SimpleTestCase):
    def test_label_and_link_per_paper_without_repeats(self):
        cell = ('<span class="ref-entry"><a href="https://doi.org/10.1/x" title="T">Smith et al., '
                '2020</a></span><span class="ref-entry"><a href="https://doi.org/10.1/x">Smith '
                'et al., 2020</a></span><a href="https://pubmed.ncbi.nlm.nih.gov/1/">Lee &amp; Ko, '
                '2019</a>')
        self.assertEqual(_references_text(cell),
                         'Smith et al., 2020 <https://doi.org/10.1/x>; '
                         'Lee & Ko, 2019 <https://pubmed.ncbi.nlm.nih.gov/1/>')


class DownloadViewTests(TestCase):
    def test_only_listed_files_are_served(self):
        self.assertEqual(self.client.get('/catalogue/downloads/settings.py').status_code, 404)
        self.assertEqual(self.client.get('/catalogue/downloads/..%2Fdb.sqlite3').status_code, 404)

    def test_comparison_needs_proteins(self):
        self.assertEqual(self.client.get('/catalogue/download/comparison/?ids=').status_code, 400)

    def test_downloads_are_a_section_of_the_about_page(self):
        self.assertRedirects(self.client.get('/catalogue/downloads/'), '/#downloads',
                             fetch_redirect_response=False)
        response = self.client.get('/')
        self.assertContains(response, 'id="downloads"')
        self.assertContains(response, 'CC BY 4.0')
