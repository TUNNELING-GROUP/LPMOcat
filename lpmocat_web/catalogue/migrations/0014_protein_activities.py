from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('catalogue', '0013_protein_second_sphere_motif'),
    ]

    operations = [
        migrations.AddField(
            model_name='protein',
            name='activities',
            field=models.TextField(blank=True, help_text='Substrates, oxidised positions and the papers reporting them (pipeline step 2j: CAZy EC numbers, Swiss-Prot, BRENDA), JSON: {"ec": {ec: {source: [ref]}}, "substrates": {substrate: {"result", "forms", "active", "inactive", "oxidation", "sources", "check"}}, "observations": [...], "notes": [...], "references": {ref: {"author", "year", "title", "pmid", "doi", "journal"}}}', null=True),
        ),
    ]
