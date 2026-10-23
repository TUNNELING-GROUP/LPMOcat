from django.contrib import admin
from .models import Protein, Clade, Substrate, Organism, Protein_files, RegioSelecitivity, Electrostatics, SurfaceAccesability, Flexibility

admin.site.register(
    [
        Protein,
        Clade,
        Substrate,
        Organism,
        Protein_files,
        RegioSelecitivity,
        Electrostatics,
        SurfaceAccesability,
        Flexibility
    ]
)

# Register your models here.
