from django.urls import path, include
from . import views
from .views import activeUsersAPI

urlpatterns = [
    # Main pages
    path("", views.index, name="index"),  # /catalogue/ - Browse Catalogue
    
    # Result/Comparison page
    path("compare/proteins/<str:protein_ids>/", views.compare_proteins, name="compare_proteins"),
    # Lazy per-section plot data (fetched on demand from the plot-bundle cache)
    path("compare/plots/", views.compare_plots, name="compare_plots"),
    path("compare/coloring/", views.compare_coloring, name="compare_coloring"),
    # sequence logo and similarity heatmap of the compare page: descriptive similarity / conservation for the selected set
    path("compare/analysis/", views.compare_analysis, name="compare_analysis"),
    # Selection overview: the catalogue overview charts and the tree for the selected set
    path("compare/overview/", views.compare_overview, name="compare_overview"),

    # Downloads
    path('download/pdb-structures/', views.download_pdb_structures, name="download_pdb_structures"),
    path('download/alignment/', views.download_alignment, name="download_alignment"),
    path('download/comparison/', views.download_comparison, name="download_comparison"),
    path('downloads/', views.downloads_page, name="downloads"),
    path('downloads/<str:name>', views.download_file, name="download_file"),
    
    # API endpoints (needed for functionality)
    path('api/pdb-structures/', views.get_pdb_structures, name="get_pdb_structures"),
    path('sequence_search/', views.sequence_search, name='sequence_search'),
    path("active-users/", activeUsersAPI, name="active_users_api"),
]
