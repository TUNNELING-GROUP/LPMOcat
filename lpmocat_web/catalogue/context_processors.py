"""Template context shared by every page."""
from .downloads import release_label


def data_release(request):
    """data_release: 'Data release 2026-10-21' (build_downloads' release.json), '' before it ran."""
    return {'data_release': release_label()}
