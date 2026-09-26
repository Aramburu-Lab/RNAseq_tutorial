"""Read one nf-core/fetchngs outdir, in the layout of the pinned release as its code writes it.

See D-006 for the layout facts and their sources.
"""

# The one definition. An upgrade is one deliberate change: this, a new fixture, and the
# versions file, which 1.13.0 renames.
FETCHNGS_VERSION = "1.12.0"
