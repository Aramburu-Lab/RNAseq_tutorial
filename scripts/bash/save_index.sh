# /// codebase-meta
# exempt: run by hand on Kebnekaise after an aligned nf-core/rnaseq run; copies its indexes, runs no step
# ///
# Copies the STAR and Salmon indexes one aligned nf-core/rnaseq run built (rnaseq.sbatch with
# save_reference) into a shared folder, with provenance.tsv naming what they were built from, so
# later runs on the same genome and annotation pass them in (NFCORE_INDEX) and build neither
# (Rytis 2026-10-23). The run's results stay as they are: the indexes are copied, never moved.
# A folder that exists already is refused, never overwritten.
#
# Usage:
#   bash "$CTX_REPO/scripts/bash/save_index.sh" "$RUN" "$Project_d/Rytis/Data/Reference/Human/Gencode/v50_primary_assembly/nf-core-rnaseq-3.27.0"
set -euo pipefail

RUN="${1:?usage: save_index.sh <run folder> <new index folder>}"
DEST="${2:?usage: save_index.sh <run folder> <new index folder>}"
die() {
    printf 'save_index.sh: %s\n' "$*" >&2
    exit 1
}
R="$RUN/nfcore_rnaseq/results"
[[ -e "$DEST" ]] && die "$DEST exists; an index folder is never overwritten"
[[ -s "$R/genome/index/star/SA" ]] || die "no STAR index at $R/genome/index/star; was the run saved with save_reference?"
[[ -d "$R/genome/index/salmon" ]] || die "no Salmon index at $R/genome/index/salmon"
command -v md5sum >/dev/null || die "no md5sum on PATH"
mkdir -p "$DEST.part"
python3 - "$R" "$DEST.part/provenance.tsv" <<'PY'
import glob, hashlib, json, re, sys
results, out = sys.argv[1:]
params = json.load(open(sorted(glob.glob(f"{results}/pipeline_info/params_*.json"))[-1]))
versions = open(f"{results}/pipeline_info/nf_core_rnaseq_software_mqc_versions.yml").read()

def tool(process, name):
    # The version a build step reports, only where this run ran that step.
    m = re.search(rf"^{process}:\n(?:  .*\n)*?  {name}: (\S+)", versions, re.M)
    if not m:
        sys.exit(f"save_index.sh: this run did not build its own {process}; nothing saved")
    return m.group(1)

def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()

release = re.search(r"nf-core/rnaseq: v?(\d+\.\d+\.\d+)", versions).group(1)
rows = {"fasta": params["fasta"], "fasta_md5": md5(params["fasta"]), "gtf": params["gtf"],
        "gtf_md5": md5(params["gtf"]), "gencode": str(params.get("gencode", False)).lower(),
        "rnaseq": release, "star": tool("STAR_GENOMEGENERATE", "star"),
        "salmon": tool("SALMON_INDEX", "salmon"), "source": results}
with open(out, "w") as f:
    f.writelines(f"{k}\t{v}\n" for k, v in rows.items())
print("".join(f"{k}\t{v}\n" for k, v in rows.items()), end="")
PY
cp -a "$R/genome/index/star" "$DEST.part/star"
cp -a "$R/genome/index/salmon" "$DEST.part/salmon"
mv "$DEST.part" "$DEST"
echo "saved: $DEST"
du -sh "$DEST"
