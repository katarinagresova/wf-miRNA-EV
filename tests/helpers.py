"""Helpers for the unit tests: write small input tables by hand and run the workflow scripts through their command
line, as the rules do."""

import gzip
import random
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "scripts"

# Arms of the test hairpins: precursor, mature, start, end. mir-c1 and mir-c2 are paralogs with the same mature name.
ARMS = [
    ("mir-a", "miR-a-5p", 5, 27),
    ("mir-a", "miR-a-3p", 45, 67),
    ("mir-b", "miR-b-5p", 5, 27),
    ("mir-c1", "miR-c-5p", 5, 27),
    ("mir-c2", "miR-c-5p", 5, 27),
]
# Step 5 settings of config/config.yaml
STEP5 = {
    "min-len": 16,
    "max-len": 28,
    "start-tol": 3,
    "end-tol": 5,
    "contaminant-min-len": 18,
    "end-len": 25,
    "end-len-exact": 18,
}
SMALL_REFS = ["mirbase", "rrna", "trna", "ncrna"]
COUNT = 10  # reads of each collapsed sequence


def run_script(name, *args):
    """Run workflow/scripts/<name>; fail the test with the script output if it fails."""
    p = subprocess.run([sys.executable, SCRIPTS / name, *map(str, args)], capture_output=True, text=True)
    assert p.returncode == 0, p.stdout + p.stderr
    return p.stdout


def write_tsv(path, rows, header=None):
    """Write rows as TSV, with gzip if the name ends with .gz."""
    lines = ([header] if header else []) + list(rows)
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as fh:
        fh.writelines("\t".join(map(str, r)) + "\n" for r in lines)
    return path


def sequence(rank, length):
    """A fixed random sequence for each read."""
    return "".join(random.Random(rank).choices("ACGT", k=length))


def options(settings):
    return [x for k, v in settings.items() for x in (f"--{k}", v)]


def classify(tmp_path, reads, hits=None, contaminant=(), genome=(), ends=()):
    """Run quantify_mirna.py (step 5) for the database mirbase, with the ARMS hairpins.

    reads: rank -> read length (COUNT reads each). hits: reference -> rows (rank, reference, offset, length,
    mismatches); contaminant: the same rows plus the strand; genome: rows (rank, mismatches, classes); ends: rows
    (end, index, reference, mismatches). Reference names are <id>|<class>|<name>. A table without rows is an empty
    file, as bowtie writes it.
    """
    d = tmp_path / "step5"
    d.mkdir()
    with gzip.open(d / "collapsed.fa.gz", "wt") as fh:
        for rank, length in sorted(reads.items()):
            fh.write(f">s{rank}_x{COUNT}\n{sequence(rank, length)}\n")
    hits = hits or {}
    hit_files = [write_tsv(d / f"{ref}.hits.tsv.gz", hits.get(ref, [])) for ref in SMALL_REFS]
    write_tsv(d / "mirbase_arms.tsv", ARMS, header=["precursor", "mature", "start", "end"])
    run_script(
        "quantify_mirna.py",
        "--reads", d / "collapsed.fa.gz",
        "--hits", *hit_files,
        "--contaminant", write_tsv(d / "contaminant.hits.tsv.gz", contaminant),
        "--genome", write_tsv(d / "genome.classes.tsv.gz", genome),
        "--ends", write_tsv(d / "unmapped_ends.hits.tsv.gz", ends),
        "--arms", d / "mirbase_arms.tsv",
        *options(STEP5),
        "--top-n", 10,
        "--out-reads", d / "read_classes.tsv.gz",
        "--out-composition", d / "composition.tsv",
        "--out-counts", d / "mirna_mirbase.tsv",
        "--out-mirna-reads", d / "mirna_mirbase.reads.tsv.gz",
        "--out-contaminants", d / "contaminants.tsv",
        "--out-contaminant-rules", d / "contaminant_rules.tsv",
        "--out-ends", d / "unmapped_ends.tsv",
        "--out-top", d / "top_hits.tsv",
    )
    return SimpleNamespace(
        dir=d,
        classes=pd.read_csv(d / "read_classes.tsv.gz", sep="\t", index_col="rank")["class_mirbase"],
        counts=pd.read_csv(d / "mirna_mirbase.tsv", sep="\t", index_col="mature"),
    )
