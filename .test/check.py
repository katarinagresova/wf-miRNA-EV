"""Check the results of the test data set against .test/data/truth.tsv (section 8.2 of PLAN.md).

- Each insert of trim.min_len nt or longer gets its expected class in step 5, for each miRNA database.
- Each shorter insert is not in the collapsed reads.
- The off-target check gives each selected rRNA insert the result rRNA or tied.

Run in the environment of workflow/envs/py.yaml, after the workflow:

    python .test/check.py .test
"""

import gzip
import sys
from pathlib import Path

import pandas as pd


def collapsed(path):
    """Sequence -> rank."""
    with gzip.open(path, "rt") as fh:
        return {next(fh).rstrip(): int(header[2 : header.index("_x")]) for header in fh}


def main():
    root = Path(sys.argv[1])
    truth = pd.read_csv(root / "data/truth.tsv", sep="\t", keep_default_na=False)
    samples = pd.read_csv(root / "config/samples.tsv", sep="\t")["sample"]
    dbs = [c.removeprefix("expected_") for c in truth.columns if c.startswith("expected_")]
    errors, checked = [], 0
    for sample in samples:
        ranks = collapsed(root / f"results/collapse/{sample}/collapsed.fa.gz")
        classes = pd.read_csv(root / f"results/classify/{sample}/read_classes.tsv.gz", sep="\t", index_col="rank")
        t = truth[truth[sample] > 0]
        for row in t.itertuples():
            if row.expected_mirbase == "removed":
                checked += 1
                if row.sequence in ranks:
                    errors.append(f"{sample}: {row.detail} ({row.sequence}) is in the collapsed reads")
                continue
            if row.sequence not in ranks:
                errors.append(f"{sample}: {row.detail} ({row.sequence}) is not in the collapsed reads")
                continue
            for db in dbs:
                checked += 1
                expected, got = getattr(row, f"expected_{db}"), classes.loc[ranks[row.sequence], f"class_{db}"]
                if got != expected:
                    errors.append(f"{sample} {db}: {row.origin} {row.detail} ({row.sequence}): {got}, expected {expected}")
        for db in dbs:
            off = pd.read_csv(root / f"results/mirna/{sample}/off_target_{db}.tsv.gz", sep="\t")
            result = off.groupby("sequence")["result"].agg(lambda r: ",".join(sorted(set(r))))
            for row in t[t["off_target"].str.split(",").apply(lambda d: db in d)].itertuples():
                checked += 1
                got = result.get(row.sequence, "no miRNA-only hit")
                if got not in ["rRNA", "tied"]:
                    errors.append(f"{sample} {db}: off-target {row.detail} ({row.sequence}): {got}, expected rRNA or tied")
    print(f"{checked} checks, {len(errors)} errors")
    print("\n".join(errors))
    sys.exit(1 if errors else 0)


if __name__ == "__main__":
    main()
