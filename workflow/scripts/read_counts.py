"""Raw read counts, read lengths and Q30 of each library."""

import argparse
from pathlib import Path

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--stats", nargs="+", required=True, help="seqkit stats -T -a output, one per sample")
    p.add_argument("--samples", required=True, help="config/samples.tsv")
    p.add_argument("--output", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    samples = pd.read_csv(args.samples, sep="\t", dtype=str).set_index("sample")

    ours = []
    for path in args.stats:
        st = pd.read_csv(path, sep="\t")
        ours.append(
            {
                "sample": Path(path).stem,
                "reads": int(st.loc[0, "num_seqs"]),
                "min_len": int(st.loc[0, "min_len"]),
                "max_len": int(st.loc[0, "max_len"]),
                "Q30_pct": float(st.loc[0, "Q30(%)"]),
            }
        )
    df = samples[["sid", "run", "exclude"]].join(pd.DataFrame(ours).set_index("sample"), how="right")

    df.to_csv(args.output, sep="\t")


if __name__ == "__main__":
    main()
