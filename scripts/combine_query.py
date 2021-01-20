import os
from script_utils import show_output
import pandas as pd


def main(s):
    """
    wrapped into function lest module be called by each subprocess
    """

    c = s.config
    w = s.wildcards
    i = s.input
    o = s.output
    cc = c["ascat"]

    show_output(
        f"Collecting chromosomal position data for {w.sample}_{w.type}",
        time=True,
    )
    pos_dfs = []
    for pos_file in i:
        pos_df = pd.read_csv(pos_file, sep="\t", compression="gzip", header=None)
        pos_dfs.append(pos_df)

    allpos_df = pd.concat(pos_dfs)
    allpos_df.to_csv(o.pos, sep="\t", index=False, compression="gzip", header=None)

    show_output(
        f"Finished combining data for {w.sample}_{w.type}. Written to {o.pos}",
        color="success",
    )


if __name__ == "__main__":
    main(snakemake)