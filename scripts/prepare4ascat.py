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
        f"Combining {i.tumor} and {i.normal}",
        time=True,
    )

    # get tumor and normal posfile into df
    tumor_df = pd.read_csv(i.tumor, sep="\t")
    normal_df = pd.read_csv(i.tumor, sep="\t")

    merged_df = tumor_df.merge(normal_df, on=["chrom", "pos"], suffixes=("_T", "_N"))
    merged_df["pos"] = merged_df["pos"].astype(int)

    results_df.to_csv(o.pos, sep="\t", index=False)

    show_output(
        f"Finished combining data for {w.sample}_{w.type}. Written to {o.pos}",
        color="success",
    )


if __name__ == "__main__":
    main(snakemake)