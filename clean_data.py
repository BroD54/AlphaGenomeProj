"""
Filter the ClinVar web-download TSV at data/raw.tsv down to a clean APC variant table
for AlphaGenome.

Usage:
    python clean_data.py

Outputs:
    apc_variants.csv  all exact "Pathogenic" / "Benign" variants that pass cleaning
"""

import io
import pandas as pd

INPUT_PATH = "data/raw.tsv"
OUTPUT_PATH = "data/apc_variants.csv"
GENE = "APC"
CHROM = "5"
MAIN_LABELS = {"Pathogenic", "Benign"}


def load(path):
    # ClinVar's TSV has a trailing extra tab on every data row (always an empty final
    # field); strip it so field counts match the header and pandas doesn't mis-align columns.
    with open(path) as f:
        text = "\n".join(line.rstrip("\t\n") for line in f)
    df = pd.read_csv(io.StringIO(text), sep="\t", dtype=str, keep_default_na=False)
    df.columns = df.columns.str.strip()
    df = df.apply(lambda col: col.str.strip())
    return df.drop_duplicates(subset="VariationID")


def parse_spdi(spdi):
    """NC_000005.10:112707231:G:A -> (112707232, 'G', 'A'). SPDI position is 0-based."""
    try:
        _, pos0, ref, alt = spdi.split(":")
        return int(pos0) + 1, ref, alt
    except (ValueError, AttributeError):
        return None, None, None


def clean(df):
    n0 = len(df)
    df = df[df["Germline classification"].isin(MAIN_LABELS)]              # exact Pathogenic/Benign only
    df = df[df["Gene(s)"] == GENE]                                        # only APC, not multi-gene CNVs
    df = df[df["Variant type"] == "single nucleotide variant"]
    df = df[df["GRCh38Chromosome"] == CHROM]

    # Drop low-confidence / conflicting review statuses
    rs = df["Germline review status"].str.lower()
    df = df[~rs.str.contains("no assertion|conflicting|no classification")]

    # ref / alt / 1-based position from Canonical SPDI
    parsed = df["Canonical SPDI"].apply(parse_spdi)
    df = df.assign(pos=[p[0] for p in parsed], ref=[p[1] for p in parsed], alt=[p[2] for p in parsed])
    df = df[df["ref"].str.len().eq(1) & df["alt"].str.len().eq(1)]

    # Sanity check: SPDI-derived position should match ClinVar's GRCh38Location
    mismatch = df["pos"].astype(str) != df["GRCh38Location"]
    if mismatch.any():
        print(f"WARNING: {mismatch.sum()} rows where SPDI position != GRCh38Location; dropping them")
        df = df[~mismatch]

    print(f"Rows after cleaning: {len(df)} (from {n0})")
    return df


def to_output(df):
    out = pd.DataFrame({
        "chrom": "chr" + CHROM,
        "pos": df["pos"].astype(int),
        "ref": df["ref"],
        "alt": df["alt"],
        "label": df["Germline classification"],
        "hgvs": [f"NC_000005.10:g.{p}{r}>{a}" for p, r, a in zip(df["pos"], df["ref"], df["alt"])],
        "name": df["Name"],
        "variation_id": df["VariationID"],
        "rsid": df["dbSNP ID"],
        "consequence": df["Molecular consequence"],
    })
    return out.sort_values("pos").reset_index(drop=True)


def main():
    df = clean(load(INPUT_PATH))

    print("\nLabel counts after cleaning:")
    print(df["Germline classification"].value_counts().to_string(), "\n")

    out = to_output(df)
    out.to_csv(OUTPUT_PATH, index=False)

    blank = (df["Molecular consequence"] == "").sum()
    if blank:
        print(f"\nNOTE: {blank} variants have a blank Molecular consequence in the download.")

    print(f"\nWrote apc_variants.csv ({len(out)} rows)")


if __name__ == "__main__":
    main()
