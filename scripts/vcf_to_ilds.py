#!/usr/bin/env python3
import argparse

import numpy as np
import pandas as pd
from cyvcf2 import VCF


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert VCF to haplotype-style CSV with syn/nonsyn variants."
    )
    parser.add_argument(
        "--input",
        default="UBannotated.vcf",
        help="Input VCF annotated with ANN field (default: UBannotated.vcf).",
    )
    parser.add_argument(
        "--output",
        default="UB_haplotypes.csv",
        help="Output CSV filename (default: UB_haplotypes.csv).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    vcf = VCF(args.input)
    samples = vcf.samples

    rows = []
    for var in vcf:
        ann = var.INFO.get("ANN")
        if not ann:
            continue

        ann_fields = ann.split(",")[0].split("|")
        effect = ann_fields[1]
        gene_id = ann_fields[3]

        if effect in ("synonymous_variant",):
            site_type = "syn"
        elif effect in ("missense_variant", "stop_gained", "stop_lost"):
            site_type = "nonsyn"
        else:
            continue

        gts = []
        for gt in var.genotypes:
            if gt[0] is None or gt[1] is None:
                gts.append(np.nan)
            else:
                gts.append(1 if (gt[0] + gt[1]) > 0 else 0)

        rows.append([var.CHROM, gene_id, var.POS, site_type] + gts)

    cols = ["contig", "gene_id", "site_pos", "site_type"] + samples
    df = pd.DataFrame(rows, columns=cols)
    df.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
