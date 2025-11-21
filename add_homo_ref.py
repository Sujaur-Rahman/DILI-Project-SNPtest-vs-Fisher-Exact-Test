#!/usr/bin/python3
import argparse, gzip

TARGET = "./.:.:.:.:."  # exact string to replace

def make_replacement(format_str):
    """Return a sample string with GT=0/0 and the rest '.' matching FORMAT length."""
    n_fields = len(format_str.split(":"))
    if n_fields <= 1:
        return "0/0"
    return "0/0" + ":" + ":".join(["."] * (n_fields - 1))

def convert(in_path, out_path):
    with gzip.open(in_path, "rt") as fin, gzip.open(out_path, "wt") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 9:
                fmt = parts[8]
                repl = make_replacement(fmt)
                for i in range(9, len(parts)):
                    if parts[i] == TARGET:
                        parts[i] = repl
                line = "\t".join(parts) + "\n"
            fout.write(line)

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Replace ./.:.:.:.:. with 0/0 (keeping FORMAT length).")
    ap.add_argument("-i", "--input", default="DILI.cases.chr1.common.vcf.gz", help="Input VCF.gz")
    ap.add_argument("-o", "--output", default="DILI.cases.chr1.shared.00.vcf.gz", help="Output VCF.gz")
    args = ap.parse_args()
    convert(args.input, args.output)
