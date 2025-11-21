#!/usr/bin/python3
import gzip

A = "1KG.EUR.chr1.vcf.gz"
B = "DILI.cases.chr1.vcf.gz"
OUT_A = A.replace(".vcf.gz", ".common.vcf.gz")
OUT_B = B.replace(".vcf.gz", ".common.vcf.gz")

def ids_from_vcf(path):
    s = set()
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"): 
                continue
            c = line.rstrip("\n").split("\t")
            s.add(f"{c[0]}:{c[1]}:{c[3]}:{c[4]}")
    return s

common = ids_from_vcf(A) & ids_from_vcf(B)

def filter_vcf(in_path, out_path, keep):
    with gzip.open(in_path, "rt") as inp, gzip.open(out_path, "wt") as outp:
        for line in inp:
            if line.startswith("#"):
                outp.write(line); continue
            c = line.rstrip("\n").split("\t")
            uid = f"{c[0]}:{c[1]}:{c[3]}:{c[4]}"
            if uid in keep:
                outp.write(line)

filter_vcf(A, OUT_A, common)
filter_vcf(B, OUT_B, common)
print(f"Kept {len(common)} common IDs.\nWrote:\n  {OUT_A}\n  {OUT_B}")
