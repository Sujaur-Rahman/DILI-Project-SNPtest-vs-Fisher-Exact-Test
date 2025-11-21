#mkdir -p ../results
#cd ../results

SNPTEST=/media/iscsi13T/home/mohammad/new_analysis_wes_4datasets/snptest_tool/snptest_v2.5.2_linux_x86_64_dynamic/snptest_v2.5.2

"$SNPTEST" \
  -data ./DILI.cases.chr1.gen.gz ./DILI.cases.chr1.bin.sample \
        ./1KG.EUR.chr1.gen.gz   ./1KG.EUR.chr1.bin.sample \
  -o snptest.chr1.out \
  -pheno bin1 \
  -frequentist 1 \
  -method score \
  -use_raw_phenotypes \
  -missing_code NA
