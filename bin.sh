# Controls: bin1 = 0
awk 'NR==1{print $0,"bin1"; next} NR==2{print $0,"B"; next} NR>2{print $0,"0"}' \
  1KG.EUR.chr1.samples > 1KG.EUR.chr1.bin.sample

# Cases: bin1 = 1
awk 'NR==1{print $0,"bin1"; next} NR==2{print $0,"B"; next} NR>2{print $0,"1"}' \
  DILI.cases.chr1.samples > DILI.cases.chr1.bin.sample
