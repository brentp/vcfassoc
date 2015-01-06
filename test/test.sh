set -ex
# steal
NVCF=test.auto_dom.no_parents.2.vcf
VCF=https://raw.githubusercontent.com/arq5x/gemini/master/test/$NVCF
wget $VCF

python vcfassoc.py $NVCF test/covs.txt "affected ~ I(genotype > 0)" --groups family | cut -f 1-9
python vcfassoc.py $NVCF test/covs.txt "affected ~ genotype" --groups family | cut -f 1-9
rm $NVCF

RVCF=https://raw.githubusercontent.com/arq5x/gemini/master/test/test.auto_rec.vcf

python vcfassoc.py $RVCF test/covs.txt "affected ~ I(genotype == 2)" --groups family | cut -f 1-10
