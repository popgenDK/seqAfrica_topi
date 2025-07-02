BFILE="/maps/projects/seqafrica/people/mdn487/Topi/goatBased/Goat_topi_variable_sites_mergeblack_wildebeest-blue_wildebeest-hartebeest"
OUTDIR="/projects/alab/people/vzw531/Topi/ibs"
OUT="Goat_topi_variable_biallelic_noindels_polymorphic_goodSites_outgroup"
PLINK="/projects/seqafrica/apps/modules/software/plink/1.90b6.21/bin/plink"

$PLINK -bfile $BFILE --distance square ibs allele-ct --out $OUTDIR/$OUT --allow-extra-chr --chr-set 29
