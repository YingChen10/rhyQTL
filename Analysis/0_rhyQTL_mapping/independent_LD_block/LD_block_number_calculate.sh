dir="/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/13_Results"
indir="${dir}/00_all_rhyQTL"
outdir="${dir}/09_LD_block"

# generate the header for the LD block number file
echo -e "Tissue\tNumber" > ${outdir}/LD_block_number.txt

# format ldetect data which was downloaded from https://bitbucket.org/nygcresearch/ldetect-data/downloads/
ldetect_data="/workspace/rsrch1/ychen/Projects/common_data/ldetect-data/ldetect/EUR/fourier_ls-all.bed"
formated_ldetect_data="/workspace/rsrch1/ychen/Projects/common_data/ldetect-data/ldetect/EUR/sorted_fourier_ls-all.bed"
tail -n +2 ${ldetect_data}| sort -k1,1V -k2,2n --field-separator=$'\t' | awk 'OFS="\t" {print $1, $2, $3,$1"_"$2"_"$3}' > ${formated_ldetect_data}

tissues=("Adipose-Subcutaneous" "Adipose-Visceral_Omentum" "AdrenalGland" "Artery-Aorta" "Artery-Coronary"
"Artery-Tibial" "Brain-Amygdala" "Brain-Anteriorcingulatecortex_BA24" "Brain-Caudate_basalganglia"
"Brain-CerebellarHemisphere" "Brain-Cerebellum" "Brain-Cortex" "Brain-FrontalCortex_BA9" "Brain-Hippocampus"
"Brain-Hypothalamus" "Brain-Nucleusaccumbens_basalganglia" "Brain-Putamen_basalganglia" "Brain-Spinalcord_cervicalc-1"
"Brain-Substantianigra" "Breast-MammaryTissue" "Colon-Sigmoid" "Colon-Transverse" "Esophagus-GastroesophagealJunction"
"Esophagus-Mucosa" "Esophagus-Muscularis" "Heart-AtrialAppendage" "Heart-LeftVentricle" "Liver" "Lung" "MinorSalivaryGland"
"Muscle-Skeletal" "Nerve-Tibial" "Ovary" "Pancreas" "Pituitary" "Prostate" "Skin-NotSunExposed_Suprapubic"
"Skin-SunExposed_Lowerleg" "SmallIntestine-TerminalIleum" "Spleen" "Stomach" "Testis" "Thyroid" "Uterus" "Vagina")

mkdir -p "${outdir}/rhySNP_bed"
mkdir -p "${outdir}/rhySNP_intersect_bed"

for tissue in "${tissues[@]}"
do
    echo "Running for tissue: $tissue"

# generate rhyQTL bed file
tail -n +2 ${indir}/${tissue}.txt|awk '{split($3, arr, "_"); print arr[1], arr[2], arr[2]+1, $1, $2, $3}'|sort -k1,1V -k2,2n --field-separator=$'\t' |awk 'OFS="\t" {print $1, $2, $3, $4, $5, $6}' > ${outdir}/rhySNP_bed/${tissue}.bed

# intersct with LD block data
bedtools intersect -a ${formated_ldetect_data} -b "${outdir}/rhySNP_bed/${tissue}.bed" -wa -wb > "${outdir}/rhySNP_intersect_bed/${tissue}.bed"

# count the number of LD blocks
cut -f 4,8,9 "${outdir}/rhySNP_intersect_bed/${tissue}.bed" | sort | uniq | wc -l | awk -v tissue="$tissue" '{print tissue "\t" $1}' >> ${outdir}/LD_block_number.txt

done
