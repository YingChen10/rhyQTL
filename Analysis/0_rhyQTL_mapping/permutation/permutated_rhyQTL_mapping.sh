#!/bin/bash

# Define parameters
tissues=("Adipose-Subcutaneous" "Adipose-Visceral_Omentum" "AdrenalGland" 
"Artery-Aorta" "Artery-Coronary" "Artery-Tibial" "Brain-Amygdala" "Brain-Anteriorcingulatecortex_BA24" 
"Brain-Caudate_basalganglia" "Brain-CerebellarHemisphere" "Brain-Cerebellum" "Brain-Cortex" 
"Brain-FrontalCortex_BA9" "Brain-Hippocampus" "Brain-Hypothalamus" "Brain-Nucleusaccumbens_basalganglia" 
"Brain-Putamen_basalganglia" "Brain-Spinalcord_cervicalc-1" "Brain-Substantianigra" "Breast-MammaryTissue" 
"Colon-Sigmoid" "Colon-Transverse" "Esophagus-GastroesophagealJunction" "Esophagus-Mucosa" 
"Esophagus-Muscularis" "Heart-AtrialAppendage" "Heart-LeftVentricle" "Liver" "Lung" "MinorSalivaryGland" 
"Muscle-Skeletal" "Nerve-Tibial" "Ovary" "Pancreas" "Pituitary" "Prostate" "Skin-NotSunExposed_Suprapubic" 
"Skin-SunExposed_Lowerleg" "SmallIntestine-TerminalIleum" "Spleen" "Stomach" "Testis" "Thyroid" "Uterus" "Vagina")

files=("HLA_split_pos_aa" "HLA_split_pos_ab" "HLA_split_pos_ac" \
       "split_pos_aa" "split_pos_ab" "split_pos_ac" "split_pos_ad" \
       "split_pos_ae" "split_pos_af" "split_pos_ag" "split_pos_ah" \
       "split_pos_ai" "split_pos_aj" "split_pos_ak" "split_pos_al" \
       "split_pos_am" "split_pos_an" "split_pos_ao" "split_pos_ap" \
       "split_pos_aq" "split_pos_ar" "split_pos_as" "split_pos_at" \
       "split_pos_au" "split_pos_av" "split_pos_aw" "split_pos_ax" \
       "split_pos_ay" "split_pos_az" "split_pos_ba" "split_pos_bb" \
       "split_pos_bc" "split_pos_bd" "split_pos_be" "split_pos_bf" \
       "split_pos_bg" "split_pos_bh" "split_pos_bi" "split_pos_bj" \
       "split_pos_bk" "split_pos_bl" "split_pos_bm" "split_pos_bn" \
       "split_pos_bo" "split_pos_bp" "split_pos_bq" "split_pos_br" \
       "split_pos_bs" "split_pos_bt")

# Scripts to run
scripts=("0_0_genotype_permutation.R" "0_1_rhythm_regression.R" "0_2_rhythm_compare.R" "0_3_rhythm_compare_multiple_times.R" "0_4_hanova.R")

# Loop through combinations of tissue, files, and scripts
for tissue in "${tissues[@]}"; do
  for file in "${files[@]}"; do
    for script in "${scripts[@]}"; do
      echo "Running: Rscript $script $tissue $file"
      Rscript "$script" "$tissue" "$file"
    done
  done
done