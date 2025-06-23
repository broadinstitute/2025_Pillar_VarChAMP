import seaborn as sns

clinvar_category = ['1_Pathogenic', '2_Benign', '3_Conflicting', '4_VUS', '5_Others', '6_No_ClinVar']
clinvar_palette_clinvar = sns.color_palette("Set2")
clinvar_palette_clinvar[0], clinvar_palette_clinvar[1], clinvar_palette_clinvar[2], clinvar_palette_clinvar[3], clinvar_palette_clinvar[4], clinvar_palette_clinvar[5], clinvar_palette_clinvar[7] = \
clinvar_palette_clinvar[1], clinvar_palette_clinvar[4], clinvar_palette_clinvar[3], clinvar_palette_clinvar[2], clinvar_palette_clinvar[0], clinvar_palette_clinvar[7], clinvar_palette_clinvar[5]
palette_dict = {
    "clinvar_clnsig_clean": dict(zip(clinvar_category, clinvar_palette_clinvar[:7]))
}