# organizing previous code for search space calculation
# 19.12.2021
# Lin Chung-wen

## setup environment
library(openxlsx)

## load data
horf <- read.xlsx("~/Documents/INET-work/references/HuRI_binaryPPI/ORF_searchSpace.xlsx", sheet = 1)
horf_ad_sym <- unique(horf[horf$in_ad_space == 1, "ensembl_geneName"])
horf_ad_id <- unique(horf[horf$in_ad_space == 1, "orf_id"])
horf_db_sym <- unique(horf[horf$in_db_space == 1, "ensembl_geneName"])
horf_db_id <- unique(horf[horf$in_db_space == 1, "orf_id"])
horf_union_sym <- unique(horf$ensembl_geneName)
horf_union_id <- unique(horf$orf_id)

fbfg <- read.xlsx("~/Documents/INET-work/references/HuRI_binaryPPI/ORF_searchSpace.xlsx", sheet = 2)
fbfg_sym <- unique(fbfg$ensembl_geneName)
fbfg_id <- unique(fbfg$ORF_id)

intersect_space <- unique(horf_union_id[horf_union_id %in% fbfg_id])

union_space <- unique(c(horf_union_id, fbfg_id))

viralORF <- read.xlsx("/Volumes/GoogleDrive/My Drive/Paper_VirHostome_CoV2/04_Supplementary Information/Supplementary_Table_1.xlsx", sheet = "1a - Viral ORFs", startRow = 4)
viralORF2 <- viralORF[c(1:as.numeric(table(!is.na(viralORF$Index))["TRUE"])), ]
his3_preyAD_viral <- viralORF2$Viral.protein[viralORF2$Y2HHIS3 == 1]
his3_baitDB_viral <- viralORF2$Viral.protein[viralORF2$Y2HHIS3 == 1]
gfp_viral <- viralORF2$Viral.protein[viralORF2$Y2HGFP == 1]

intersect_viral <- viralORF2[viralORF2$Y2HHIS3 == 1 & viralORF2$Y2HGFP == 1, "Viral.protein"]

union_viral <- unique(c(his3_preyAD_viral, gfp_viral))

## do math
gfp_screen <- length(gfp_viral) * length(fbfg_id)
his3_ad_screen <- length(his3_preyAD_viral) * length(horf_db_id)
his3_db_screen <- (length(his3_baitDB_viral) - 1) * length(horf_ad_id)
intersect_screen <- length(intersect_viral) * length(intersect_space)
union_screen <- length(union_viral) * length(union_space)

gfp_viral_screen <- length(gfp_viral) * length(gfp_viral)
his3_viral_screen <- length(his3_preyAD_viral) * (length(his3_baitDB_viral) - 1)
intersect_viral_screen <- length(intersect_viral) * length(intersect_viral)
union_viral_screen <- length(union_viral) * length(union_viral)

message(rep("=", 60), "\n",
    "GFP search space: ", length(fbfg_id), "\n",
    "HIS3 AD preys search space: ", length(horf_ad_id), "\n",
    "HIS3 DB baits search space: ", length(horf_db_id), "\n",
    "HIS3 union search space: ", length(horf_union_id), "\n",
    "Intersect search space: ", length(intersect_space), "\n",
    "Union search space: ", length(union_space), "\n",

    "GFP viral ORF: ", length(gfp_viral), "\n",
    "HIS3 AD viral ORF: ", length(his3_preyAD_viral), "\n",
    "HIS3 DB viral ORF: ", length(his3_baitDB_viral) - 1, "\n",
    rep("=", 60), "\n",
    "GFP screening space: ", gfp_screen, "\n",
    "HIS3 screening space, vPrey (AD): ", his3_ad_screen, "\n",
    "HIS3 screening space, vBait (DB): ", his3_db_screen, "\n",
    "Intersect screening space: ", intersect_screen, "\n",
    "Union screening space: ", union_screen, "\n",
    "Fraction of viral ORF vs Human ORF: ", round(intersect_screen / union_screen, 4) * 100, "%", "\n",
    rep("=", 60), "\n",
    "GFP viral ORF screening space: ", gfp_viral_screen, "\n",
    "HIS3 viral ORF screening space: ", his3_viral_screen, "\n",
    "Intersect viral ORF screening space: ", intersect_viral_screen, "\n",
    "Union viral ORF screening space: ", union_viral_screen, "\n",
    "Fraction of intersect viral ORF / union viral ORF: ", round(intersect_viral_screen / union_viral_screen, 4) * 100, "%", "\n",
    rep("=", 60), "\n")
