library(msa)
setwd("~/workplace/hCoV")
orf1ab_nt <- readDNAStringSet("orf1ab_nt.fa")
orf1ab_pp <- readAAStringSet("orf1ab_pp.fa")
aa <- readAAStringSet("aa.fa")
nt <- readDNAStringSet("nt.fa")

######################################################################
# SARS2 vs 229E
# nsp5
nsp5 <- c("SARS2-nsp5", "h229E-nsp5")
nsp5_nt <- nt[nsp5, ]
nsp5_aa <- aa[nsp5, ]
nsp5_nt_msa <- msa(nt[nsp5, ], method = "ClustalW")
nsp5_aa_msa <- msa(aa[nsp5, ], method = "ClustalW")

msaPrettyPrint(nsp5_nt_msa, 
    alFile = "SARS2-229e_nsp5_nt_msa.fasta",
    shadingMode = "identical", 
    askForOverwrite = FALSE)
msaPrettyPrint(nsp5_aa_msa,
    alFile = "SARS2-229e_nsp5_aa_msa.fasta",
    shadingMode = "identical", 
    askForOverwrite = FALSE)

# nsp8
nsp8 <- c("SARS2-nsp8", "h229E-nsp8")
nsp8_nt <- nt[nsp8, ]
nsp8_aa <- aa[nsp8, ]
nsp8_nt_msa <- msa(nt[nsp8, ], method = "ClustalW")
nsp8_aa_msa <- msa(aa[nsp8, ], method = "ClustalW")

msaPrettyPrint(nsp8_nt_msa, 
    alFile = "SARS2-229e_nsp8_nt_msa.fasta",
    shadingMode = "identical", 
    askForOverwrite = FALSE)
msaPrettyPrint(nsp8_aa_msa, 
    alFile = "SARS2-229e_nsp8_aa_msa.fasta",
    shadingMode = "identical", 
    askForOverwrite = FALSE)

# nsp12
nsp12 <- c("SARS2-nsp12", "h229E-nsp12")
nsp12_nt <- nt[nsp12, ]
nsp12_aa <- aa[nsp12, ]
nsp12_nt_msa <- msa(nt[nsp12, ], method = "ClustalW")
nsp12_aa_msa <- msa(aa[nsp12, ], method = "ClustalW")

msaPrettyPrint(nsp12_nt_msa, 
    alFile = "SARS2-229e_nsp12_nt_msa.fasta",
    shadingMode = "identical", 
    askForOverwrite = FALSE)
msaPrettyPrint(nsp12_aa_msa, 
    alFile = "SARS2-229e_nsp12_aa_msa.fasta",
    shadingMode = "identical", 
    askForOverwrite = FALSE)

######################################################################
# SARS vs 229E
# nsp5
nsp5 <- c("SARS-nsp5", "h229E-nsp5")
nsp5_nt <- nt[nsp5, ]
nsp5_aa <- aa[nsp5, ]
nsp5_nt_msa <- msa(nt[nsp5, ], method = "ClustalW")
nsp5_aa_msa <- msa(aa[nsp5, ], method = "ClustalW")

msaPrettyPrint(nsp5_nt_msa, 
    alFile = "SARS-229e_nsp5_nt_msa.fasta",
    shadingMode = "identical", 
    askForOverwrite = FALSE)
msaPrettyPrint(nsp5_aa_msa,
    alFile = "SARS-229e_nsp5_aa_msa.fasta",
    shadingMode = "identical", 
    askForOverwrite = FALSE)

# nsp8
nsp8 <- c("SARS-nsp8", "h229E-nsp8")
nsp8_nt <- nt[nsp8, ]
nsp8_aa <- aa[nsp8, ]
nsp8_nt_msa <- msa(nt[nsp8, ], method = "ClustalW")
nsp8_aa_msa <- msa(aa[nsp8, ], method = "ClustalW")

msaPrettyPrint(nsp8_nt_msa, 
    alFile = "SARS-229e_nsp8_nt_msa.fasta",
    shadingMode = "identical", 
    askForOverwrite = FALSE)
msaPrettyPrint(nsp8_aa_msa, 
    alFile = "SARS-229e_nsp8_aa_msa.fasta",
    shadingMode = "identical", 
    askForOverwrite = FALSE)

# nsp12
nsp12 <- c("SARS-nsp12", "h229E-nsp12")
nsp12_nt <- nt[nsp12, ]
nsp12_aa <- aa[nsp12, ]
nsp12_nt_msa <- msa(nt[nsp12, ], method = "ClustalW")
nsp12_aa_msa <- msa(aa[nsp12, ], method = "ClustalW")

msaPrettyPrint(nsp12_nt_msa, 
    alFile = "SARS-229e_nsp12_nt_msa.fasta",
    shadingMode = "identical", 
    askForOverwrite = FALSE)
msaPrettyPrint(nsp12_aa_msa, 
    alFile = "SARS-229e_nsp12_aa_msa.fasta",
    shadingMode = "identical", 
    askForOverwrite = FALSE)
######################################################################

select <- c("SARS2", "h229E")
orf1ab_nt_msa <- msa(orf1ab_nt[select, ])
msaPrettyPrint(orf1ab_nt_msa, 
    shadingMode = "identical", 
    askForOverwrite = FALSE,
    code = c("\\shadingmode[allmatchspecial]{similar}", "\\fingerprint{2500}", "\\shownumbering{right}", "\\showruler{1}{top}", "\\shadingcolors{grays}", "\\feature{bottom}{1}{8895..9801}{<=>}{h229E-nsp5}", "\\feature{top}{1}{9789..10708}{<=>}{SARS2-nsp5}"))

h229e_sub <- c(nt["h229E-nsp5", ], orf1ab_nt["h229E", ])
h229e_sub_msa <- msa(h229e_sub)
msaPrettyPrint(h229e_sub_msa, 
    shadingMode = "similar", 
    askForOverwrite = FALSE,
    code = c("\\shadingmode[allmatchspecial]{similar}", "\\fingerprint{2300}", "\\shownumbering{right}", "\\showruler{1}{top}", "\\shadingcolors{grays}", "\\feature{bottom}{1}{8895..9801}{,-,}{nsp5}"))

SARS2_sub <- c(nt["SARS2-nsp5", ], orf1ab_nt["SARS2", ])
SARS2_sub_msa <- msa(SARS2_sub)

# writeXStringSet(nsp5_nt, file = "nsp5_nt.fa", append = FALSE,
                    #  compress = FALSE, format SARS2-229e_= "fasta")