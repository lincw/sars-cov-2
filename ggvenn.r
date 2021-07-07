# trying different venn diagram package
# Lin Chung-wen

if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")

library(ggvenn)
x <- list(c1 = c("FDX2","ICAM1","ICAM3","ICAM4","ICAM5","IFNAR2","KAT7","LZTFL1","OAS1","OAS2","OAS3","RAVER1","TAC4","TYK2"),
    c2 = c("FDX2","FOXP4","ICAM1","ICAM3","ICAM4","ICAM5","IFNAR2","KANSL1","LINC02210-CRHR1","LRRC37A","LRRC37A2","LZTFL1","MAPT","NSF","OAS1","OAS2","OAS3","PLEKHM1","RAVER1","SPPL2C","STH","TMEM65","TYK2"),
    c3 = c("LZTFL1","NUCB1","NXPE3","OAS1","OAS2","OAS3","PLEKHA4","PPP1R15A","RPL24","SLC6A20")
)
x_ggvenn <- ggvenn(x, stroke_size = 0.5, set_name_size = 4)

