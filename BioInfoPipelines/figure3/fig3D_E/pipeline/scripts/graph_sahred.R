####################################LIBRARIES###################################
library(data.table)
library(ggplot2)
################################################################################
#####################################INPUTS#####################################
ori <- ""
folder <- "data/"
file <- paste0(ori,folder,"MA1MA2FgAlp_frequencies.csv")
################################################################################
#####################################OUTPUTS####################################
out_patterns <- paste0(ori,folder,"pattern_histogram.pdf")
out_freqs <- paste0(ori,folder,"frequencies.pdf")
out_freqs_norm <- paste0(ori,folder,"frequencies_norm.pdf")
csv_freqs_norm <- paste0(ori,folder,"frequencies_norm.csv")
out_freqs_normS <- paste0(ori,folder,"frequencies_norm_Singles.pdf")
out_freqs_normV <- paste0(ori,folder,"frequencies_norm_Various.pdf")
csvOut <- paste0(ori,folder,"histogram.csv")
################################################################################
##############################CONSTANTS & FUNCTIONS#############################

################################################################################
#######################################CODE#####################################
freqs <- fread(file,header=TRUE)
colnames(freqs)[1] <- "ID"
freqs[,`:=`(code = vapply(strsplit(ID,"_"), `[`, 2, FUN.VALUE=character(1)) )]
freqs[,`:=`(code2 = vapply(strsplit(code,"G"), `[`, 1, FUN.VALUE=character(1)) )]
freqs[,`:=`(Chr = paste("Chr",code2,sep="_") )]
freqs$Chr[freqs$Chr == "Chr_NA"] <- "Mit"
freqs_melt <- melt(freqs[,c(1,2,3,4,5,8)],id.vars=c("ID","Chr"))
freqs_melt[Chr=="Chr_03"&value>0,]
######
#PLOT#
p_freqs <- ggplot(freqs_melt,aes(x=variable,y=value,group=ID)) +
  geom_point() + geom_line() +
  theme_bw() + guides(color = "none") +
  theme(axis.text = element_text(size=10), axis.title =element_text(size=14),
    axis.text.x = element_text(angle=90)) +
  ylab("Frquency in population") + xlab("S. cerevisiae populations") +
  facet_wrap(Chr ~ .,ncol=4)
ggsave(out_freqs,plot=p_freqs,height=24,width=26,units="cm",dpi=300)

######
abs <- freqs
abs$Tamaulipas <- as.integer(abs$Tamaulipas > 0.000)
abs$MMMG <- as.integer(abs$MMMG > 0.000)
abs$FG <- as.integer(abs$FG > 0.000)
abs$Alpechin <- as.integer(abs$Alpechin > 0.000)
abs[,`:=`(pattern = paste0(Tamaulipas,MMMG,FG,Alpechin) )]
######
#PLOT#
p_abs <- ggplot(abs,aes(x=pattern)) +
  geom_histogram(stat="count") +
  theme_bw() + ggtitle("Codification : MA1|MA2|FG|Alp") +
  theme(axis.text = element_text(size=10), axis.title =element_text(size=12),
    axis.text.x = element_text(angle=90)) +
  ylab("Count") + xlab("Observed Patterns")
ggsave(out_patterns,plot=p_abs,height=12,width=14,units="cm",dpi=300)
######

abs <- abs[abs$Chr != "Mit",]
norml <- as.data.frame(table(abs$pattern))
#Alp
cccu <- nrow(abs[abs$Alpechin==1,"ID"])
#FG
ccuc <- nrow(abs[abs$FG==1,"ID"])
#MA1
uccc <- nrow(abs[abs$Tamaulipas==1,"ID"])
#MA2
cucc <- nrow(abs[abs$MMMG==1,"ID"])
#MA1 U MA2
uucc <- min(nrow(abs[abs$Tamaulipas==1,"ID"]),nrow(abs[abs$MMMG==1,"ID"]))
#MA1 U FG
ucuc <- min(nrow(abs[abs$Tamaulipas==1,"ID"]),nrow(abs[abs$FG==1,"ID"]))
#MA1 U Alp
uccu <- min(nrow(abs[abs$Tamaulipas==1,"ID"]),nrow(abs[abs$Alpechin==1,"ID"]))
#MA2 U Fg
cuuc <- min(nrow(abs[abs$MMMG==1,"ID"]),nrow(abs[abs$FG==1,"ID"]))
#MA1 U Alp
cucu <- min(nrow(abs[abs$MMMG==1,"ID"]),nrow(abs[abs$Alpechin==1,"ID"]))
#FG U Alp
ccuu <- min(nrow(abs[abs$FG==1,"ID"]),nrow(abs[abs$Alpechin==1,"ID"]))
#MA1 U MA2 U FG
uuuc <- min(nrow(abs[abs$Tamaulipas==1,"ID"]), nrow(abs[abs$MMMG==1,"ID"]),nrow(abs[abs$FG==1,"ID"]))
#MA1 U MA2 U Alp
uucu <- min(nrow(abs[abs$Tamaulipas==1,"ID"]), nrow(abs[abs$MMMG==1,"ID"]),nrow(abs[abs$Alpechin==1,"ID"]))
#MA2 U FG U Alp
cuuu <- min(nrow(abs[abs$MMMG==1,"ID"]), nrow(abs[abs$FG==1,"ID"]),nrow(abs[abs$Alpechin==1,"ID"]))
# MA1 U MA2 U FG U Alp
uuuu <- nrow(abs[abs$FG==1,"ID"])

#
facts <- c(cccu,ccuc,ccuu,cucc,cucu,cuuc,cuuu,uccc,uccu,ucuc,uucc,uucu,uuuc,uuuu)
facts_order <- c("1000","0100","0010","0001","1100","1010","1001","0110","0101","0011","1110","1101","0111","1111")

norml$normfact <- facts
norml$nor <- norml$Freq/norml$normfact
norml$Var1 <- factor(norml$Var1,levels=facts_order)
######
#CSV#
fwrite(norml,csvOut)
#####
#PLOT#
p_norm <- ggplot(norml,aes(x=Var1,y=nor)) +
  geom_histogram(stat="identity") +
  theme_bw() + ggtitle("Codification : MA1|MA2|FG|Alp") +
  theme(axis.text = element_text(size=10), axis.title =element_text(size=12),
    axis.text.x = element_text(angle=90)) +
  ylab("Frequency relative to maximum number of possible shared genes") + xlab("Observed Patterns")
p_norm
ggsave(out_freqs_norm,plot=p_norm,height=14,width=16,units="cm",dpi=300)
fwrite(norml,csv_freqs_norm)

normlS <- norml[norml$Var1 == "0001" | norml$Var1 == "0010" | norml$Var1 == "0100" | norml$Var1 == "1000",]
p_norm_sing <- ggplot(normlS,aes(x=Var1,y=nor)) +
  geom_histogram(stat="identity") +
  theme_bw() + ggtitle("MA1|MA2|FG|Alp") +
  theme(axis.text = element_text(size=10), axis.title =element_text(size=12),
    axis.text.x = element_text(angle=90)) +
  ylab("Frequency relative to maximum number of possible shared genes") + xlab("Observed Patterns")

ggsave(out_freqs_normS,plot=p_norm_sing,height=14,width=5.3,units="cm",dpi=300)


normlV <- norml[norml$Var1 != "0001" & norml$Var1 != "0010" & norml$Var1 != "0100" & norml$Var1 != "1000",]
p_norm_sing <- ggplot(normlV,aes(x=Var1,y=nor)) +
  geom_histogram(stat="identity") +
  theme_bw() + ggtitle("MA1|MA2|FG|Alp") +
  theme(axis.text = element_text(size=10), axis.title =element_text(size=12),
    axis.text.x = element_text(angle=90)) +
  ylab("Frequency relative to maximum number of possible shared genes") + xlab("Observed Patterns")

ggsave(out_freqs_normV,plot=p_norm_sing,height=14,width=11.7,units="cm",dpi=300)

######

################################################################################
