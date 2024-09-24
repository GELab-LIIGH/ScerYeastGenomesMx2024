library(data.table)
library(ggplot2)
#################################PARSE ARGUMENTS################################
argsVal <- commandArgs(trailingOnly = T)
in_file <- argsVal[1] #Input file
freq_bars <- argsVal[2] #Output file
################################################################################

df = fread(in_file)
df_melt = melt(df,id.vars = c("ID","PhyloGroup_SACE469"))
df_ag = aggregate(.~PhyloGroup_SACE469+variable,data=df_melt[,c(-1)],FUN=sum)
df_ag = df_ag[df_ag$value!=0,]
pop_sizes = c()
pop_groups = unique(df$PhyloGroup_SACE469)
for(g in pop_groups){
  pop_sizes = c(pop_sizes,sum(df$PhyloGroup_SACE469==g))
}
pop = data.frame(PhyloGroup_SACE469 = pop_groups,size=pop_sizes)
df_freq = merge(df_ag,pop,by="PhyloGroup_SACE469")
df_freq$freq = df_freq$value/df_freq$size
heights <- aggregate(.~PhyloGroup_SACE469+freq,data=df_freq[,c(1,2,5)],FUN=length)
csvFile <- paste0(strsplit(freq_bars,split="pdf")[[1]],"csv")
fwrite(heights,file = csvFile)
bars <- ggplot(df_freq,aes(freq)) +
    geom_bar(width=0.02) +
    theme(axis.text.y = element_text(size=3),axis.text.x=element_text(angle = 90, size=3, hjust=1, vjust = 0.5))+
    theme(strip.text.x = element_text(size = 4))+
    theme(axis.title =element_text(size=6)) + theme_bw() +
    ylab("Introgressed genes") + xlab("Frequency") +
    facet_grid(PhyloGroup_SACE469 ~ .,scales="free_y")
#####
ggsave(freq_bars,plot=bars,width = 5, height = 9, units = "cm",bg = "white", dpi = 300)
