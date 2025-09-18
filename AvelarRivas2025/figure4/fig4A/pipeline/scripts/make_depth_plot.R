library(data.table)
library(zoo)
library(ggplot2)

args=(commandArgs(TRUE))

if(length(args)==0){
    stop("No arguments supplied.")
}else{
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
}

if( getDTthreads() < 10 || getDTthreads() > 10){
  setDTthreads(threads = 10)
}else{print(paste("WARNING: data.table package is running with ", getDTthreads(), " threads.", sep=''))}

plus_10 <- function(y){
    n = length(y)
    pt = 0
    for(x in y){
        if(x > 10){
            pt = pt + 1
        }
    }
    return(format(round((pt/n)*100, 4),nsmall=4))
}

PerContigPATH = paste(STATS,SAMPLE,"_perContig.csv",sep="")
PerRefPATH = paste(STATS,SAMPLE,"_perReference.csv",sep="")
WindowsPATH = paste(STATS,SAMPLE,"_windows.csv",sep="")

depth <- fread(FILE)
names(depth) <- c("Chr", "locus", "depth")

depth[,`:=`(Ref = vapply(strsplit(Chr,"_"), `[`, 1, FUN.VALUE=character(1)) )]
depth2plot <- depth[, .(window.start = rollapply(locus, width=SIZE, by=STEP, FUN=min, align="left"),window.end = rollapply(locus, width=SIZE, by=STEP, FUN=max, align="left"),coverage.mean = rollapply(depth, width=SIZE, by=STEP, FUN=mean, align="left")), .(Chr)]
write.csv(depth2plot,WindowsPATH,row.names=FALSE)

depth2plot[,`:=`(Ref = vapply(strsplit(Chr,"_"), `[`, 1, FUN.VALUE=character(1)) )]
depth2plot[,`:=`(chr = vapply(strsplit(Chr,"_"), `[`, 5, FUN.VALUE=character(1)) )]
depth2plot[,`:=`(Chr = paste(Ref,chr,sep="_chr") )]

gc()

depth_contig <- depth[,.(Mean = mean(as.double(depth)),Min = min(as.double(depth)),
  Max = max(as.double(depth)),l_quantile = quantile(as.double(depth),probs=0.25),Median = median(as.double(depth)),
  u_quantile = quantile(as.double(depth),probs=0.75),Above_10X = plus_10(as.double(depth))),.(Chr)]
for(i in 1:length(depth_contig$Chr)){
    depth_contig$Chr[i] = paste(strsplit(depth_contig$Chr,"_")[[i]][1],(strsplit(depth_contig$Chr,"_")[[i]][5]),sep="_")
}
write.csv(depth_contig,PerContigPATH,row.names=FALSE)

genome_Stats <- depth[,.(Reads = sum(as.double(depth)),Mean = mean(as.double(depth)),Min = min(as.double(depth)),
  Max = max(as.double(depth)),l_quantile = quantile(as.double(depth),probs=0.25),Median = median(as.double(depth)),
  u_quantile = quantile(as.double(depth),probs=0.75),Above_10X = plus_10(as.double(depth))),.(Ref)]
write.csv(genome_Stats,PerRefPATH,row.names=FALSE)
