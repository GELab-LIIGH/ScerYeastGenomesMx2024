pcadata=read.table("mds_data.csv", sep=",", header=TRUE)

# rehacer el PCA ####
head(pcadata)

library(ggplot2)

ggplot(pcadata, aes(x=C1,y=C2,col=PhyloGroup_SACE469))+
  geom_point()
  
justplotthese=c("Tequila","MMMG","FG","SAM2","SAM1","SAM3","NortAmOak","21. Ecuadorean",
                "23. North American oak","10. French Guiana human","9. Mexican agave","Tamaulipas")

justplotthese=c("Tequila","MMMG","FG","SAM2","SAM1","SAM3","NortAmOak","21. Ecuadorean",
                "23. North American oak","10. French Guiana human","9. Mexican agave","Tamaulipas")



justplottheseStrains=c("XB075c7","YMX005635","YMX005636","YMX005638","YMX005641","YMX005648") # MixedOrigin y Wine
justplottheseStrains=c("XA124c1","XA125c5","XA126c1","XA126c5","YMX005637") # NorthAmericanOak


pcadata$Group <- ifelse(pcadata$PhyloGroup_SACE469 %in% justplotthese, 
                             pcadata$PhyloGroup_SACE469,"others")

pcadata2=pcadata
pcadata2$Group <- ifelse(pcadata$PhyloGroup_SACE469 %in% justplotthese, 
                        pcadata$PhyloGroup_SACE469, ifelse(pcadata$FID %in% justplottheseStrains, "YMXOther", "others") )


library(viridis)
my_palette <- viridis(length(justplotthese))
my_palette[length(justplotthese)+1] = "#BBBBBBBB"

my_palette = c("#ED9B52","#B05A28","#579357","#973894","#1C72B8","#EFCDE1","#212221",
               "#B05A28","#F2B7A1","#FFFFB3","#973894","#CACA30")

#pdf("MDSplot_fig1A.pdf",7.5,5)
ggplot(pcadata, aes(x=C1,y=C2,col=Group))+
  geom_point(size = ifelse(pcadata$Group != "others", 4, 1)) +
  scale_color_manual(values = my_palette) +  
  theme_minimal()
#dev.off()
  
my_palette2 = c("#ED9B52","#B05A28","#579357","#973894","#1C72B8","#EFCDE1","#212221",
               "#B05A28","#F2B7A1","#FFFFB3","#973894","#CACA30", "#CCCC30")

ggplot(pcadata2, aes(x=C1,y=C2,col=Group))+
  geom_point(size = ifelse(pcadata2$Group != "others", 4, 1)) +
  scale_color_manual(values = my_palette2) +  
  theme_minimal()


#pdf("MDSplot_fig1A_2.pdf",3.75,2.5)
ggplot(pcadata2, aes(x = C1, y = C2, fill = Group)) +
  geom_point(size = ifelse(pcadata2$Group != "others", 4, 1), 
             shape = 21, 
             color = ifelse(pcadata2$ThisStudy == 1, "black", 'gray'), 
             #color = ifelse(pcadata2$ThisStudy == 1, "black", ifelse(pcadata2$ThisStudy == 1, "pink","gray")), 
             stroke = ifelse(pcadata2$Group != "others", 1, 0)) +
  scale_fill_manual(values = my_palette2) +
  theme_minimal()
#dev.off()


justplotthese=c("Tequila","MMMG","FG","SAM2","SAM1","SAM3","NortAmOak","21. Ecuadorean",
                "23. North American oak","10. French Guiana human","9. Mexican agave","Tamaulipas")

justplotthese = unique(pcadata[!is.na(pcadata$ToColorMds),"Phylogenetic_Group"])

pcadata3=pcadata
pcadata3$Group <- ifelse(pcadata$Phylogenetic_Group %in% justplotthese, 
                         pcadata$Phylogenetic_Group, 
                         ifelse(pcadata$FID %in% justplottheseStrains, "YMXOther", "others") )



ggplot(pcadata, aes(x=C1,y=C2,col=Phylogenetic_Group))+
  geom_point(size = ifelse(pcadata$ToColorMds == 1, 4, 1)) +
  #scale_color_manual(values = my_palette) +  
  theme_minimal()

##!##!##!##!##!##!##!##!##!##!#!##!
####  con los colores correctos #####
##!##!##!##!##!##!##!##!##!##!#!##!

library(ggplot2)
pcadata=read.table("mds_data_desde250328.csv", sep=",", header=TRUE)
pcadata=read.table("mds_data_hasta250328.csv", sep=",", header=TRUE)
my_palette2 = c("#ED9B52","#973894","#1C72B8","#695333","#99CB99",
                "#346734","#F4AEE4","#AB6948","#CACA30","#9B7D3D",
                "#99CB99", #99CB00 es NAOak de peter. la otra North_America_oak son las que yo agregué a ese grupo
                "#1C72B8","#99CB99", "#878787")

my_palette2 = c("#ED9B52", #FG
                "#973894", #MA1
                "#1C72B8", #MA2
                "#695333", #MixOri
                "#346734", #NorthAmOak
                "#F4AEE4", #SAM2
                "#CACA30", #Tequila
                "#AB6948", #WB3
                "#ABCB30", #Wine14
                "#9B7D3D", #MA2
                "#99CB99", #99CB00 es NAOak de peter. la otra North_America_oak son las que yo agregué a ese grupo
                "#1C72B8",
                "#99CB99", 
                "#878787")


pdf("MDSplot_fig1A_5.pdf",5.35,2.6)
ggplot(pcadata, aes(x = C1, y = C2, fill = colorMDS)) +
  geom_point(size = ifelse(pcadata$colorMDS != "others", 2, 1.2),
             shape = 21, 
             color = ifelse(pcadata$ThisStudy == 1, "black", 'gray'), 
             #color = ifelse(pcadata2$ThisStudy == 1, "black", ifelse(pcadata2$ThisStudy == 1, "pink","gray")), 
             stroke = ifelse(pcadata$colorMDS != "others", .5, 0)) +
  scale_fill_manual(values = my_palette2) +
  theme_minimal()+
  guides(fill = guide_legend(override.aes = list(size = 2)))+
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    panel.border = element_rect(color = "black", fill = NA)
  )
dev.off()


x

# rehacer el MDS con 4 cepas menos y  ####
##!##!##!##!##!##!##!##!##!##!#!##!
####  con los colores correctos #####
##!##!##!##!##!##!##!##!##!##!#!##!
pcadata2=read.table("240604_MDS_SetR_denoised.csv", sep=",", header=TRUE)
pcadata2=read.table("250328_MDS_SetR_denoised.csv", sep=",", header=TRUE)
pcadata2=read.table("mds_SetR_denoised_data_raw.csv", sep=",", header=TRUE)
pcadata=read.table("mds_data_desde250328.csv", sep=",", header=TRUE)

head(pcadata)
head(pcadata2)
library(ggplot2)
pca=pcadata[c("FID","Phylogenetic_Group","ToColorMds","colorMDS","ThisStudy")]
#merged_df <- merge.data.frame(Strict_origins, new_Mgd, by.x = "Row_Name", by.y="Row_Names", all.x=TRUE)

pcad2=merge.data.frame(pcadata2,pca,by.x="FID",by.y="FID",all.x=TRUE)
dim(pcad2)

my_palette2 = c("#ED9B52",#FG
                "#973894",#MA1
                "#1C72B8",#MA2
                "#695333",#MixOri
                "#99CB99",#NorthAmOak
                "#F4AEE4",#SAM2
                "#CACA30",#Tequila
                "#AB1948",#WB3
                "#9B7D3D",#Wine14
                "#1C72B8",#MA2 #99CB00 es NAOak de peter. la otra North_America_oak son las que yo agregué a ese grupo
                "#AAAAAA",# Others
                "#1C12B8",#
                "#878787")#"#99CB99",

pdf("MDSplot_fig1A_6.pdf",5.55,2.8)
ggplot(pcad2, aes(x = C1, y = C2, fill = colorMDS)) +
  geom_point(size = ifelse(pcad2$colorMDS != "others", 2, 1.2),
             shape = 21, 
             color = ifelse(pcad2$ThisStudy == 1, "black", 'gray'), 
             #color = ifelse(pcadata2$ThisStudy == 1, "black", ifelse(pcadata2$ThisStudy == 1, "pink","gray")), 
             stroke = ifelse(pcad2$colorMDS != "others", .5, 0)) +
  scale_fill_manual(values = my_palette2) +
  theme_minimal()+
  guides(fill = guide_legend(override.aes = list(size = 2)))+
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    panel.border = element_rect(color = "black", fill = NA)
  )
dev.off()

