####################################################################################
##########     CH donor 5 solo ####################################
#########################################################################################
#########################################################################################



#Other tissues
#1) upload the cpgs of the skeletal tissue from the paper  and also other tissues
Skeletal.muscle.tissue.DNA.methylation <- read.delim("~/Documents/Newcastle Data/common files/different tissues DNA methylation paper/Skeletal muscle tissue DNA methylation.txt", dec=",")#1403
hPSC <- read.delim("~/Documents/Newcastle Data/common files/different tissues DNA methylation paper/hPSC.txt")#1985
Kidney <- read.delim("~/Documents/Newcastle Data/common files/different tissues DNA methylation paper/Kidney.txt", dec=",")#1007
Pancreas <- read.delim("~/Documents/Newcastle Data/common files/different tissues DNA methylation paper/Pancreas.txt", dec=",")#3333
Brain <- read.delim("~/Documents/Newcastle Data/common files/different tissues DNA methylation paper/Brain.txt", dec=",")#only a few hyper
Heart <- read.delim("~/Documents/Newcastle Data/common files/different tissues DNA methylation paper/Heart.txt", dec=",")# only a few hyper
Thymus <- read.delim("~/Documents/Newcastle Data/common files/different tissues DNA methylation paper/Thymus.txt", dec=",")
CH_donor5_solo_hypocpgs <- read.table("~/Documents/Newcastle Data/Analysis of Matt things/Fusion/CH_donor5_solo_hypocpgs.txt", header=T, quote="\"")
CH_donor5_solo_hypercpgs <- read.table("~/Documents/Newcastle Data/Analysis of Matt things/Fusion/CH_donor5_solo_hypercpgs.txt", header=T, quote="\"")



# A) Hypomethylated cpgs in CH_5_solo tissue

head(CH_donor5_solo_hypocpgs)
CH_5_solo_tissue_hypocpgs<-as.character(CH_donor5_solo_hypocpgs[,1] )
length(CH_5_solo_tissue_hypocpgs)#6307

#with the OB
intersect(CH_5_solo_tissue_hypocpgs,OB_array2_cpglist_hypo)#
sapply(intersect(CH_5_solo_tissue_hypocpgs,OB_array2_cpglist_hypo),cpg2gen)
intersect(CH_5_solo_tissue_hypocpgs,OB_array2_cpglist_hyper)#12
sapply(intersect(CH_5_solo_tissue_hypocpgs,OB_array2_cpglist_hyper),cpg2gen)

#With the OA
intersect(CH_5_solo_tissue_hypocpgs,OA_cpglist_hypo)#
sapply(intersect(CH_5_solo_tissue_hypocpgs,OA_cpglist_hypo),cpg2gen)
intersect(CH_5_solo_tissue_hypocpgs,OA_cpglist_hyper)#
sapply(intersect(CH_5_solo_tissue_hypocpgs,OA_cpglist_hyper),cpg2gen)

#with the CH
intersect(CH_5_solo_tissue_hypocpgs,CH_cpglist_hypo)#
sapply(intersect(CH_5_solo_tissue_hypocpgs,CH_cpglist_hypo),cpg2gen)
intersect(CH_5_solo_tissue_hypocpgs,CH_cpglist_hyper)#
sapply(intersect(CH_5_solo_tissue_hypocpgs,CH_cpglist_hyper),cpg2gen)
#with the HK
#23
#36
# A) Hypermethylated cpgs in CH_5_solo tissue

head(CH_donor5_solo_hypercpgs)
CH_5_solo_tissue_hypercpgs<-as.character(CH_donor5_solo_hypercpgs[,1] )
length(CH_5_solo_tissue_hypercpgs)#740

#with OB
intersect(CH_5_solo_tissue_hypercpgs,OB_array2_cpglist_hypo)#1
sapply(intersect(CH_5_solo_tissue_hypercpgs,OB_array2_cpglist_hypo),cpg2gen)
intersect(CH_5_solo_tissue_hypercpgs,OB_array2_cpglist_hyper)#0
sapply(intersect(CH_5_solo_tissue_hypercpgs,OB_array2_cpglist_hyper),cpg2gen)

#With the OA
intersect(CH_5_solo_tissue_hypercpgs,OA_cpglist_hypo)#0
sapply(intersect(CH_5_solo_tissue_hypercpgs,OA_cpglist_hypo),cpg2gen)
intersect(CH_5_solo_tissue_hypercpgs,OA_cpglist_hyper)#1
sapply(intersect(CH_5_solo_tissue_hypercpgs,OA_cpglist_hyper),cpg2gen)

#with the CH
intersect(CH_5_solo_tissue_hypercpgs,CH_cpglist_hypo)#10
sapply(intersect(CH_5_solo_tissue_hypercpgs,CH_cpglist_hypo),cpg2gen)
intersect(CH_5_solo_tissue_hypercpgs,CH_cpglist_hyper)#24
sapply(intersect(CH_5_solo_tissue_hypercpgs,CH_cpglist_hyper),cpg2gen)
#with the HK
#30
#29

#Lenghts of the overlapings
#HypoCPGs

#with the OB
length(intersect(CH_5_solo_tissue_hypocpgs,OB_array2_cpglist_hypo))#571
length(intersect(CH_5_solo_tissue_hypocpgs,OB_array2_cpglist_hyper))#12
#With the OA
length(intersect(CH_5_solo_tissue_hypocpgs,OA_cpglist_hypo))#459
length(intersect(CH_5_solo_tissue_hypocpgs,OA_cpglist_hyper))#191
#with the CH
length(intersect(CH_5_solo_tissue_hypocpgs,CH_cpglist_hypo))#2210
length(intersect(CH_5_solo_tissue_hypocpgs,CH_cpglist_hyper))#0
#with the H-K
length(intersect(CH_5_solo_tissue_hypocpgs,HvsK_cpglist_hypo))#208
length(intersect(CH_5_solo_tissue_hypocpgs,HvsK_cpglist_hyper))#50



#Hyper CPGs
length(intersect(CH_5_solo_tissue_hypercpgs,OB_array2_cpglist_hypo))#38
length(intersect(CH_5_solo_tissue_hypercpgs,OB_array2_cpglist_hyper))#23
#With the OA
length(intersect(CH_5_solo_tissue_hypercpgs,OA_cpglist_hypo))#43
length(intersect(CH_5_solo_tissue_hypercpgs,OA_cpglist_hyper))#20
#with the CH
length(intersect(CH_5_solo_tissue_hypercpgs,CH_cpglist_hypo))#0
length(intersect(CH_5_solo_tissue_hypercpgs,CH_cpglist_hyper))#106
#with the H-K
length(intersect(CH_5_solo_tissue_hypercpgs,HvsK_cpglist_hypo))#9
length(intersect(CH_5_solo_tissue_hypercpgs,HvsK_cpglist_hyper))#8



# trasncriptogram generation for  the CH_5_solo tissue and also for the clock
#########################################################################################calculation of the transcription factors associated with these cpgs
#Firts is necessary to use another script to make this chuck of code
TFnames4cpg<-function (x){
  B<-InfiniumMethylation[ which(names(InfiniumMethylation)==x) ]
  A<-TXF_Ruddy
  D<-as.data.frame(subsetByOverlaps(A,B))
  df<-as.data.frame (rownames(D))
  colnames(df)<-c("TXF")
  df<-t(df)
  return(df) 
}

Frecuency_CH_5_solo_tissue_hypocpgs<-table((unlist(lapply(CH_5_solo_tissue_hypocpgs,TFnames4cpg))))
Frecuency_CH_5_solo_tissue_hypocpgs
Frecuency_CH_5_solo_tissue_hypocpgs<-as.data.frame(Frecuency_CH_5_solo_tissue_hypocpgs)
Frecuency_CH_5_solo_tissue_hypocpgs_TXF<-Frecuency_CH_5_solo_tissue_hypocpgs
write.table(Frecuency_CH_5_solo_tissue_hypocpgs_TXF, file = "~/Documents/Newcastle Data/common files/Frecuency_CH_5_solo_tissue_hypocpgs_TXF.txt", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"))

Frecuency_CH_5_solo_tissue_hypercpgs<-table((unlist(lapply(CH_5_solo_tissue_hypercpgs,TFnames4cpg))))
Frecuency_CH_5_solo_tissue_hypercpgs
Frecuency_CH_5_solo_tissue_hypercpgs<-as.data.frame(Frecuency_CH_5_solo_tissue_hypercpgs)
Frecuency_CH_5_solo_tissue_hypercpgs_TXF<-Frecuency_CH_5_solo_tissue_hypercpgs
write.table(Frecuency_CH_5_solo_tissue_hypercpgs_TXF, file = "~/Documents/Newcastle Data/common files/Frecuency_CH_5_solo_tissue_hypercpgs_TXF.txt", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"))

length(CH_5_solo_tissue_hypocpgs)#6915
length(CH_5_solo_tissue_hypercpgs)#8324
CH_5_solo_tissue_TXF_hypohyper<-merge(Frecuency_CH_5_solo_tissue_hypocpgs_TXF, Frecuency_CH_5_solo_tissue_hypercpgs_TXF, by.x="Var1",by.y="Var1",all.x=TRUE,all.y=TRUE)
CH_5_solo_tissue_TXF_hypohyper
#and now we transform the NA into 0s
CH_5_solo_tissue_TXF_hypohyper_ceros<-CH_5_solo_tissue_TXF_hypohyper
CH_5_solo_tissue_TXF_hypohyper_ceros[is.na(CH_5_solo_tissue_TXF_hypohyper_ceros)] <- 0
CH_5_solo_tissue_TXF_hypohyper_ceros
colnames(CH_5_solo_tissue_TXF_hypohyper_ceros)<-c("TXF","CH_5_solo_hypomet","CH_5_solo_hypermet")
colnames(CH_5_solo_tissue_TXF_hypohyper_ceros)
CH_5_solo_tissue_TXF_hypohyper_ceros

plot(CH_5_solo_tissue_TXF_hypohyper_ceros[1:50,2], type="l",col="red")
lines(CH_5_solo_tissue_TXF_hypohyper_ceros[1:50,3], type="l",col=19)
#I use the work previously done with all the cell fates
Allfatesandtissues_TXF_hypohyper_ceros
Allfatesandtissues_TXF_hypohyper<-merge(Allfatesandtissues_TXF_hypohyper_ceros,CH_5_solo_tissue_TXF_hypohyper_ceros, by.x="TXF.genome.cpgs",by.y="TXF",all.x=TRUE,all.y=TRUE)
Allfatesandtissues_TXF_hypohyper_ceros<-Allfatesandtissues_TXF_hypohyper
Allfatesandtissues_TXF_hypohyper_ceros[is.na(Allfatesandtissues_TXF_hypohyper_ceros)] <- 0
Allfatesandtissues_TXF_hypohyper_ceros
Allfatesandtissues_TXF_hypohyper_ceros[,46]<-(Allfatesandtissues_TXF_hypohyper_ceros[,44]*100)/6307
Allfatesandtissues_TXF_hypohyper_ceros[,47]<-(Allfatesandtissues_TXF_hypohyper_ceros[,45]*100)/740
Allfatesandtissues_TXF_hypohyper_ceros

write.table(Allfatesandtissues_TXF_hypohyper_ceros, file = "~/Documents/Newcastle Data/common files/Frequency CH_OB_HIPOA/Allfatesandtissues_TXF_hypohyper_ceros_CH_5_solo.txt", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"))





plot.new()
barplot(Allfatesandtissues_TXF_hypohyper_ceros[,48],type="n", ylim=c(-10,10), col="red", ,xlab="TXFs", col.lab="black",ylab="Percentage %", col.lab="black")
axis(1, at=1:161, lab=Allfatesandtissues_TXF_hypohyper_ceros[,1])

lines(Allfatesandtissues_TXF_hypohyper_ceros[,48],  col="red")
lines(Allfatesandtissues_TXF_hypohyper_ceros[,49],  col="green4")
lines(Allfatesandtissues_TXF_hypohyper_ceros[,50],  col="yellow")
lines(Allfatesandtissues_TXF_hypohyper_ceros[,51],  col="orange")
lines(Allfatesandtissues_TXF_hypohyper_ceros[,52],  col="purple")
lines(Allfatesandtissues_TXF_hypohyper_ceros[,53],  col="brown")
lines(Allfatesandtissues_TXF_hypohyper_ceros[,54],  col="grey")
lines(Allfatesandtissues_TXF_hypohyper_ceros[,55],  col="yellow4")
lines(Allfatesandtissues_TXF_hypohyper_ceros[,56],  col="orange3")


Allfatesandtissues_TXF_hypohyper_ceros[,48]
Allfatesandtissues_TXF_hypohyper_ceros[,49]
Allfatesandtissues_TXF_hypohyper_ceros[,50]
Allfatesandtissues_TXF_hypohyper_ceros[,51]
Allfatesandtissues_TXF_hypohyper_ceros[,52]
Allfatesandtissues_TXF_hypohyper_ceros[,53]
Allfatesandtissues_TXF_hypohyper_ceros[,54]
Allfatesandtissues_TXF_hypohyper_ceros[,55]
Allfatesandtissues_TXF_hypohyper_ceros[,56]



Allfatesandtissues_TXF_hypohyper_ceros[,48]<-log2(Allfatesandtissues_TXF_hypohyper_ceros[,13]/Allfatesandtissues_TXF_hypohyper_ceros[,14])
Allfatesandtissues_TXF_hypohyper_ceros[,49]<-log2(Allfatesandtissues_TXF_hypohyper_ceros[,18]/Allfatesandtissues_TXF_hypohyper_ceros[,19])
Allfatesandtissues_TXF_hypohyper_ceros[,50]<-log2(Allfatesandtissues_TXF_hypohyper_ceros[,22]/Allfatesandtissues_TXF_hypohyper_ceros[,23])
Allfatesandtissues_TXF_hypohyper_ceros[,51]<-log2(Allfatesandtissues_TXF_hypohyper_ceros[,26]/Allfatesandtissues_TXF_hypohyper_ceros[,27])
Allfatesandtissues_TXF_hypohyper_ceros[,52]<-log2(Allfatesandtissues_TXF_hypohyper_ceros[,34]/Allfatesandtissues_TXF_hypohyper_ceros[,35])
Allfatesandtissues_TXF_hypohyper_ceros[,53]<-log2(Allfatesandtissues_TXF_hypohyper_ceros[,38]/Allfatesandtissues_TXF_hypohyper_ceros[,39])
Allfatesandtissues_TXF_hypohyper_ceros[,54]<-log2(Allfatesandtissues_TXF_hypohyper_ceros[,9]/Allfatesandtissues_TXF_hypohyper_ceros[,10])
Allfatesandtissues_TXF_hypohyper_ceros[,55]<-log2(Allfatesandtissues_TXF_hypohyper_ceros[,8]/Allfatesandtissues_TXF_hypohyper_ceros[,7])
Allfatesandtissues_TXF_hypohyper_ceros[,56]<-log2(Allfatesandtissues_TXF_hypohyper_ceros[,46]/Allfatesandtissues_TXF_hypohyper_ceros[,47])





Allfatesandtissues_TXF_hypohyper_ceros[,57]<-log2(Allfatesandtissues_TXF_hypohyper_ceros[,46]/Allfatesandtissues_TXF_hypohyper_ceros[,47])
Allfatesandtissues_TXF_hypohyper_ceros[,58]<-log2(Allfatesandtissues_TXF_hypohyper_ceros[,46]/Allfatesandtissues_TXF_hypohyper_ceros[,47])








plot(ratioCH_donor5_solo,type="o", ylim=c(-5,5), col="red", ,xlab="TXFs", col.lab="black",ylab="Percentage %", col.lab="black")
axis(1, at=1:161, lab=Allfatesandtissues_TXF_hypohyper_ceros[,1])
lines(Allfatesandtissues_TXF_hypohyper_ceros[,46],  col="red")
lines(Allfatesandtissues_TXF_hypohyper_ceros[,47],  col="green4")
lines(ratioCH_donor5_solo,  col="black")

#Bar plot for individual Transcription factors
TXF2plot<-"SUZ12"
numberofTXF<-match(TXF2plot,Allfatesandtissues_TXF_hypohyper_ceros[,1])
numberofTXF
barplot(as.numeric(Allfatesandtissues_TXF_hypohyper_ceros[numberofTXF,c(2,13,14,18,19,22,23,26,27,34,35,38,39,9,10, 8,7,46,47)]),
        names.arg=c("Array","OB","OB","Adipose","Adipose","Muscle","Muscle","hPSC","hPSC","Pancreas","Pancreas","Thymus", "Thymus","OA", "OA", "CH","CH","CH5","CH5"),
        col=c("red", "white","black","white","black","white","black","white","black","white","black","white","black","white","blue","white","green","white","purple"), 
        main=TXF2plot, ylab="Percentage %" )
