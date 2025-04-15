source("/bioware/bin/SingleCell/Func.R")
setwd("/mnt/date/Project/yupeng/iCell8/Zongna_mouse_TAC/clean_data")
load("TAC_sc.RData")
#the data should have a umi matrix file and a cell info file 

mus_ensembl_gene<-read.table("/mnt/date/genomelib/annotation/Mus_musculus.GRCm38.75.info")
cell_umi_clean<-cell_umi_combine[rownames(cell_umi_combine) %in% mus_ensembl_gene[,2],]
rownames(cell_umi_clean)<-mus_ensembl_gene[match(rownames(cell_umi_clean),mus_ensembl_gene[,2]),1]
TAC_new_umi<-cell_umi_clean

TAC_new_umi_nomit<- TAC_new_umi[-grep("^mt-",rownames(TAC_new_umi)),]
TAC_preTsne<-preTSNE(TAC_new_umi_nomit,TAC_cell_info,"TAC_180411",regress_factor="nUMI")
TAC_tsne<-doTSNE_2(TAC_preTsne,"TAC_180411",1:10,1.5)

TAC_tsne$seurat@meta.data$sample<- TAC_cell_info[rownames(TAC_tsne$seurat@meta.data),"animal"]
TAC_tsne$seurat@meta.data$group<- sub("CM|NCM","",TAC_tsne$seurat@meta.data$group)
TAC_tsne$seurat@meta.data$group[TAC_tsne$seurat@meta.data$group ==""] <-"PI"
write_marker_smp(TAC_tsne,"TAC_tsne_marker")

outlier_cell<-rownames(TAC_tsne$seurat@meta.data)[TAC_tsne$seurat@meta.data$res.1.5 %in% c(6,11,16,20)]
TAC_clean_umi_nomit <- TAC_new_umi_nomit[,-which(colnames(TAC_new_umi_nomit) %in% outlier_cell)]
 
TAC_clean_preTsne<-preTSNE(TAC_clean_umi_nomit,TAC_cell_info,"TAC_clean",regress_factor="nUMI")
TAC_clean_tsne<-doTSNE_2(TAC_clean_preTsne,"TAC_clean",1:10,1.5)

TAC_clean_tsne$seurat@meta.data$sample<- TAC_cell_info[rownames(TAC_clean_tsne$seurat@meta.data),"animal"]
TAC_clean_tsne$seurat@meta.data$group<- sub("CM|NCM","",TAC_clean_tsne$seurat@meta.data$group)
TAC_clean_tsne$seurat@meta.data$group[TAC_clean_tsne$seurat@meta.data$group ==""] <-"PI"
write_marker(TAC_clean_tsne,"TAC_clean_res15",marker=c("Actb","Gapdh","Eif3d","Eif5","Rpl5","Rpl8"))
save(TAC_clean_preTsne,TAC_clean_tsne,TAC_clean_umi_nomit,TAC_cell_info,file="TAC_clean_20180411.RData")

#
CM：0,2,3,5,9,11,13,15,16
EC：4,8,14
FB：6,7,17,20,22
MP：1,12,19
Granulocytes：18
T：10,21

CM_cell<-rownames(subset(TAC_clean_tsne$seurat@meta.data, res.1.5 %in% c(0,2,3,5,9,11,13,15,16)))
EC_cell<-rownames(subset(TAC_clean_tsne$seurat@meta.data, res.1.5 %in% c(4,8,14)))
FB_cell<- rownames(subset(TAC_clean_tsne$seurat@meta.data, res.1.5 %in% c(6,7,17,20,22)))
MP_cell<- rownames(subset(TAC_clean_tsne$seurat@meta.data, res.1.5 %in% c(1,12,19)))
Granulo_cell<- rownames(subset(TAC_clean_tsne$seurat@meta.data, res.1.5 %in% c(18)))
T_cell <- rownames(subset(TAC_clean_tsne$seurat@meta.data, res.1.5 %in% c(10,21)))

CM_umi<- TAC_clean_tsne$seurat@raw.data[,CM_cell]
EC_umi<- TAC_clean_tsne$seurat@raw.data[,EC_cell]
FB_umi<- TAC_clean_tsne$seurat@raw.data[,FB_cell]
MP_umi<- TAC_clean_tsne$seurat@raw.data[,MP_cell]
Granulo_umi<- TAC_clean_tsne$seurat@raw.data[,Granulo_cell]
T_umi<- TAC_clean_tsne$seurat@raw.data[,T_cell]


CM_monocle<- do_monocle(CM_umi,TAC_cell_info)
EC_monocle<- do_monocle(EC_umi,TAC_cell_info)
FB_monocle<- do_monocle(FB_umi,TAC_cell_info)
MP_monocle<- do_monocle(MP_umi,TAC_cell_info)
Granulo_monocle<- do_monocle(Granulo_umi,TAC_cell_info)
T_monocle<- do_monocle(T_umi,TAC_cell_info)


CM_preTSNE<- preTSNE(CM_umi,TAC_cell_info,"CM",regress_factor="nUMI")
#CM_TSNE<- doTSNE_2(CM_preTSNE,"CM",1:8,1)

CM_TSNE_p8<-CM_TSNE
CM_TSNE<- doTSNE_2(CM_preTSNE,"CM",1:10,1)

EC_preTSNE<- preTSNE(EC_umi,TAC_cell_info,"EC",regress_factor="nUMI")
EC_TSNE<- doTSNE_2(EC_preTSNE,"EC",1:10,1)

FB_preTSNE<- preTSNE(FB_umi,TAC_cell_info,"FB",regress_factor="nUMI")
FB_TSNE<- doTSNE_2(FB_preTSNE,"FB",1:10,1)

MP_preTSNE<- preTSNE(MP_umi,TAC_cell_info,"MP",regress_factor="nUMI")
MP_TSNE<- doTSNE_2(MP_preTSNE,"MP",1:10,1)

Granulo_preTSNE<- preTSNE(Granulo_umi,TAC_cell_info,"Granulo",regress_factor="nUMI")
Granulo_TSNE<- doTSNE_2(Granulo_preTSNE,"Granulo",1:10,1)

T_preTSNE<- preTSNE(T_umi,TAC_cell_info,"T",regress_factor="nUMI")
T_TSNE<- doTSNE_2(T_preTSNE,"T",1:10,1)

CM_TSNE$seurat@meta.data$sample<- TAC_cell_info[rownames(CM_TSNE$seurat@meta.data),"animal"]
CM_TSNE$seurat@meta.data$group<- sub("CM|NCM","",CM_TSNE$seurat@meta.data$group)
CM_TSNE$seurat@meta.data$group[CM_TSNE$seurat@meta.data$group ==""] <-"PI"

EC_TSNE$seurat@meta.data$sample<- TAC_cell_info[rownames(EC_TSNE$seurat@meta.data),"animal"]
EC_TSNE$seurat@meta.data$group<- sub("CM|NCM","",EC_TSNE$seurat@meta.data$group)
EC_TSNE$seurat@meta.data$group[EC_TSNE$seurat@meta.data$group ==""] <-"PI"

FB_TSNE$seurat@meta.data$sample<- TAC_cell_info[rownames(FB_TSNE$seurat@meta.data),"animal"]
FB_TSNE$seurat@meta.data$group<- sub("CM|NCM","",FB_TSNE$seurat@meta.data$group)
FB_TSNE$seurat@meta.data$group[FB_TSNE$seurat@meta.data$group ==""] <-"PI"

MP_TSNE$seurat@meta.data$sample<- TAC_cell_info[rownames(MP_TSNE$seurat@meta.data),"animal"]
MP_TSNE$seurat@meta.data$group<- sub("CM|NCM","",MP_TSNE$seurat@meta.data$group)
MP_TSNE$seurat@meta.data$group[MP_TSNE$seurat@meta.data$group ==""] <-"PI"

Granulo_TSNE$seurat@meta.data$sample<- TAC_cell_info[rownames(Granulo_TSNE$seurat@meta.data),"animal"]
Granulo_TSNE$seurat@meta.data$group<- sub("CM|NCM","",Granulo_TSNE$seurat@meta.data$group)
Granulo_TSNE$seurat@meta.data$group[Granulo_TSNE$seurat@meta.data$group ==""] <-"PI"

T_TSNE$seurat@meta.data$sample<- TAC_cell_info[rownames(T_TSNE$seurat@meta.data),"animal"]
T_TSNE$seurat@meta.data$group<- sub("CM|NCM","",T_TSNE$seurat@meta.data$group)
T_TSNE$seurat@meta.data$group[T_TSNE$seurat@meta.data$group ==""] <-"PI"

CM_TSNE$seurat@meta.data$percent.mito<- TAC_cell_info[rownames(CM_TSNE$seurat@meta.data),"mito.perc"]
EC_TSNE$seurat@meta.data$percent.mito<- TAC_cell_info[rownames(EC_TSNE$seurat@meta.data),"mito.perc"]
FB_TSNE$seurat@meta.data$percent.mito<- TAC_cell_info[rownames(FB_TSNE$seurat@meta.data),"mito.perc"]
MP_TSNE$seurat@meta.data$percent.mito<- TAC_cell_info[rownames(MP_TSNE$seurat@meta.data),"mito.perc"]
Granulo_TSNE$seurat@meta.data$percent.mito<- TAC_cell_info[rownames(Granulo_TSNE$seurat@meta.data),"mito.perc"]
T_TSNE$seurat@meta.data$percent.mito<- TAC_cell_info[rownames(T_TSNE$seurat@meta.data),"mito.perc"]


write_marker_smp(CM_TSNE,"CM_tsne_marker",class_level=c("SHAM","TAC2W","TAC5W","TAC8W","TAC11W"),axis.tick=c("0","2","5","8","11"))
write_marker_smp(EC_TSNE,"EC_tsne_marker",class_level=c("SHAM","TAC2W","TAC5W","TAC8W","TAC11W"),axis.tick=c(0,2,5,8,11))
write_marker_smp(FB_TSNE,"FB_tsne_marker",class_level=c("SHAM","TAC2W","TAC5W","TAC8W","TAC11W"),axis.tick=c(0,2,5,8,11))
write_marker_smp(MP_TSNE,"MP_tsne_marker",class_level=c("SHAM","TAC2W","TAC5W","TAC8W","TAC11W"),axis.tick=c(0,2,5,8,11))
write_marker_smp(Granulo_TSNE,"Granulo_tsne_marker",class_level=c("SHAM","TAC2W","TAC5W","TAC8W","TAC11W"),axis.tick=c(0,2,5,8,11))
write_marker_smp(T_TSNE,"T_tsne_marker",class_level=c("SHAM","TAC2W","TAC5W","TAC8W","TAC11W"),axis.tick=c(0,2,5,8,11))
save(CM_TSNE,EC_TSNE,FB_TSNE,MP_TSNE,Granulo_TSNE, T_TSNE, file="TAC_subCell_tsne_0413.RData")

create_cluster_sum <- function (x,column,cell_type){
	
	data_df<-data.frame(cluster=x$seurat@ident,class=x$seurat@meta.data[,as.character(column)],value=1)
	
	cluster_sum<-aggregate(value~cluster, data_df,sum)
	group_sum<-aggregate(value~cluster+class, data_df,sum)
	group_percent<-data.frame(group_sum,percent=0)
	class_percent<-data.frame(group_sum,percent=0)

	total_class_sum<-aggregate(value~class, data_df,sum)
	total_sum<-sum(data_df$value)

	cluster_total<-max(as.numeric(data_df$cluster))

	total_class_df<-data.frame(cluster=factor("total"),total_class_sum,percent=round(total_class_sum[,"value"]/total_sum,2)*100)

	for (i in 1:nrow(group_percent)){

			cluster0=group_percent[i,"cluster"]
			class0=group_percent[i,"class"]

			group_percent[i,"percent"]=round(group_percent[i,"value"]/cluster_sum[cluster_sum$cluster==cluster0,"value"]*100,2)
			class_percent[i,"percent"]=round(class_percent[i,"value"]/total_class_sum[total_class_sum$class==class0,"value"]*100,2)

	}

	#fill lines of cluster/class combination with values of 0

	for (clust in names(table(class_percent$cluster))){

		for (cond in names(table(class_percent$class))){

			if(nrow(subset(class_percent,cluster == clust & class == cond))==0){

				class_percent<-rbind(class_percent,c(clust,cond,0,0))
			}
		}

	}
	
	class_percent$percent<-as.numeric(class_percent$percent)
    class_percent$class<-factor(class_percent$class)
	class_percent$cell_type<-cell_type
	return(class_percent)

}

CM_cluster_sum<-create_cluster_sum(CM_TSNE,"condition","CM")
EC_cluster_sum<-create_cluster_sum(EC_TSNE,"condition","EC")
FB_cluster_sum<-create_cluster_sum(FB_TSNE,"condition","FB")
MP_cluster_sum<-create_cluster_sum(MP_TSNE,"condition","MP")
Granulo_cluster_sum<-create_cluster_sum(Granulo_TSNE,"condition","Granulo")
T_cluster_sum<-create_cluster_sum(T_TSNE,"condition","T")

All_cluster_sum<-rbind(CM_cluster_sum,EC_cluster_sum,FB_cluster_sum,MP_cluster_sum,Granulo_cluster_sum,T_cluster_sum)

All_cluster_sum$cluster_id<-paste0(All_cluster_sum$cell_type,All_cluster_sum$cluster)

All_cluster_sum<-All_cluster_sum[order(All_cluster_sum[,6],All_cluster_sum[,2]),]

All_cluster_time_df<-data.frame(SHAM=subset(All_cluster_sum,class=="SHAM",select="percent"),TAC2W=subset(All_cluster_sum,class=="TAC2W",select="percent"),TAC5W=subset(All_cluster_sum,class=="TAC5W",select="percent"),TAC8W=subset(All_cluster_sum,class=="TAC8W",select="percent"),TAC11W=subset(All_cluster_sum,class=="TAC11W",select="percent"),row.names=subset(All_cluster_sum,class=="SHAM",select="cluster_id")[,1])

colnames(All_cluster_time_df)<-c("SHAM","TAC2W","TAC5W","TAC8W","TAC11W")
#All_cluster_time_df_select<- All_cluster_time_df[apply(All_cluster_time_df,1,max)>8,]

All_cluster_time_df[All_cluster_time_df==0]<-0.01

All_cluster_time_df_ratio<- data.frame(s1=All_cluster_time_df[,2]/All_cluster_time_df[,1],s2=All_cluster_time_df[,3]/All_cluster_time_df[,2],s3=All_cluster_time_df[,4]/All_cluster_time_df[,3],s4=All_cluster_time_df[,5]/All_cluster_time_df[,4])

All_cluster_time_df_ratio<-log2(All_cluster_time_df_ratio)
rownames(All_cluster_time_df_ratio)<- rownames(All_cluster_time_df)

ncm_cluster_time_df<- All_cluster_time_df[regexec("CM",rownames(All_cluster_time_df))==-1,]
ncm_cluster_time_df_ratio <- All_cluster_time_df_ratio[regexec("CM",rownames(All_cluster_time_df_ratio))==-1,]

TAC2W_ratio_df<- data.frame(ncm_cluster_time_df[,1:2],ratio=ncm_cluster_time_df_ratio[,1])
TAC5W_ratio_df<- data.frame(ncm_cluster_time_df[,2:3],ratio=ncm_cluster_time_df_ratio[,2])
TAC8W_ratio_df<- data.frame(ncm_cluster_time_df[,3:4],ratio=ncm_cluster_time_df_ratio[,3])
TAC11W_ratio_df<- data.frame(ncm_cluster_time_df[,4:5],ratio=ncm_cluster_time_df_ratio[,4])

save(All_cluster_time_df,All_cluster_time_df_ratio, TAC2W_ratio_select,TAC5W_ratio_select,TAC8W_ratio_select,TAC11W_ratio_select,file="TAC_ratio_df.RData")

select_change_group <- function(x){

	select<- intersect(which(abs(x[,3])>median(abs(x[,3]))),which(apply(x[,1:2],1,max)>10))
	
	data<-x[select,]
	
	data<-data[order(data[,3],decreasing=T),]

	return(data)
}

TAC2W_ratio_select<- select_change_group(TAC2W_ratio_df)
TAC5W_ratio_select<- select_change_group(TAC5W_ratio_df)
TAC8W_ratio_select<- select_change_group(TAC8W_ratio_df)
TAC11W_ratio_select<- select_change_group(TAC11W_ratio_df)

CM_marker<-CM_TSNE$markers
EC_marker<-EC_TSNE$markers
FB_marker<-FB_TSNE$markers
MP_marker<-MP_TSNE$markers
Granulo_marker<-Granulo_TSNE$markers
T_marker<-T_TSNE$markers

CM_marker$cluster<-paste0("CM",CM_marker$cluster)
EC_marker$cluster<-paste0("EC",EC_marker$cluster)
FB_marker$cluster<-paste0("FB",FB_marker$cluster)
MP_marker$cluster<-paste0("MP",MP_marker$cluster)
T_marker$cluster<-paste0("T",T_marker$cluster)
Granulo_marker$cluster<-paste0("GN",Granulo_marker$cluster)

all_marker<-rbind(CM_marker,EC_marker,FB_marker,MP_marker,T_marker,Granulo_marker)
all_marker_select<-subset(all_marker,p_val_adj<0.01)

##### download STRING-db data 
#10090.protein.links.full.v10.5.txt.gz
cd /mnt/date/genomelib/IntAct/STRING-db
awk '($10>0 || $11>0 )' 10090.protein.links.full.v10.5.txt >mus_string_valid_ppi

#cd /mnt/date/genomelib/annotation
#perl extract_protein_id.pl gencode/gencode.vM15.annotation.gtf > gencode.v15.Mm.protein.info

perl 1.pl gencode.v15.Mm.protein.info mus_string_valid_ppi >mus_string_valid_ppi.symbol

goa:("transcription factor") locations:(location:"Nucleus [SL-0191]") AND reviewed:yes AND organism:"Mus musculus (Mouse) [10090]"
1490

#select membrane proteins

annotation:(type:transmem) goa:("plasma membrane [5886]") AND reviewed:yes AND organism:"Mus musculus (Mouse) [10090]"
1535

locations:(location:"Secreted [SL-0243]") AND reviewed:yes AND organism:"Mus musculus (Mouse) [10090]"
1614

awk -F "\t" '{print $5}' uniprot-Secreted |awk '($1){print $1}' |grep -v Gene >uniprot_mouse_sec
awk -F "\t" '{print $5}' uniprot-type_transmem |awk '($1){print $1}' |grep -v Gene >uniprot_mouse_mem
awk -F "\t" '{print $5}' uniprot_goa_tf  |awk '($1){print $1}' |grep -v Gene >uniprot_mouse_tf

sec<-read.table("uniprot_mouse_sec")
mem<-read.table("uniprot_mouse_mem")
sec_uniq<-as.matrix(sec[!sec[,1] %in% mem[,1],1])
write.table(sec_uniq,file="uniprot_mouse_sec",quote=F,row.names=F,col.names=F)

perl uniq_intact.pl mus_string_valid_ppi.symbol >mus_string_valid_ppi.symbol.uniq &
perl uniq_intact.pl string_valid_ppi.symbol > string_valid_ppi.symbol.uniq &

mus_interaction_data<-read.delim("mus_string_valid_ppi.symbol.uniq",sep="\t",header=F)

mus_sec<-read.table("uniprot_mouse_sec")
mus_mem<-read.table("uniprot_mouse_mem")

mus_sec<-as.matrix(unique(mus_sec[,1]))
mus_mem<-as.matrix(unique(mus_mem[,1]))
all_gene<- rbind(mus_sec,mus_mem)

mus_interact_pair <- mus_interaction_data[mus_interaction_data[,1] %in% all_gene[,1] & mus_interaction_data[,2] %in% all_gene[,1],]
write.table(mus_sec,file="uniprot_mouse_sec",row.names=F,col.names=F,quote=F)
write.table(mus_mem,file="uniprot_mouse_mem",row.names=F,col.names=F,quote=F)
write.table(mus_interact_pair,file="mouse_sec_mem_interact_pair",row.names=F,col.names=F,quote=F)

hum_interaction_data<-read.delim("string_valid_ppi.symbol.uniq",sep="\t",header=F)

hum_sec<-read.table("human_sec.txt")
hum_mem<-read.table("human_mem.txt")

hum_sec<-as.matrix(unique(hum_sec[,1]))
hum_mem<-as.matrix(unique(hum_mem[,1]))
hum_all_gene<- rbind(hum_sec,hum_mem)

hum_interact_pair <- hum_interaction_data[hum_interaction_data[,1] %in% hum_all_gene[,1] & hum_interaction_data[,2] %in% hum_all_gene[,1],]
write.table(hum_sec,file="human_sec.txt",row.names=F,col.names=F,quote=F)
write.table(hum_mem,file="human_mem.txt",row.names=F,col.names=F,quote=F)
write.table(hum_interact_pair,file="hum_sec_mem_interact_pair",row.names=F,col.names=F,quote=F)


mus_interact_pair<-read.table("mouse_interact_pair",header=F)

LA_cm_sp_marker<-subset(N_CM_LA_sp_marker,p_val_adj<0.05)
LV_cm_sp_marker<-subset(N_CM_LV_sp_marker,p_val_adj<0.05)
immCM_sp_marker<-subset(N_CM_C8_sp_marker,p_val_adj<0.05)

EC_sp_marker<-subset(ec_marker,p_val_adj<0.05)
Fibroblast_sp_marker<-subset(fibro_marker,p_val_adj<0.05)
Macrophage_sp_marker<-subset(macrophage_marker,p_val_adj<0.05)
VSMC_sp_marker<-subset(vsmc_marker,p_val_adj<0.05)
tcell_sp_marker<-subset(t_cell_marker,p_val_adj<0.05)

LA_cm_ligand<- LA_cm_sp_marker[LA_cm_sp_marker$gene %in% hum_sec[,1],"gene"]
LA_cm_receptor<- LA_cm_sp_marker[LA_cm_sp_marker$gene %in% hum_mem[,1],"gene"]

LA_cm_pair<-data.frame(cell="LA_CM",gene=c(LA_cm_ligand,LA_cm_receptor),x=c(1:length(c(LA_cm_ligand,LA_cm_receptor))),type=c(rep("ligand",length(LA_cm_ligand)),rep("receptor",length(LA_cm_receptor))))


LV_cm_ligand<- LV_cm_sp_marker[LV_cm_sp_marker$gene %in% hum_sec[,1],"gene"]
LV_cm_receptor<- LV_cm_sp_marker[LV_cm_sp_marker$gene %in% hum_mem[,1],"gene"]

LV_cm_pair<-data.frame(cell="LV_CM",gene=c(LV_cm_ligand,LV_cm_receptor),x=c(1:length(c(LV_cm_ligand,LV_cm_receptor))),type=c(rep("ligand",length(LV_cm_ligand)),rep("receptor",length(LV_cm_receptor))))


immCM_ligand<- immCM_sp_marker[immCM_sp_marker$gene %in% hum_sec[,1],"gene"]
immCM_receptor<- immCM_sp_marker[immCM_sp_marker$gene %in% hum_mem[,1],"gene"]

immCM_pair<-data.frame(cell="immCM",gene=c(immCM_ligand,immCM_receptor),x=c(1:length(c(immCM_ligand,immCM_receptor))),type=c(rep("ligand",length(immCM_ligand)),rep("receptor",length(immCM_receptor))))

EC_ligand<- EC_sp_marker[EC_sp_marker$gene %in% hum_sec[,1],"gene"]
EC_receptor<- EC_sp_marker[EC_sp_marker$gene %in% hum_mem[,1],"gene"]

EC_pair<-data.frame(cell="EC",gene=c(EC_ligand,EC_receptor),x=c(1:length(c(EC_ligand,EC_receptor))),type=c(rep("ligand",length(EC_ligand)),rep("receptor",length(EC_receptor))))

Fibroblast_ligand<- Fibroblast_sp_marker[Fibroblast_sp_marker$gene %in% hum_sec[,1],"gene"]
Fibroblast_receptor<- Fibroblast_sp_marker[Fibroblast_sp_marker$gene %in% hum_mem[,1],"gene"]

Fibroblast_pair<-data.frame(cell="Fibroblast",gene=c(Fibroblast_ligand,Fibroblast_receptor),x=c(1:length(c(Fibroblast_ligand,Fibroblast_receptor))),type=c(rep("ligand",length(Fibroblast_ligand)),rep("receptor",length(Fibroblast_receptor))))

VSMC_ligand<- VSMC_sp_marker[VSMC_sp_marker$gene %in% hum_sec[,1],"gene"]
VSMC_receptor<- VSMC_sp_marker[VSMC_sp_marker$gene %in% hum_mem[,1],"gene"]

VSMC_pair<-data.frame(cell="VSMC",gene=c(VSMC_ligand,VSMC_receptor),x=c(1:length(c(VSMC_ligand,VSMC_receptor))),type=c(rep("ligand",length(VSMC_ligand)),rep("receptor",length(VSMC_receptor))))


Macrophage_ligand<- Macrophage_sp_marker[Macrophage_sp_marker$gene %in% hum_sec[,1],"gene"]
Macrophage_receptor<- Macrophage_sp_marker[Macrophage_sp_marker$gene %in% hum_mem[,1],"gene"]

Macrophage_pair<-data.frame(cell="Macrophage",gene=c(Macrophage_ligand,Macrophage_receptor),x=c(1:length(c(Macrophage_ligand,Macrophage_receptor))),type=c(rep("ligand",length(Macrophage_ligand)),rep("receptor",length(Macrophage_receptor))))

tcell_ligand<- tcell_sp_marker[tcell_sp_marker$gene %in% hum_sec[,1],"gene"]
tcell_receptor<- tcell_sp_marker[tcell_sp_marker$gene %in% hum_mem[,1],"gene"]

tcell_pair<-data.frame(cell="T cell",gene=c(tcell_ligand,tcell_receptor),x=c(1:length(c(tcell_ligand,tcell_receptor))),type=c(rep("ligand",length(tcell_ligand)),rep("receptor",length(tcell_receptor))))


all_gene_pair<-data.frame(rbind(LA_cm_pair,LV_cm_pair,immCM_pair,EC_pair,Fibroblast_pair,VSMC_pair,Macrophage_pair,tcell_pair))

all_gene_pair$cell<-factor(all_gene_pair$cell,levels=c("LA_CM","LV_CM","immCM","EC","Fibroblast","VSMC","Macrophage","T cell"))

all_gene_pair$color<-"red"
all_gene_pair[all_gene_pair$type=="receptor","color"] <-"green"

pair_combine<-function(x,y){

	total_col<- ncol(x)+ncol(y)
	
	data<-data.frame()
	
	for(i in 1:nrow(x)){
	
		for(j in 1:nrow(y)){
		
			vec<-cbind(x[i,],y[j,])
			data<-rbind(data,vec)
		}
	
	}
	return(data)
}


all_cell_interaction<-data.frame()
for(i in 1:nrow(ligand_receptor_pair)){

	l<-as.character(ligand_receptor_pair[i,1])
	r<-as.character(ligand_receptor_pair[i,2])
	
	l_matrix<-all_gene_pair[all_gene_pair$gene ==l,]
	r_matrix<-all_gene_pair[all_gene_pair$gene ==r,]
	
	if(nrow(l_matrix)>0 & nrow(r_matrix)>0){
		l_r_matrix<-pair_combine(l_matrix,r_matrix)
		all_cell_interaction<-rbind(all_cell_interaction,l_r_matrix)
	}
}



pData(CM_monocle)$cluster <- CM_TSNE$seurat@meta.data[rownames(pData(CM_monocle)),"res.1"]
pData(EC_monocle)$cluster <- EC_TSNE$seurat@meta.data[rownames(pData(EC_monocle)),"res.1"]
pData(FB_monocle)$cluster <- FB_TSNE$seurat@meta.data[rownames(pData(FB_monocle)),"res.1"]
pData(MP_monocle)$cluster <- MP_TSNE$seurat@meta.data[rownames(pData(MP_monocle)),"res.1"]
pData(Granulo_monocle)$cluster <- Granulo_TSNE$seurat@meta.data[rownames(pData(Granulo_monocle)),"res.1"]
pData(T_monocle)$cluster <- T_TSNE$seurat@meta.data[rownames(pData(T_monocle)),"res.1"]


calculate_state_sum <-function(x,title){
	
	x$value=1
	state_sum<-aggregate(value ~ State+condition,x,sum)
	ggplot(state_sum,aes(x=State,y=value,fill=condition))+geom_bar(stat='identity')+ggtitle(title)
	
}

calculate_tsne_percent <-function(x,title){
	
	x$value=1
	cluster_sum<-aggregate(value ~ State+cluster,x,sum)
	state_sum<-aggregate(value ~ State,x,sum)
	cluster_percent<-data.frame(cluster_sum,percent=0)
	
	 for (i in 1:nrow(cluster_percent)){

		cluster0=cluster_percent[i,"cluster"]
		state0=cluster_percent[i,"State"]

		cluster_percent[i,"percent"]=round(cluster_percent[i,"value"]/state_sum[state_sum$State==state0,"value"]*100,2)
		
     }

	ggplot(cluster_percent,aes(x=State,y=percent,fill=cluster))+geom_bar(stat='identity')+ggtitle(title)
	
}


pdf("TAC_monocle_latest_complex.pdf")

#pData(CM_monocle)$condition<- factor(pData(CM_monocle)$condition,levels=c("SHAM","TAC2W","TAC5W","TAC8W","TAC11W"))
#calculate_state_sum(pData(CM_monocle),"CM")
#calculate_tsne_percent(pData(CM_monocle),"CM")
#CM_monocle<-orderCells(CM_monocle,root_state=1)
plot_complex_cell_trajectory(CM_monocle,color_by="condition")+ggtitle("CM")
plot_complex_cell_trajectory(CM_monocle,color_by="State")+ggtitle("CM")
plot_complex_cell_trajectory(CM_monocle,color_by="cluster")+ggtitle("CM")

CM_monocle_df<-pData(CM_monocle)
CM_monocle_df$condition<-factor(CM_monocle_df$condition,levels=c("SHAM","TAC2W","TAC5W","TAC8W","TAC11W"))
p<-ggplot(CM_monocle_df,aes(Pseudotime,colour=(condition)))+geom_density(alpha = 0.1)+scale_colour_brewer(palette="Set1")+ggtitle("CM")
p

#EC_monocle<-orderCells(EC_monocle,root_state=6)
pData(EC_monocle)$condition<- factor(pData(EC_monocle)$condition,levels=c("SHAM","TAC2W","TAC5W","TAC8W","TAC11W"))
calculate_state_sum(pData(EC_monocle),"EC")
calculate_tsne_percent(pData(EC_monocle),"EC")
plot_complex_cell_trajectory(EC_monocle,color_by="condition")+ggtitle("EC")
plot_complex_cell_trajectory(EC_monocle,color_by="State")+ggtitle("EC")
plot_complex_cell_trajectory(EC_monocle,color_by="cluster")+ggtitle("EC")

EC_monocle_df<-pData(EC_monocle)
EC_monocle_df$condition<-factor(EC_monocle_df$condition,levels=c("SHAM","TAC2W","TAC5W","TAC8W","TAC11W"))
p<-ggplot(EC_monocle_df,aes(Pseudotime,colour=(condition)))+geom_density(alpha = 0.1)+scale_colour_brewer(palette="Set1")+ggtitle("EC")
p

FB_monocle<-orderCells(FB_monocle,root_state=4)
calculate_state_sum(pData(FB_monocle),"FB")
calculate_tsne_percent(pData(FB_monocle),"FB")

pData(FB_monocle)$condition<- factor(pData(FB_monocle)$condition,levels=c("SHAM","TAC2W","TAC5W","TAC8W","TAC11W"))

plot_complex_cell_trajectory(FB_monocle,color_by="condition")+ggtitle("FB")
plot_complex_cell_trajectory(FB_monocle,color_by="State")+ggtitle("FB")
plot_complex_cell_trajectory(FB_monocle,color_by="cluster")+ggtitle("FB")

FB_monocle_df<-pData(FB_monocle)
FB_monocle_df$condition<-factor(FB_monocle_df$condition,levels=c("SHAM","TAC2W","TAC5W","TAC8W","TAC11W"))
p<-ggplot(FB_monocle_df,aes(Pseudotime,colour=(condition)))+geom_density(alpha = 0.1)+scale_colour_brewer(palette="Set1")+ggtitle("FB")
p

MP_monocle<-orderCells(MP_monocle,root_state=3)
calculate_state_sum(pData(MP_monocle),"MP")
calculate_tsne_percent(pData(MP_monocle),"MP")

pData(MP_monocle)$condition<- factor(pData(MP_monocle)$condition,levels=c("SHAM","TAC2W","TAC5W","TAC8W","TAC11W"))

plot_complex_cell_trajectory(MP_monocle,color_by="condition")+ggtitle("MP")
plot_complex_cell_trajectory(MP_monocle,color_by="State")+ggtitle("MP")
plot_complex_cell_trajectory(MP_monocle,color_by="cluster")+ggtitle("MP")

MP_monocle_df<-pData(MP_monocle)
MP_monocle_df$condition<-factor(MP_monocle_df$condition,levels=c("SHAM","TAC2W","TAC5W","TAC8W","TAC11W"))
p<-ggplot(MP_monocle_df,aes(Pseudotime,colour=(condition)))+geom_density(alpha = 0.1)+scale_colour_brewer(palette="Set1")+ggtitle("MP")
p

#T_monocle<-orderCells(T_monocle,root_state=1)
calculate_state_sum(pData(CM_monocle),"T cell")
calculate_tsne_percent(pData(T_monocle),"T cell")

pData(T_monocle)$condition<- factor(pData(T_monocle)$condition,levels=c("SHAM","TAC2W","TAC5W","TAC8W","TAC11W"))

plot_complex_cell_trajectory(T_monocle,color_by="condition")+ggtitle("T cell")
plot_complex_cell_trajectory(T_monocle,color_by="State")+ggtitle("T cell")
plot_complex_cell_trajectory(T_monocle,color_by="cluster")+ggtitle("T cell")

T_monocle_df<-pData(T_monocle)
T_monocle_df$condition<-factor(T_monocle_df$condition,levels=c("SHAM","TAC2W","TAC5W","TAC8W","TAC11W"))
p<-ggplot(T_monocle_df,aes(Pseudotime,colour=(condition)))+geom_density(alpha = 0.1)+scale_colour_brewer(palette="Set1")+ggtitle("T_cell")
p

#Granulo_monocle<-orderCells(Granulo_monocle,root_state=4)
calculate_state_sum(pData(Granulo_monocle),"Granulo")
calculate_tsne_percent(pData(Granulo_monocle),"Granulo")

pData(Granulo_monocle)$condition<- factor(pData(Granulo_monocle)$condition,levels=c("SHAM","TAC2W","TAC5W","TAC8W","TAC11W"))

plot_complex_cell_trajectory(Granulo_monocle,color_by="condition")+ggtitle("Granulocyte")
plot_complex_cell_trajectory(Granulo_monocle,color_by="State")+ggtitle("Granulocyte")
plot_complex_cell_trajectory(Granulo_monocle,color_by="cluster")+ggtitle("Granulo")

Granulo_monocle_df<-pData(Granulo_monocle)
Granulo_monocle_df$condition<-factor(Granulo_monocle_df$condition,levels=c("SHAM","TAC2W","TAC5W","TAC8W","TAC11W"))
p<-ggplot(Granulo_monocle_df,aes(Pseudotime,colour=(condition)))+geom_density(alpha = 0.1)+scale_colour_brewer(palette="Set1")+ggtitle("Granulo_cell")
p
dev.off()

#EC4,7; FB 4,5,7;  MP 0,3; T2; Granulo 0 是心肌
CM_cell<-rownames(subset(TAC_clean_tsne$seurat@meta.data, res.1.5 %in% c(0,2,3,5,9,11,13,15,16)))

EC_cm_cell<- rownames(subset(EC_TSNE$seurat@meta.data,res.1 %in% c(4,7)))
FB_cm_cell<- rownames(subset(FB_TSNE$seurat@meta.data,res.1 %in% c(4,5,7)))
MP_cm_cell<- rownames(subset(MP_TSNE$seurat@meta.data,res.1 %in% c(0,3)))
T_cm_cell<- rownames(subset(T_TSNE$seurat@meta.data,res.1 %in% c(2)))
Granulo_cm_cell<- rownames(subset(Granulo_TSNE$seurat@meta.data,res.1 %in% c(0)))

CM_combine_cell<-c(CM_cell,EC_cm_cell,FB_cm_cell,MP_cm_cell,T_cm_cell,Granulo_cm_cell)

CM_combine_umi<- TAC_clean_tsne$seurat@raw.data[,CM_combine_cell]
CM_combine_preTSNE<- preTSNE(CM_combine_umi,TAC_cell_info,"CM_combine",regress_factor="nUMI")
CM_combine_TSNE<- doTSNE_2(CM_combine_preTSNE,"CM_combine",1:10,1.5)
write_marker_smp(CM_combine_TSNE,"CM_combine")

cell_info<-TAC_clean_tsne$seurat@meta.data
CM_combine_TSNE$seurat@meta.data$newID<- cell_info[rownames(CM_combine_TSNE$seurat@meta.data),"newID"]
CM_combine_TSNE$seurat@meta.data$value=1
CM_combine_ID<- aggregate(value~res.1.5+newID,CM_combine_TSNE$seurat@meta.data,sum)
CM_combine_ID<-CM_combine_ID[order(CM_combine_ID[,1:2]),]

CM_combine_TSNE$seurat@ident<- factor(CM_combine_TSNE$seurat@meta.data$newID)
names(CM_combine_TSNE$seurat@ident) <- rownames(CM_combine_TSNE$seurat@meta.data)


pdf("CM_combine_violin.pdf")
draw_horizon_violin(CM_combine_TSNE,FB_marker,ylabel="FB")
draw_horizon_violin(CM_combine_TSNE,Macrophage_marker,ylabel="Macrophage")
draw_horizon_violin(CM_combine_TSNE,CM_marker,ylabel="CM")
draw_horizon_violin(CM_combine_TSNE,EC_marker,ylabel="EC")
draw_horizon_violin(CM_combine_TSNE,immune_marker,ylabel="Other")
dev.off()

##180426 
setwd("/mnt/date/Project/yupeng/iCell8/Zongna_mouse_TAC/clean_data/new")

CM_EC_FB_cell<-c(CM_cell,EC_cell,FB_cell)
CM_EC_FB_umi<- TAC_clean_tsne$seurat@raw.data[,CM_EC_FB_cell]
CM_EC_FB_preTSNE<- preTSNE(CM_EC_FB_umi,TAC_cell_info,"CM_EC_FB",regress_factor="nUMI")
CM_EC_FB_TSNE<- doTSNE_2(CM_EC_FB_preTSNE,"CM_EC_FB",1:10,1.5)
write_marker_smp(CM_EC_FB_TSNE,"CM_EC_FB")


cell_info<- TAC_clean_tsne$seurat@meta.data
CM_EC_FB_TSNE$seurat@meta.data$newID<- cell_info[rownames(CM_EC_FB_TSNE$seurat@meta.data),"newID"]
CM_EC_FB_TSNE$seurat@meta.data$value=1
CM_EC_FB_combine<-aggregate(value ~res.1.5+newID,CM_EC_FB_TSNE$seurat@meta.data,sum)

#CM3 高表达 FB 类似的细胞外基质, CM5 高表达EC marker 如Cdh5, Kdr, Gpihbp1. 但两者都高表达CM marker 所以更可能是CM 
#FB3 相对高表达Cdh5, Kdr， 但仍然与FB 细胞聚在一起。 目前认为是FB.

pdf("CM_EC_FB_com_TSNE_violin.pdf")
draw_horizon_violin(CM_EC_FB_TSNE,FB_marker,ylabel="FB")
draw_horizon_violin(CM_EC_FB_TSNE,Macrophage_marker,ylabel="Macrophage")
draw_horizon_violin(CM_EC_FB_TSNE,CM_marker,ylabel="CM")
draw_horizon_violin(CM_EC_FB_TSNE,EC_marker,ylabel="EC")
draw_horizon_violin(CM_EC_FB_TSNE,immune_marker,ylabel="Other")
dev.off()



simpleCap <- function(x) {
  
  x<-tolower(x)
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
      sep="", collapse=" ")
}


marker<-read.table("zongna_marker",header=F)
marker<-apply(marker,1,simpleCap)


FB_marker<-marker[1:16]
Macrophage_marker<-marker[17:23]
CM_marker<-marker[24:33]
EC_marker<- marker[34:43]
immune_marker<-marker[44:53]

pdf("CM_combine_violin.pdf")
draw_horizon_violin(CM_combine_TSNE,FB_marker,ylabel="FB")
draw_horizon_violin(CM_combine_TSNE,Macrophage_marker,ylabel="Macrophage")
draw_horizon_violin(CM_combine_TSNE,CM_marker,ylabel="CM")
draw_horizon_violin(CM_combine_TSNE,EC_marker,ylabel="EC")
draw_horizon_violin(CM_combine_TSNE,immune_marker,ylabel="Other")
dev.off()



pdf("CM_EC_FB_TSNE_violin.pdf")
draw_horizon_violin(CM_TSNE,FB_marker,ylabel="FB")
draw_horizon_violin(CM_TSNE,Macrophage_marker,ylabel="Macrophage")
draw_horizon_violin(CM_TSNE,CM_marker,ylabel="CM")
draw_horizon_violin(CM_TSNE,EC_marker,ylabel="EC")
draw_horizon_violin(CM_TSNE,immune_marker,ylabel="Other")

draw_horizon_violin(EC_TSNE,FB_marker,ylabel="FB")
draw_horizon_violin(EC_TSNE,Macrophage_marker,ylabel="Macrophage")
draw_horizon_violin(EC_TSNE,CM_marker,ylabel="CM")
draw_horizon_violin(EC_TSNE,EC_marker,ylabel="EC")
draw_horizon_violin(EC_TSNE,immune_marker,ylabel="Other")

draw_horizon_violin(FB_TSNE,FB_marker,ylabel="FB")
draw_horizon_violin(FB_TSNE,Macrophage_marker,ylabel="Macrophage")
draw_horizon_violin(FB_TSNE,CM_marker,ylabel="CM")
draw_horizon_violin(FB_TSNE,EC_marker,ylabel="EC")
draw_horizon_violin(FB_TSNE,immune_marker,ylabel="Other")

dev.off()


CM_select_cell<- rownames(subset(CM_combine_TSNE$seurat@meta.data,res.1.5 %in% c(2,4)))
CM_select_cell_umi<- TAC_clean_tsne$seurat@raw.data[,CM_select_cell]
CM_select_preTSNE<- preTSNE(CM_select_cell_umi,TAC_cell_info,"CM_select",regress_factor="nUMI")
CM_select_TSNE<- doTSNE_2(CM_select_preTSNE,"CM_select",1:10,1.5)
write_marker_smp(CM_select_TSNE,"CM_select")

CM_select_TSNE$seurat@meta.data$newID<-sub_cluster_id[rownames(CM_select_TSNE$seurat@meta.data),11]
CM_select_TSNE$seurat@meta.data$value=1


CM_select_ID<- aggregate(value~res.1.5+newID,CM_select_TSNE$seurat@meta.data,sum)
CM_select_ID<-CM_select_ID[order(CM_select_ID[,1:2]),]



TAC_clean_tsne$seurat@meta.data$newID <-""

CM_TSNE$seurat@meta.data$newID<- paste0("CM",CM_TSNE$seurat@meta.data$res.1)
EC_TSNE$seurat@meta.data$newID<- paste0("EC",EC_TSNE$seurat@meta.data$res.1)
FB_TSNE$seurat@meta.data$newID<- paste0("FB",FB_TSNE$seurat@meta.data$res.1)
MP_TSNE$seurat@meta.data$newID<- paste0("MP",MP_TSNE$seurat@meta.data$res.1)
Granulo_TSNE$seurat@meta.data$newID<- paste0("GN",Granulo_TSNE$seurat@meta.data$res.1)
T_TSNE$seurat@meta.data$newID<- paste0("T",T_TSNE$seurat@meta.data$res.1)

sub_cluster_id <- data.frame (rbind(CM_TSNE$seurat@meta.data,EC_TSNE$seurat@meta.data,FB_TSNE$seurat@meta.data,MP_TSNE$seurat@meta.data,Granulo_TSNE$seurat@meta.data,T_TSNE$seurat@meta.data))

TAC_clean_tsne$seurat@meta.data[rownames(sub_cluster_id),"newID"] <-sub_cluster_id[,"newID"] 

TAC_clean_tsne$seurat@ident <- factor(TAC_clean_tsne$seurat@meta.data$newID)

