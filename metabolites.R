#Species indicator analysis----

library(labdsv)

setwd("/Users/elizabethmallott/Dropbox/Projects/Gut_microbiome/howler_metabolites/")
metab<-read.csv("metabolites.csv",header=T)

spec_indic<-indval(metab[,4:285],metab[,3])
gr<-spec_indic$maxcls[spec_indic$pval<=0.05]
iv<-spec_indic$indcls[spec_indic$pval<=0.05]
pv<-spec_indic$pval[spec_indic$pval<=0.05]
fr<-apply(metab[,4:285]>0, 2, sum)[spec_indic$pval<=0.05]
indvalsummary<-data.frame(group=gr,indval=iv, pvalue=pv, freq=fr)
indvalsummary<-indvalsummary[order(indvalsummary$group, -indvalsummary$indval),]
indvalsummary
write.csv(indvalsummary,"indicator_species_results.csv")

#Permanovas----

library(vegan)

metab<-read.csv("metabolites.csv",header=T)

bray_metab<-vegdist(metab[,4:285])
adonis(bray_metab~Season,strata=metab$Individaul, data=metab, permutations=5000)

##Pairwise comparisons----
library(pairwiseAdonis)

pairwise.adonis(bray_metab, factors = metab$Season, perm = 5000, 
                p.adjust.m='holm')

pairwise.adonis <- function(x,y,group, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni')
{
  library(vegan)
  
  co = combn(unique(as.character(group)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[group %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[group %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    
    ad = adonis(x1 ~ group[group %in% c(co[1,elem],co[2,elem])] ,strata=x$y, permutations=10000);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
  return(pairw.res)
} 

attach(metab)
pairwise.adonis(x=metab[,4:285],y=Individual,group=Season)

#Used metaboanalyst for heavy lifting things - removed all samples with >90% missing values

#Percent abundance of compounds found in >50% of samples----
setwd("/Users/elizabethmallott/Dropbox/Projects/Gut_microbiome/howler_metabolites/")
metab50<-read.csv("metabolites_name_match_50plus_r.csv",header=T)
library(vegan)
metab50_relab<-decostand(metab50[,4:169],method="total",MARGIN=1)
write.csv(metab50_relab,"metabolites_50plus_relab.csv")

#Metabolite networks----
library(ccrepe)
library(tidyverse)
setwd("/Users/elizabethmallott/Dropbox/Projects/Gut_microbiome/howler_metabolites/")

metab50<-read.csv("metabolites_name_match_50plus_r.csv",header=T)
metadata = metab50 %>% select(1:3)

metab50_relab_wmeta = bind_cols(metadata, metab50_relab)

metab50_relab_wmeta_rainy = metab50_relab_wmeta %>% filter(Season == "R")
metab50_relab_wmeta_int = metab50_relab_wmeta %>% filter(Season == "I")
metab50_relab_wmeta_dry = metab50_relab_wmeta %>% filter(Season == "D")

metab50_relab_rainy = metab50_relab_wmeta_rainy[,4:169]
metab50_relab_int = metab50_relab_wmeta_int[,4:169]
metab50_relab_dry = metab50_relab_wmeta_dry[,4:169]

output_all = ccrepe(metab50_relab)
output_rainy = ccrepe(metab50_relab_rainy)
output_int = ccrepe(metab50_relab_int)
output_dry = ccrepe(metab50_relab_dry)

pvalues_all = rownames_to_column(as.data.frame(output_all[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname) %>% drop_na(pvalue)
qvalues_all = rownames_to_column(as.data.frame(output_all[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname) %>% drop_na(qvalue)
zvalues_all = rownames_to_column(as.data.frame(output_all[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname) %>% drop_na(zvalue)
rhovalue_all = rownames_to_column(as.data.frame(output_all[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname) %>% drop_na(rhovalue)

output_table = bind_cols(list(pvalues_all$rowname, pvalues_all$column, pvalues_all$pvalue, 
                           qvalues_all$qvalue, zvalues_all$zvalue, rhovalue_all$rhovalue))

write.csv(output_table, "metab_network_all.csv")
write.csv(pvalues_all, "metab_network_all_pvalues.csv")
write.csv(qvalues_all, "metab_network_all_qvalues.csv")
write.csv(zvalues_all, "metab_network_all_zvalues.csv")
write.csv(rhovalue_all, "metab_network_all_rhovalues.csv")

pvalues_rainy = rownames_to_column(as.data.frame(output_rainy[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname) %>% drop_na(pvalue)
qvalues_rainy = rownames_to_column(as.data.frame(output_rainy[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname) %>% drop_na(qvalue)
zvalues_rainy = rownames_to_column(as.data.frame(output_rainy[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname) %>% drop_na(zvalue)
rhovalue_rainy = rownames_to_column(as.data.frame(output_rainy[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname) %>% drop_na(rhovalue)

output_table_rain = bind_cols(list(pvalues_rainy$rowname, pvalues_rainy$column, pvalues_rainy$pvalue, 
                              qvalues_rainy$qvalue, zvalues_rainy$zvalue, rhovalue_rainy$rhovalue))

write.csv(output_table_rain, "metab_network_rainy.csv")
write.csv(pvalues_rainy, "metab_network_rainy_pvalues.csv")
write.csv(qvalues_rainy, "metab_network_rainy_qvalues.csv")
write.csv(zvalues_rainy, "metab_network_rainy_zvalues.csv")
write.csv(rhovalue_rainy, "metab_network_rainy_rhovalues.csv")

pvalues_int = rownames_to_column(as.data.frame(output_int[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname) %>% drop_na(pvalue)
qvalues_int = rownames_to_column(as.data.frame(output_int[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname) %>% drop_na(qvalue)
zvalues_int = rownames_to_column(as.data.frame(output_int[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname) %>% drop_na(zvalue)
rhovalue_int = rownames_to_column(as.data.frame(output_int[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname) %>% drop_na(rhovalue)

output_table_int = bind_cols(list(pvalues_int$rowname, pvalues_int$column, pvalues_int$pvalue, 
                              qvalues_int$qvalue, zvalues_int$zvalue, rhovalue_int$rhovalue))

write.csv(output_table_int, "metab_network_int.csv")
write.csv(pvalues_int, "metab_network_int_pvalues.csv")
write.csv(qvalues_int, "metab_network_int_qvalues.csv")
write.csv(zvalues_int, "metab_network_int_zvalues.csv")
write.csv(rhovalue_int, "metab_network_int_rhovalues.csv")

pvalues_dry = rownames_to_column(as.data.frame(output_dry[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname) %>% drop_na(pvalue)
qvalues_dry = rownames_to_column(as.data.frame(output_dry[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname) %>% drop_na(qvalue)
zvalues_dry = rownames_to_column(as.data.frame(output_dry[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname) %>% drop_na(zvalue)
rhovalue_dry = rownames_to_column(as.data.frame(output_dry[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname) %>% drop_na(rhovalue)

output_table_dry = bind_cols(list(pvalues_dry$rowname, pvalues_dry$column, pvalues_dry$pvalue, 
                              qvalues_dry$qvalue, zvalues_dry$zvalue, rhovalue_dry$rhovalue))

write.csv(output_table_dry, "metab_network_dry.csv")
write.csv(pvalues_dry, "metab_network_dry_pvalues.csv")
write.csv(qvalues_dry, "metab_network_dry_qvalues.csv")
write.csv(zvalues_dry, "metab_network_dry_zvalues.csv")
write.csv(rhovalue_dry, "metab_network_dry_rhovalues.csv")

#16S Permanovas ----

asvs = read.csv("16S_ASV_table_wtax_t.csv", header = T)
metadata = read.csv("methanol_metadata.csv", header = T)
metadata_16S = read.csv("16S_metadata.csv", header = T)

library(tidyverse)

asvs_meta = inner_join(asvs, metadata, by = "Sample")
asvs_16S = inner_join(asvs, metadata_16S, by = "Sample")

library(vegan)

asvs_meta_ra = decostand(asvs_meta[,2:1891], method = "total")
asvs_16s_ra = decostand(asvs_16S[2:1891], method = "total")

bray_asvs_all = vegdist(asvs_16s_ra)
adonis(bray_asvs_all ~ Season, strata=asvs_16S$Individual, 
       data=asvs_16S, permutations=5000)

bray_asvs_metab = vegdist(asvs_meta_ra)
adonis(bray_asvs_metab ~ Season, strata=asvs_meta$Individual, 
       data=asvs_meta, permutations=5000)

write.csv(asvs_meta_ra, "16S_ASV_table_wmetabolite_ra.csv")
write.csv(asvs_16s_ra, "16S_ASV_table_ra.csv")

#16S Networks ----
library(ccrepe) 

#Use ASV table with ASVs present in at least 25% of samples

asvs_monthly = read.csv("16S_ASV_table_ra_monthly_wmetab_11plus.csv", header = T)

asvs_monthly_rainy = asvs_monthly %>% filter(Season == "Rain")
asvs_monthly_int = asvs_monthly %>% filter(Season == "Intermediate")
asvs_monthly_dry = asvs_monthly %>% filter(Season == "Dry")

asvs_monthly_rainy = asvs_monthly_rainy[,7:115]
asvs_monthly_int = asvs_monthly_int[,7:115]
asvs_monthly_dry = asvs_monthly_dry[,7:115]
asvs_monthly_all = asvs_monthly[7:115]

output_all = ccrepe(asvs_monthly_all, min.subj = 5, errthresh = 1e-01)
output_rainy = ccrepe(asvs_monthly_rainy, min.subj = 5, errthresh = 1e-01)
output_int = ccrepe(asvs_monthly_int, min.subj = 5, errthresh = 1e-01)
output_dry = ccrepe(asvs_monthly_dry, min.subj = 5, errthresh = 1e-01)

pvalues_all = rownames_to_column(as.data.frame(output_all[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname) %>% drop_na(pvalue)
qvalues_all = rownames_to_column(as.data.frame(output_all[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname) %>% drop_na(qvalue)
zvalues_all = rownames_to_column(as.data.frame(output_all[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname) %>% drop_na(zvalue)
rhovalue_all = rownames_to_column(as.data.frame(output_all[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname) %>% drop_na(rhovalue)

output_table = bind_cols(list(pvalues_all$rowname, pvalues_all$column, pvalues_all$pvalue, 
                              qvalues_all$qvalue, zvalues_all$zvalue, rhovalue_all$rhovalue))

write.csv(output_table, "asvs_network_all.csv")
write.csv(pvalues_all, "asvs_network_all_pvalues.csv")
write.csv(qvalues_all, "asvs_network_all_qvalues.csv")
write.csv(zvalues_all, "asvs_network_all_zvalues.csv")
write.csv(rhovalue_all, "asvs_network_all_rhovalues.csv")

pvalues_rainy = rownames_to_column(as.data.frame(output_rainy[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname)
qvalues_rainy = rownames_to_column(as.data.frame(output_rainy[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname)
zvalues_rainy = rownames_to_column(as.data.frame(output_rainy[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname)
rhovalue_rainy = rownames_to_column(as.data.frame(output_rainy[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname)

output_table_rain = bind_cols(list(pvalues_rainy$rowname, pvalues_rainy$column, pvalues_rainy$pvalue, 
                                   qvalues_rainy$qvalue, zvalues_rainy$zvalue, rhovalue_rainy$rhovalue))

write.csv(output_table_rain, "asvs_network_rainy.csv")
write.csv(pvalues_rainy, "asvs_network_rainy_pvalues.csv")
write.csv(qvalues_rainy, "asvs_network_rainy_qvalues.csv")
write.csv(zvalues_rainy, "asvs_network_rainy_zvalues.csv")
write.csv(rhovalue_rainy, "asvs_network_rainy_rhovalues.csv")

pvalues_int = rownames_to_column(as.data.frame(output_int[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname)
qvalues_int = rownames_to_column(as.data.frame(output_int[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname)
zvalues_int = rownames_to_column(as.data.frame(output_int[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname)
rhovalue_int = rownames_to_column(as.data.frame(output_int[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname)

output_table_int = bind_cols(list(pvalues_int$rowname, pvalues_int$column, pvalues_int$pvalue, 
                                  qvalues_int$qvalue, zvalues_int$zvalue, rhovalue_int$rhovalue))

write.csv(output_table_int, "asvs_network_int.csv")
write.csv(pvalues_int, "asvs_network_int_pvalues.csv")
write.csv(qvalues_int, "asvs_network_int_qvalues.csv")
write.csv(zvalues_int, "asvs_network_int_zvalues.csv")
write.csv(rhovalue_int, "asvs_network_int_rhovalues.csv")

pvalues_dry = rownames_to_column(as.data.frame(output_dry[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname)
qvalues_dry = rownames_to_column(as.data.frame(output_dry[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname)
zvalues_dry = rownames_to_column(as.data.frame(output_dry[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname)
rhovalue_dry = rownames_to_column(as.data.frame(output_dry[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname)

output_table_dry = bind_cols(list(pvalues_dry$rowname, pvalues_dry$column, pvalues_dry$pvalue, 
                                  qvalues_dry$qvalue, zvalues_dry$zvalue, rhovalue_dry$rhovalue))

write.csv(output_table_dry, "asvs_network_dry.csv")
write.csv(pvalues_dry, "asvs_network_dry_pvalues.csv")
write.csv(qvalues_dry, "asvs_network_dry_qvalues.csv")
write.csv(zvalues_dry, "asvs_network_dry_zvalues.csv")
write.csv(rhovalue_dry, "asvs_network_dry_rhovalues.csv")

#Genera networks ----

genera = read.csv("genera_table_relab_11plus_wmeta.csv")

genera_rainy = genera %>% filter(Season == "Rain")
genera_int = genera %>% filter(Season == "Intermediate")
genera_dry = genera %>% filter(Season == "Dry")

genera_rainy = genera_rainy[,11:61]
genera_int = genera_int[,11:61]
genera_dry = genera_dry[,11:61]
genera_all = genera[,11:61]

output_all = ccrepe(genera_all, min.subj = 5, errthresh = 1e-01)
output_rainy = ccrepe(genera_rainy, min.subj = 5, errthresh = 1e-01)
output_int = ccrepe(genera_int, min.subj = 5, errthresh = 1e-01)
output_dry = ccrepe(genera_dry, min.subj = 5, errthresh = 1e-01)

pvalues_all = rownames_to_column(as.data.frame(output_all[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname) %>% drop_na(pvalue)
qvalues_all = rownames_to_column(as.data.frame(output_all[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname) %>% drop_na(qvalue)
zvalues_all = rownames_to_column(as.data.frame(output_all[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname) %>% drop_na(zvalue)
rhovalue_all = rownames_to_column(as.data.frame(output_all[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname) %>% drop_na(rhovalue)

output_table = bind_cols(list(pvalues_all$rowname, pvalues_all$column, pvalues_all$pvalue, 
                              qvalues_all$qvalue, zvalues_all$zvalue, rhovalue_all$rhovalue))

write.csv(output_table, "genera_network_all.csv")
write.csv(pvalues_all, "genera_network_all_pvalues.csv")
write.csv(qvalues_all, "genera_network_all_qvalues.csv")
write.csv(zvalues_all, "genera_network_all_zvalues.csv")
write.csv(rhovalue_all, "genera_network_all_rhovalues.csv")

pvalues_rainy = rownames_to_column(as.data.frame(output_rainy[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname)
qvalues_rainy = rownames_to_column(as.data.frame(output_rainy[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname)
zvalues_rainy = rownames_to_column(as.data.frame(output_rainy[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname)
rhovalue_rainy = rownames_to_column(as.data.frame(output_rainy[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname)

output_table_rain = bind_cols(list(pvalues_rainy$rowname, pvalues_rainy$column, pvalues_rainy$pvalue, 
                                   qvalues_rainy$qvalue, zvalues_rainy$zvalue, rhovalue_rainy$rhovalue))

write.csv(output_table_rain, "genera_network_rainy.csv")
write.csv(pvalues_rainy, "genera_network_rainy_pvalues.csv")
write.csv(qvalues_rainy, "genera_network_rainy_qvalues.csv")
write.csv(zvalues_rainy, "genera_network_rainy_zvalues.csv")
write.csv(rhovalue_rainy, "genera_network_rainy_rhovalues.csv")

pvalues_int = rownames_to_column(as.data.frame(output_int[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname)
qvalues_int = rownames_to_column(as.data.frame(output_int[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname)
zvalues_int = rownames_to_column(as.data.frame(output_int[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname)
rhovalue_int = rownames_to_column(as.data.frame(output_int[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname)

output_table_int = bind_cols(list(pvalues_int$rowname, pvalues_int$column, pvalues_int$pvalue, 
                                  qvalues_int$qvalue, zvalues_int$zvalue, rhovalue_int$rhovalue))

write.csv(output_table_int, "genera_network_int.csv")
write.csv(pvalues_int, "genera_network_int_pvalues.csv")
write.csv(qvalues_int, "genera_network_int_qvalues.csv")
write.csv(zvalues_int, "genera_network_int_zvalues.csv")
write.csv(rhovalue_int, "genera_network_int_rhovalues.csv")

pvalues_dry = rownames_to_column(as.data.frame(output_dry[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname)
qvalues_dry = rownames_to_column(as.data.frame(output_dry[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname)
zvalues_dry = rownames_to_column(as.data.frame(output_dry[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname)
rhovalue_dry = rownames_to_column(as.data.frame(output_dry[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname)

output_table_dry = bind_cols(list(pvalues_dry$rowname, pvalues_dry$column, pvalues_dry$pvalue, 
                                  qvalues_dry$qvalue, zvalues_dry$zvalue, rhovalue_dry$rhovalue))

write.csv(output_table_dry, "genera_network_dry.csv")
write.csv(pvalues_dry, "genera_network_dry_pvalues.csv")
write.csv(qvalues_dry, "genera_network_dry_qvalues.csv")
write.csv(zvalues_dry, "genera_network_dry_zvalues.csv")
write.csv(rhovalue_dry, "genera_network_dry_rhovalues.csv")

#Metab micro mantel tests ----
metab_monthly = read.csv("metabolites_name_match_50plus_monthly_r_with16S_mantel.csv")
asvs_monthly = read.csv("16S_ASV_table_ra_monthly_wmetab_11plus_mantel.csv", header = T)
genus_monthly = read.csv("genera_table_ra_monthly_wmetab_11plus_mantel.csv", header = T)
family_monthly = read.csv("family_table_ra_monthly_wmetab_11plus_mantel.csv", header = T)
phyla_monthly = read.csv("phyla_table_ra_monthly_wmetab_11plus_mantel.csv", header = T)

metab_monthly_relab = decostand(metab_monthly[,4:169], method="total", MARGIN=1)
metab_monthly_relab_meta = bind_cols(metab_monthly[,1:3], metab_monthly_relab)

bray_metab_monthly = vegdist(metab_monthly_relab)
bray_asvs_monthly = vegdist(asvs_monthly[,7:115])
mantel(bray_metab_monthly, bray_asvs_monthly, 
       method = "spearman", permutations = 4999)
bray_genus_monthly = vegdist(genus_monthly[,7:102])
mantel(bray_metab_monthly, bray_genus_monthly, 
       method = "spearman", permutations = 4999)
bray_family_monthly = vegdist(family_monthly[,7:60])
mantel(bray_metab_monthly, bray_family_monthly, 
       method = "spearman", permutations = 4999)
bray_phyla_monthly = vegdist(phyla_monthly[,7:15])
mantel(bray_metab_monthly, bray_phyla_monthly, 
       method = "spearman", permutations = 4999)

##Scatterplots ----

library(graph4lg)
library(cowplot)

asvs_scatter = scatter_dist(bray_asvs_monthly, bray_metab_monthly, method = "lm") + 
  labs(y = "ASVs", x = "Metabolites")
genus_scatter = scatter_dist(bray_genus_monthly, bray_metab_monthly, method = "lm") + 
  labs(y = "Genera", x = "Metabolites")
family_scatter = scatter_dist(bray_family_monthly, bray_metab_monthly, method = "lm") + 
  labs(y = "Families", x = "Metabolites")
phyla_scatter = scatter_dist(bray_phyla_monthly, bray_metab_monthly, method = "lm") + 
  labs(y = "Phyla", x = "Metabolites")

tiff("microbe_metab_distance_scatterplots.tif", width = 8, height = 6, units = "in", 
     res = 300)
plot_grid(asvs_scatter, genus_scatter, family_scatter, phyla_scatter,
          nrow = 2, ncol = 2, align = "hv", axis = "lb")
dev.off()

#Metab micro network comparisons -----

##Monthly combined ----

metab_monthly = read.csv("metabolites_name_match_50plus_monthly_r_with16S.csv")
genera_monthly = read.csv("genera_table_relab_11plus_monthly_wmetab.csv", header = T)
asvs_monthly = read.csv("16S_ASV_table_ra_monthly_wmetab_11plus.csv", header = T)
family_monthly = read.csv("family_table_relab_11plus_monthly_wmetab.csv", header = T)
phyla_monthly = read.csv("phyla_table_relab_11plus_monthly_wmetab.csv", header = T)

metab_monthly_relab = decostand(metab_monthly[,4:169], method="total", MARGIN=1)
metab_monthly_relab_meta = bind_cols(metab_monthly[,1:3], metab_monthly_relab)

output_monthly_all_interact = ccrepe(x = metab_monthly_relab, 
                                     y = asvs_monthly[,7:115],
                                     min.subj = 5, errthresh = 1e-01)

output_monthly_all_genus_interact = ccrepe(x = metab_monthly_relab, 
                                     y = genera_monthly[,6:56],
                                     min.subj = 5, errthresh = 1e-01)

output_monthly_all_family_interact = ccrepe(x = metab_monthly_relab, 
                                           y = family_monthly[,6:35],
                                           min.subj = 5, errthresh = 1e-01)

output_monthly_all_phyla_interact = ccrepe(x = metab_monthly_relab, 
                                           y = phyla_monthly[,6:15],
                                           min.subj = 5, errthresh = 1e-01)

metab_monthly_relab_wmeta_rainy = metab_monthly_relab_meta %>% 
  filter(Season == "Rain")
metab_monthly_wmeta_int = metab_monthly_relab_meta %>% 
  filter(Season == "Intermediate")
metab_monthly_wmeta_dry = metab_monthly_relab_meta %>% 
  filter(Season == "Dry")

asvs_monthly_rainy = asvs_monthly %>% 
  filter(Season == "Rain")
asvs_monthly_int = asvs_monthly %>% 
  filter(Season == "Intermediate")
asvs_monthly_dry = asvs_monthly %>% 
  filter(Season == "Dry")

genera_monthly_rainy = genera_monthly %>% 
  filter(Season == "Rain")
genera_monthly_int = genera_monthly %>% 
  filter(Season == "Intermediate")
genera_monthly_dry = genera_monthly %>% 
  filter(Season == "Dry")

output_monthly_rainy_interact = ccrepe(x = metab_monthly_relab_wmeta_rainy[,4:169], 
                                       y = asvs_monthly_rainy[,7:115],
                                       min.subj = 5, errthresh = 1e-01)

output_monthly_dry_interact = ccrepe(x = metab_monthly_wmeta_dry[,4:169], 
                                     y = asvs_monthly_dry[,7:115],
                                     min.subj = 5, errthresh = 1e-01)

output_monthly_int_interact = ccrepe(x = metab_monthly_wmeta_int[,4:169], 
                                     y = asvs_monthly_int[,7:115],
                                     min.subj = 5, errthresh = 1e-01)

output_monthly_rainy_genus_interact = ccrepe(x = metab_monthly_relab_wmeta_rainy[,4:169], 
                                       y = genera_monthly_rainy[,6:56],
                                       min.subj = 5, errthresh = 1e-01)

output_monthly_dry_genus_interact = ccrepe(x = metab_monthly_wmeta_dry[,4:169], 
                                     y = genera_monthly_dry[,6:56],
                                     min.subj = 5, errthresh = 1e-01)

output_monthly_int_genus_interact = ccrepe(x = metab_monthly_wmeta_int[,4:169], 
                                     y = genera_monthly_int[,6:56],
                                     min.subj = 5, errthresh = 1e-01)

pvalues_monthly_all_interact = rownames_to_column(as.data.frame(
  output_monthly_all_interact[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname) %>% drop_na(pvalue)
qvalues_monthly_all_interact = rownames_to_column(as.data.frame(
  output_monthly_all_interact[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname) %>% drop_na(qvalue)
zvalues_monthly_all_interact = rownames_to_column(as.data.frame(
  output_monthly_all_interact[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname) %>% drop_na(zvalue)
rhovalue_monthly_all_interact = rownames_to_column(as.data.frame(
  output_monthly_all_interact[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname) %>% drop_na(rhovalue)

output_table = bind_cols(list(pvalues_monthly_all_interact$rowname, 
                              pvalues_monthly_all_interact$column, 
                              pvalues_monthly_all_interact$pvalue, 
                              qvalues_monthly_all_interact$qvalue, 
                              zvalues_monthly_all_interact$zvalue, 
                              rhovalue_monthly_all_interact$rhovalue))

write.csv(output_table, "interact_network.csv")
write.csv(pvalues_monthly_all_interact, "interact_network_pvalues.csv")
write.csv(qvalues_monthly_all_interact, "interact_network_qvalues.csv")
write.csv(zvalues_monthly_all_interact, "interact_network_zvalues.csv")
write.csv(rhovalue_monthly_all_interact, "interact_network_rhovalues.csv")

pvalues_monthly_all_genus_interact = rownames_to_column(as.data.frame(
  output_monthly_all_genus_interact[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname) %>% drop_na(pvalue)
qvalues_monthly_all_genus_interact = rownames_to_column(as.data.frame(
  output_monthly_all_genus_interact[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname) %>% drop_na(qvalue)
zvalues_monthly_all_genus_interact = rownames_to_column(as.data.frame(
  output_monthly_all_genus_interact[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname) %>% drop_na(zvalue)
rhovalue_monthly_all_genus_interact = rownames_to_column(as.data.frame(
  output_monthly_all_genus_interact[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname) %>% drop_na(rhovalue)

output_table = bind_cols(list(pvalues_monthly_all_genus_interact$rowname, 
                              pvalues_monthly_all_genus_interact$column, 
                              pvalues_monthly_all_genus_interact$pvalue, 
                              qvalues_monthly_all_genus_interact$qvalue, 
                              zvalues_monthly_all_genus_interact$zvalue, 
                              rhovalue_monthly_all_genus_interact$rhovalue))

write.csv(output_table, "interact_genera_network.csv")
write.csv(pvalues_monthly_all_interact, "interact_network_genera_pvalues.csv")
write.csv(qvalues_monthly_all_interact, "interact_network_genera_qvalues.csv")
write.csv(zvalues_monthly_all_interact, "interact_network_genera_zvalues.csv")
write.csv(rhovalue_monthly_all_interact, "interact_network_genera_rhovalues.csv")

pvalues_monthly_all_family_interact = rownames_to_column(as.data.frame(
  output_monthly_all_family_interact[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname) %>% drop_na(pvalue)
qvalues_monthly_all_family_interact = rownames_to_column(as.data.frame(
  output_monthly_all_family_interact[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname) %>% drop_na(qvalue)
zvalues_monthly_all_family_interact = rownames_to_column(as.data.frame(
  output_monthly_all_family_interact[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname) %>% drop_na(zvalue)
rhovalue_monthly_all_family_interact = rownames_to_column(as.data.frame(
  output_monthly_all_family_interact[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname) %>% drop_na(rhovalue)

output_table = bind_cols(list(pvalues_monthly_all_family_interact$rowname, 
                              pvalues_monthly_all_family_interact$column, 
                              pvalues_monthly_all_family_interact$pvalue, 
                              qvalues_monthly_all_family_interact$qvalue, 
                              zvalues_monthly_all_family_interact$zvalue, 
                              rhovalue_monthly_all_family_interact$rhovalue))

write.csv(output_table, "interact_family_network.csv")

pvalues_monthly_all_phyla_interact = rownames_to_column(as.data.frame(
  output_monthly_all_phyla_interact[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname) %>% drop_na(pvalue)
qvalues_monthly_all_phyla_interact = rownames_to_column(as.data.frame(
  output_monthly_all_phyla_interact[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname) %>% drop_na(qvalue)
zvalues_monthly_all_phyla_interact = rownames_to_column(as.data.frame(
  output_monthly_all_phyla_interact[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname) %>% drop_na(zvalue)
rhovalue_monthly_all_phyla_interact = rownames_to_column(as.data.frame(
  output_monthly_all_phyla_interact[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname) %>% drop_na(rhovalue)

output_table = bind_cols(list(pvalues_monthly_all_phyla_interact$rowname, 
                              pvalues_monthly_all_phyla_interact$column, 
                              pvalues_monthly_all_phyla_interact$pvalue, 
                              qvalues_monthly_all_phyla_interact$qvalue, 
                              zvalues_monthly_all_phyla_interact$zvalue, 
                              rhovalue_monthly_all_phyla_interact$rhovalue))

write.csv(output_table, "interact_phyla_network.csv")


pvalues_monthly_rainy_interact = rownames_to_column(as.data.frame(
  output_monthly_rainy_interact[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname) %>% drop_na(pvalue)
qvalues_monthly_rainy_interact = rownames_to_column(as.data.frame(
  output_monthly_rainy_interact[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname) %>% drop_na(qvalue)
zvalues_monthly_rainy_interact = rownames_to_column(as.data.frame(
  output_monthly_rainy_interact[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname) %>% drop_na(zvalue)
rhovalue_monthly_rainy_interact = rownames_to_column(as.data.frame(
  output_monthly_rainy_interact[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname) %>% drop_na(rhovalue)

output_table = bind_cols(list(pvalues_monthly_rainy_interact$rowname, 
                              pvalues_monthly_rainy_interact$column, 
                              pvalues_monthly_rainy_interact$pvalue, 
                              qvalues_monthly_rainy_interact$qvalue, 
                              zvalues_monthly_rainy_interact$zvalue, 
                              rhovalue_monthly_rainy_interact$rhovalue))

write.csv(output_table, "interact_network_rainy.csv")
write.csv(pvalues_monthly_rainy_interact, "interact_network_rainy_pvalues.csv")
write.csv(qvalues_monthly_rainy_interact, "interact_network_rainy_qvalues.csv")
write.csv(zvalues_monthly_rainy_interact, "interact_network_rainy_zvalues.csv")
write.csv(rhovalue_monthly_rainy_interact, "interact_network_rainy_rhovalues.csv")

pvalues_monthly_int_interact = rownames_to_column(as.data.frame(
  output_monthly_int_interact[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname) %>% drop_na(pvalue)
qvalues_monthly_int_interact = rownames_to_column(as.data.frame(
  output_monthly_int_interact[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname) %>% drop_na(qvalue)
zvalues_monthly_int_interact = rownames_to_column(as.data.frame(
  output_monthly_int_interact[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname) %>% drop_na(zvalue)
rhovalue_monthly_int_interact = rownames_to_column(as.data.frame(
  output_monthly_int_interact[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname) %>% drop_na(rhovalue)

output_table = bind_cols(list(pvalues_monthly_int_interact$rowname, 
                              pvalues_monthly_int_interact$column, 
                              pvalues_monthly_int_interact$pvalue, 
                              qvalues_monthly_int_interact$qvalue, 
                              zvalues_monthly_int_interact$zvalue, 
                              rhovalue_monthly_int_interact$rhovalue))

write.csv(output_table, "interact_network_int.csv")
write.csv(pvalues_monthly_int_interact, "interact_network_int_pvalues.csv")
write.csv(qvalues_monthly_int_interact, "interact_network_int_qvalues.csv")
write.csv(zvalues_monthly_int_interact, "interact_network_int_zvalues.csv")
write.csv(rhovalue_monthly_int_interact, "interact_network_int_rhovalues.csv")

pvalues_monthly_dry_interact = rownames_to_column(as.data.frame(
  output_monthly_dry_interact[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname) %>% drop_na(pvalue)
qvalues_monthly_dry_interact = rownames_to_column(as.data.frame(
  output_monthly_dry_interact[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname) %>% drop_na(qvalue)
zvalues_monthly_dry_interact = rownames_to_column(as.data.frame(
  output_monthly_dry_interact[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname) %>% drop_na(zvalue)
rhovalue_monthly_dry_interact = rownames_to_column(as.data.frame(
  output_monthly_dry_interact[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname) %>% drop_na(rhovalue)

output_table = bind_cols(list(pvalues_monthly_dry_interact$rowname, 
                              pvalues_monthly_dry_interact$column, 
                              pvalues_monthly_dry_interact$pvalue, 
                              qvalues_monthly_dry_interact$qvalue, 
                              zvalues_monthly_dry_interact$zvalue, 
                              rhovalue_monthly_dry_interact$rhovalue))

write.csv(output_table, "interact_network_dry.csv")
write.csv(pvalues_monthly_dry_interact, "interact_network_dry_pvalues.csv")
write.csv(qvalues_monthly_dry_interact, "interact_network_dry_qvalues.csv")
write.csv(zvalues_monthly_dry_interact, "interact_network_dry_zvalues.csv")
write.csv(rhovalue_monthly_dry_interact, "interact_network_dry_rhovalues.csv")

pvalues_monthly_rainy_genus_interact = rownames_to_column(as.data.frame(
  output_monthly_rainy_genus_interact[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname) %>% drop_na(pvalue)
qvalues_monthly_rainy_genus_interact = rownames_to_column(as.data.frame(
  output_monthly_rainy_genus_interact[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname) %>% drop_na(qvalue)
zvalues_monthly_rainy_genus_interact = rownames_to_column(as.data.frame(
  output_monthly_rainy_genus_interact[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname) %>% drop_na(zvalue)
rhovalue_monthly_rainy_genus_interact = rownames_to_column(as.data.frame(
  output_monthly_rainy_genus_interact[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname) %>% drop_na(rhovalue)

output_table = bind_cols(list(pvalues_monthly_rainy_genus_interact$rowname, 
                              pvalues_monthly_rainy_genus_interact$column, 
                              pvalues_monthly_rainy_genus_interact$pvalue, 
                              qvalues_monthly_rainy_genus_interact$qvalue, 
                              zvalues_monthly_rainy_genus_interact$zvalue, 
                              rhovalue_monthly_rainy_genus_interact$rhovalue))

write.csv(output_table, "interact_network_rainy_genus.csv")

pvalues_monthly_int_genus_interact = rownames_to_column(as.data.frame(
  output_monthly_int_genus_interact[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname) %>% drop_na(pvalue)
qvalues_monthly_int_genus_interact = rownames_to_column(as.data.frame(
  output_monthly_int_genus_interact[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname) %>% drop_na(qvalue)
zvalues_monthly_int_genus_interact = rownames_to_column(as.data.frame(
  output_monthly_int_genus_interact[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname) %>% drop_na(zvalue)
rhovalue_monthly_int_genus_interact = rownames_to_column(as.data.frame(
  output_monthly_int_genus_interact[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname) %>% drop_na(rhovalue)

output_table = bind_cols(list(pvalues_monthly_int_genus_interact$rowname, 
                              pvalues_monthly_int_genus_interact$column, 
                              pvalues_monthly_int_genus_interact$pvalue, 
                              qvalues_monthly_int_genus_interact$qvalue, 
                              zvalues_monthly_int_genus_interact$zvalue, 
                              rhovalue_monthly_int_genus_interact$rhovalue))

write.csv(output_table, "interact_network_int_genus.csv")

pvalues_monthly_dry_genus_interact = rownames_to_column(as.data.frame(
  output_monthly_dry_genus_interact[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname) %>% drop_na(pvalue)
qvalues_monthly_dry_genus_interact = rownames_to_column(as.data.frame(
  output_monthly_dry_genus_interact[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname) %>% drop_na(qvalue)
zvalues_monthly_dry_genus_interact = rownames_to_column(as.data.frame(
  output_monthly_dry_genus_interact[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname) %>% drop_na(zvalue)
rhovalue_monthly_dry_genus_interact = rownames_to_column(as.data.frame(
  output_monthly_dry_genus_interact[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname) %>% drop_na(rhovalue)

output_table = bind_cols(list(pvalues_monthly_dry_genus_interact$rowname, 
                              pvalues_monthly_dry_genus_interact$column, 
                              pvalues_monthly_dry_genus_interact$pvalue, 
                              qvalues_monthly_dry_genus_interact$qvalue, 
                              zvalues_monthly_dry_genus_interact$zvalue, 
                              rhovalue_monthly_dry_genus_interact$rhovalue))

write.csv(output_table, "interact_network_dry_genus.csv")

##Monthly each ----

metab_monthly = read.csv("metabolites_name_match_50plus_monthly_r_with16S.csv")
genera_monthly = read.csv("genera_table_relab_11plus_monthly_wmetab.csv", header = T)

metab_monthly_relab = decostand(metab_monthly[,4:169], method="total", MARGIN=1)
metab_monthly_relab_meta = bind_cols(metab_monthly[,1:3], metab_monthly_relab)

metab_monthly_relab_wmeta_jan = metab_monthly_relab_meta %>% 
  filter(Month == "January")
metab_monthly_wmeta_feb = metab_monthly_relab_meta %>% 
  filter(Month == "February")
metab_monthly_wmeta_mar = metab_monthly_relab_meta %>% 
  filter(Month == "March")
metab_monthly_relab_wmeta_apr = metab_monthly_relab_meta %>% 
  filter(Month == "April")
metab_monthly_wmeta_may = metab_monthly_relab_meta %>% 
  filter(Month == "May")
metab_monthly_wmeta_jun = metab_monthly_relab_meta %>% 
  filter(Month == "June")
metab_monthly_relab_wmeta_oct = metab_monthly_relab_meta %>% 
  filter(Month == "October")
metab_monthly_wmeta_nov = metab_monthly_relab_meta %>% 
  filter(Month == "November")

genera_monthly_jan = genera_monthly %>% 
  filter(Month == "January")
genera_monthly_feb = genera_monthly %>% 
  filter(Month == "February")
genera_monthly_mar = genera_monthly %>% 
  filter(Month == "March")
genera_monthly_apr = genera_monthly %>% 
  filter(Month == "April")
genera_monthly_may = genera_monthly %>% 
  filter(Month == "May")
genera_monthly_jun = genera_monthly %>% 
  filter(Month == "June")
genera_monthly_oct = genera_monthly %>% 
  filter(Month == "October")
genera_monthly_nov = genera_monthly %>% 
  filter(Month == "November")

output_monthly_jan_genus_interact = ccrepe(x = metab_monthly_relab_wmeta_jan[,4:169], 
                                             y = genera_monthly_jan[,6:56],
                                             min.subj = 3, errthresh = 1e-01)

output_monthly_feb_genus_interact = ccrepe(x = metab_monthly_wmeta_feb[,4:169], 
                                           y = genera_monthly_feb[,6:56],
                                           min.subj = 5, errthresh = 1e-01)

output_monthly_mar_genus_interact = ccrepe(x = metab_monthly_wmeta_mar[,4:169], 
                                           y = genera_monthly_mar[,6:56],
                                           min.subj = 5, errthresh = 1e-01)

output_monthly_apr_genus_interact = ccrepe(x = metab_monthly_relab_wmeta_apr[,4:169], 
                                           y = genera_monthly_apr[,6:56],
                                           min.subj = 3, errthresh = 1e-01)

output_monthly_may_genus_interact = ccrepe(x = metab_monthly_wmeta_may[,4:169], 
                                           y = genera_monthly_may[,6:56],
                                           min.subj = 5, errthresh = 1e-01)

output_monthly_jun_genus_interact = ccrepe(x = metab_monthly_wmeta_jun[,4:169], 
                                           y = genera_monthly_jun[,6:56],
                                           min.subj = 5, errthresh = 1e-01)

output_monthly_oct_genus_interact = ccrepe(x = metab_monthly_relab_wmeta_oct[,4:169], 
                                           y = genera_monthly_oct[,6:56],
                                           min.subj = 2, errthresh = 1e-01)

output_monthly_nov_genus_interact = ccrepe(x = metab_monthly_wmeta_nov[,4:169], 
                                           y = genera_monthly_nov[,6:56],
                                           min.subj = 5, errthresh = 1e-01)

pvalues_monthly_jan_interact = rownames_to_column(as.data.frame(
  output_monthly_jan_interact[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname) %>% drop_na(pvalue)
qvalues_monthly_jan_interact = rownames_to_column(as.data.frame(
  output_monthly_jan_interact[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname) %>% drop_na(qvalue)
zvalues_monthly_jan_interact = rownames_to_column(as.data.frame(
  output_monthly_jan_interact[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname) %>% drop_na(zvalue)
rhovalue_monthly_jan_interact = rownames_to_column(as.data.frame(
  output_monthly_jan_interact[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname) %>% drop_na(rhovalue)

output_table = bind_cols(list(pvalues_monthly_jan_interact$rowname, 
                              pvalues_monthly_jan_interact$column, 
                              pvalues_monthly_jan_interact$pvalue, 
                              qvalues_monthly_jan_interact$qvalue, 
                              zvalues_monthly_jan_interact$zvalue, 
                              rhovalue_monthly_jan_interact$rhovalue))

write.csv(output_table, "interact_network_jan.csv")

pvalues_monthly_feb_interact = rownames_to_column(as.data.frame(
  output_monthly_feb_interact[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname) %>% drop_na(pvalue)
qvalues_monthly_feb_interact = rownames_to_column(as.data.frame(
  output_monthly_feb_interact[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname) %>% drop_na(qvalue)
zvalues_monthly_feb_interact = rownames_to_column(as.data.frame(
  output_monthly_feb_interact[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname) %>% drop_na(zvalue)
rhovalue_monthly_feb_interact = rownames_to_column(as.data.frame(
  output_monthly_feb_interact[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname) %>% drop_na(rhovalue)

output_table = bind_cols(list(pvalues_monthly_feb_interact$rowname, 
                              pvalues_monthly_feb_interact$column, 
                              pvalues_monthly_feb_interact$pvalue, 
                              qvalues_monthly_feb_interact$qvalue, 
                              zvalues_monthly_feb_interact$zvalue, 
                              rhovalue_monthly_feb_interact$rhovalue))

write.csv(output_table, "interact_network_feb.csv")

pvalues_monthly_mar_interact = rownames_to_column(as.data.frame(
  output_monthly_mar_interact[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname) %>% drop_na(pvalue)
qvalues_monthly_mar_interact = rownames_to_column(as.data.frame(
  output_monthly_mar_interact[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname) %>% drop_na(qvalue)
zvalues_monthly_mar_interact = rownames_to_column(as.data.frame(
  output_monthly_mar_interact[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname) %>% drop_na(zvalue)
rhovalue_monthly_mar_interact = rownames_to_column(as.data.frame(
  output_monthly_mar_interact[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname) %>% drop_na(rhovalue)

output_table = bind_cols(list(pvalues_monthly_mar_interact$rowname, 
                              pvalues_monthly_mar_interact$column, 
                              pvalues_monthly_mar_interact$pvalue, 
                              qvalues_monthly_mar_interact$qvalue, 
                              zvalues_monthly_mar_interact$zvalue, 
                              rhovalue_monthly_mar_interact$rhovalue))

write.csv(output_table, "interact_network_mar.csv")

pvalues_monthly_apr_interact = rownames_to_column(as.data.frame(
  output_monthly_apr_interact[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname) %>% drop_na(pvalue)
qvalues_monthly_apr_interact = rownames_to_column(as.data.frame(
  output_monthly_apr_interact[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname) %>% drop_na(qvalue)
zvalues_monthly_apr_interact = rownames_to_column(as.data.frame(
  output_monthly_apr_interact[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname) %>% drop_na(zvalue)
rhovalue_monthly_apr_interact = rownames_to_column(as.data.frame(
  output_monthly_apr_interact[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname) %>% drop_na(rhovalue)

output_table = bind_cols(list(pvalues_monthly_apr_interact$rowname, 
                              pvalues_monthly_apr_interact$column, 
                              pvalues_monthly_apr_interact$pvalue, 
                              qvalues_monthly_apr_interact$qvalue, 
                              zvalues_monthly_apr_interact$zvalue, 
                              rhovalue_monthly_apr_interact$rhovalue))

write.csv(output_table, "interact_network_apr.csv")

pvalues_monthly_may_interact = rownames_to_column(as.data.frame(
  output_monthly_may_interact[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname) %>% drop_na(pvalue)
qvalues_monthly_may_interact = rownames_to_column(as.data.frame(
  output_monthly_may_interact[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname) %>% drop_na(qvalue)
zvalues_monthly_may_interact = rownames_to_column(as.data.frame(
  output_monthly_may_interact[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname) %>% drop_na(zvalue)
rhovalue_monthly_may_interact = rownames_to_column(as.data.frame(
  output_monthly_may_interact[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname) %>% drop_na(rhovalue)

output_table = bind_cols(list(pvalues_monthly_may_interact$rowname, 
                              pvalues_monthly_may_interact$column, 
                              pvalues_monthly_may_interact$pvalue, 
                              qvalues_monthly_may_interact$qvalue, 
                              zvalues_monthly_may_interact$zvalue, 
                              rhovalue_monthly_may_interact$rhovalue))

write.csv(output_table, "interact_network_may.csv")

pvalues_monthly_jun_interact = rownames_to_column(as.data.frame(
  output_monthly_jun_interact[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname) %>% drop_na(pvalue)
qvalues_monthly_jun_interact = rownames_to_column(as.data.frame(
  output_monthly_jun_interact[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname) %>% drop_na(qvalue)
zvalues_monthly_jun_interact = rownames_to_column(as.data.frame(
  output_monthly_jun_interact[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname) %>% drop_na(zvalue)
rhovalue_monthly_jun_interact = rownames_to_column(as.data.frame(
  output_monthly_jun_interact[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname) %>% drop_na(rhovalue)

output_table = bind_cols(list(pvalues_monthly_jun_interact$rowname, 
                              pvalues_monthly_jun_interact$column, 
                              pvalues_monthly_jun_interact$pvalue, 
                              qvalues_monthly_jun_interact$qvalue, 
                              zvalues_monthly_jun_interact$zvalue, 
                              rhovalue_monthly_jun_interact$rhovalue))

write.csv(output_table, "interact_network_jun.csv")

pvalues_monthly_oct_interact = rownames_to_column(as.data.frame(
  output_monthly_oct_interact[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname) %>% drop_na(pvalue)
qvalues_monthly_oct_interact = rownames_to_column(as.data.frame(
  output_monthly_oct_interact[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname) %>% drop_na(qvalue)
zvalues_monthly_oct_interact = rownames_to_column(as.data.frame(
  output_monthly_oct_interact[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname) %>% drop_na(zvalue)
rhovalue_monthly_oct_interact = rownames_to_column(as.data.frame(
  output_monthly_oct_interact[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname) %>% drop_na(rhovalue)

output_table = bind_cols(list(pvalues_monthly_oct_interact$rowname, 
                              pvalues_monthly_oct_interact$column, 
                              pvalues_monthly_oct_interact$pvalue, 
                              qvalues_monthly_oct_interact$qvalue, 
                              zvalues_monthly_oct_interact$zvalue, 
                              rhovalue_monthly_oct_interact$rhovalue))

write.csv(output_table, "interact_network_oct.csv")

pvalues_monthly_nov_interact = rownames_to_column(as.data.frame(
  output_monthly_nov_interact[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname) %>% drop_na(pvalue)
qvalues_monthly_nov_interact = rownames_to_column(as.data.frame(
  output_monthly_nov_interact[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname) %>% drop_na(qvalue)
zvalues_monthly_nov_interact = rownames_to_column(as.data.frame(
  output_monthly_nov_interact[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname) %>% drop_na(zvalue)
rhovalue_monthly_nov_interact = rownames_to_column(as.data.frame(
  output_monthly_nov_interact[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname) %>% drop_na(rhovalue)

output_table = bind_cols(list(pvalues_monthly_nov_interact$rowname, 
                              pvalues_monthly_nov_interact$column, 
                              pvalues_monthly_nov_interact$pvalue, 
                              qvalues_monthly_nov_interact$qvalue, 
                              zvalues_monthly_nov_interact$zvalue, 
                              rhovalue_monthly_nov_interact$rhovalue))

write.csv(output_table, "interact_network_nov.csv")

#Plant fecal metab network comparisons ----

library(vegan)
library(ccrepe)
library(tidyverse)

fecal = read.csv("metabolites_name_match_50plus_r_season.csv", header = T)
plant = read.csv("plant_metabolites_season.csv", header = T)

fecal_rainy = fecal %>% 
  filter(Season == "R")
fecal_int = fecal %>% 
  filter(Season == "I")
fecal_dry = fecal %>% 
  filter(Season == "D")

plant_rainy = plant %>% 
  filter(Season == "Rain")
plant_int = plant %>% 
  filter(Season == "Intermediate")
plant_dry = plant %>% 
  filter(Season == "Dry")

metab_all_relab = decostand(fecal[,3:168], method="total", MARGIN=1)
metab_rain_relab = decostand(fecal_rainy[,3:168], method="total", MARGIN=1)
metab_int_relab = decostand(fecal_int[,3:168], method="total", MARGIN=1)
metab_dry_relab = decostand(fecal_dry[,3:168], method="total", MARGIN=1)
plant_all_relab = decostand(plant[,6:285], method="total", MARGIN=1)
plant_rain_relab = decostand(plant_rainy[,6:285], method="total", MARGIN=1)
plant_int_relab = decostand(plant_int[,6:285], method="total", MARGIN=1)
plant_dry_relab = decostand(plant_dry[,6:285], method="total", MARGIN=1)


output_monthly_all_pf_interact = ccrepe(x = metab_all_relab, 
                                           y = plant_all_relab,
                                           min.subj = 5, errthresh = 1e-01)

output_monthly_rainy_pf_interact = ccrepe(x = metab_rain_relab, 
                                       y = plant_rain_relab,
                                       min.subj = 5, errthresh = 1e-01)

output_monthly_int_pf_interact = ccrepe(x = metab_int_relab, 
                                     y = plant_int_relab,
                                     min.subj = 5, errthresh = 1e-01)

output_monthly_dry_pf_interact = ccrepe(x = metab_dry_relab, 
                                     y = plant_dry_relab,
                                     min.subj = 5, errthresh = 1e-01)

pvalues_monthly_all_pf_interact = rownames_to_column(as.data.frame(
  output_monthly_all_pf_interact[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname) %>% drop_na(pvalue)
qvalues_monthly_all_pf_interact = rownames_to_column(as.data.frame(
  output_monthly_all_pf_interact[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname) %>% drop_na(qvalue)
zvalues_monthly_all_pf_interact = rownames_to_column(as.data.frame(
  output_monthly_all_pf_interact[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname) %>% drop_na(zvalue)
rhovalue_monthly_all_pf_interact = rownames_to_column(as.data.frame(
  output_monthly_all_pf_interact[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname) %>% drop_na(rhovalue)

output_table = bind_cols(list(pvalues_monthly_all_pf_interact$rowname, 
                              pvalues_monthly_all_pf_interact$column, 
                              pvalues_monthly_all_pf_interact$pvalue, 
                              qvalues_monthly_all_pf_interact$qvalue, 
                              zvalues_monthly_all_pf_interact$zvalue, 
                              rhovalue_monthly_all_pf_interact$rhovalue))

write.csv(output_table, "interact_network_pf.csv")

pvalues_monthly_rainy_pf_interact = rownames_to_column(as.data.frame(
  output_monthly_rainy_pf_interact[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname) %>% drop_na(pvalue)
qvalues_monthly_rainy_pf_interact = rownames_to_column(as.data.frame(
  output_monthly_rainy_pf_interact[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname) %>% drop_na(qvalue)
zvalues_monthly_rainy_pf_interact = rownames_to_column(as.data.frame(
  output_monthly_rainy_pf_interact[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname) %>% drop_na(zvalue)
rhovalue_monthly_rainy_pf_interact = rownames_to_column(as.data.frame(
  output_monthly_rainy_pf_interact[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname) %>% drop_na(rhovalue)

output_table = bind_cols(list(pvalues_monthly_rainy_pf_interact$rowname, 
                              pvalues_monthly_rainy_pf_interact$column, 
                              pvalues_monthly_rainy_pf_interact$pvalue, 
                              qvalues_monthly_rainy_pf_interact$qvalue, 
                              zvalues_monthly_rainy_pf_interact$zvalue, 
                              rhovalue_monthly_rainy_pf_interact$rhovalue))

write.csv(output_table, "interact_network_rainy_pf.csv")

pvalues_monthly_int_pf_interact = rownames_to_column(as.data.frame(
  output_monthly_int_pf_interact[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname) %>% drop_na(pvalue)
qvalues_monthly_int_pf_interact = rownames_to_column(as.data.frame(
  output_monthly_int_pf_interact[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname) %>% drop_na(qvalue)
zvalues_monthly_int_pf_interact = rownames_to_column(as.data.frame(
  output_monthly_int_pf_interact[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname) %>% drop_na(zvalue)
rhovalue_monthly_int_pf_interact = rownames_to_column(as.data.frame(
  output_monthly_int_pf_interact[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname) %>% drop_na(rhovalue)

output_table = bind_cols(list(pvalues_monthly_int_pf_interact$rowname, 
                              pvalues_monthly_int_pf_interact$column, 
                              pvalues_monthly_int_pf_interact$pvalue, 
                              qvalues_monthly_int_pf_interact$qvalue, 
                              zvalues_monthly_int_pf_interact$zvalue, 
                              rhovalue_monthly_int_pf_interact$rhovalue))

write.csv(output_table, "interact_network_int_pf.csv")

pvalues_monthly_dry_pf_interact = rownames_to_column(as.data.frame(
  output_monthly_dry_pf_interact[["p.values"]])) %>% 
  gather(key = column, value = pvalue, -rowname) %>% drop_na(pvalue)
qvalues_monthly_dry_pf_interact = rownames_to_column(as.data.frame(
  output_monthly_dry_pf_interact[["q.values"]])) %>% 
  gather(key = column, value = qvalue, -rowname) %>% drop_na(qvalue)
zvalues_monthly_dry_pf_interact = rownames_to_column(as.data.frame(
  output_monthly_dry_pf_interact[["z.stat"]])) %>% 
  gather(key = column, value = zvalue, -rowname) %>% drop_na(zvalue)
rhovalue_monthly_dry_pf_interact = rownames_to_column(as.data.frame(
  output_monthly_dry_pf_interact[["sim.score"]])) %>% 
  gather(key = column, value = rhovalue, -rowname) %>% drop_na(rhovalue)

output_table = bind_cols(list(pvalues_monthly_dry_pf_interact$rowname, 
                              pvalues_monthly_dry_pf_interact$column, 
                              pvalues_monthly_dry_pf_interact$pvalue, 
                              qvalues_monthly_dry_pf_interact$qvalue, 
                              zvalues_monthly_dry_pf_interact$zvalue, 
                              rhovalue_monthly_dry_pf_interact$rhovalue))

write.csv(output_table, "interact_network_dry_pf.csv")



