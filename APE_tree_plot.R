library (ape)
library(geiger)
setwd <- 'C:/Users/Chris/Desktop/PE/Tree_building_in_APE/'
tips<-read.table("CFEA_subfamily.txt", row.names=1)
tree <-read.nexus('C:/Users/Chris/Desktop/PE/Tree_building_in_APE/MCC_2.tre')
tree<-ladderize(tree)
outgroup        <- 'Scolecomorphus_vittatus'

plot(tree)
plot(tree, typ="fan", show.tip=FALSE, edge.color="black", tip.color="blue", cex=1) 
title("CFEA species assemblage")

#Tiplabels
host_colour<-rep("white", length(row.names(tips)))
names(host_colour)<-row.names(tips)
host_colour[tips=="1"]<-"yel
low"
host_colour[tips=="0"]<-"blue"

#Tiplabels subfamily
host_colour<-rep("white", length(row.names(tips)))
names(host_colour)<-row.names(tips)
host_colour[tips=="Microhylidae"]<-"pink"
host_colour[tips=="Pipidae"]<-"plum3"
host_colour[tips=="Gymnophiona"]<-"purple"
host_colour[tips=="Bufonidae"]<-"green"
host_colour[tips=="Hyperoliidae"]<-"yellow2"
host_colour[tips=="Hemisotidae"]<-"orange2"
host_colour[tips=="Brevicipitidae"]<-"orange4"
host_colour[tips=="Arthroleptidae"]<-"red2"
host_colour[tips=="Ranidae"]<-"navy"
host_colour[tips=="Ptychadenidae"]<-"blue"
host_colour[tips=="Phrynobatrachidae"]<-"blue2"
host_colour[tips=="Pyxicephalidae"]<-"cadetblue"
host_colour[tips=="Rhacophoridae"]<-"cadetblue3"
host_colour[tips=="Petropedetidae"]<-"lightblue"

tiplabels(pch=21, cex=0.2, col=host_colour[tree$tip.label], bg=host_colour[tree$tip.label], adj=c(0,0.5))

#Tree without ASR
par(mfrow=c(1,1))

plot(tree, show.tip.label=F, type="fan", cex=3.0) #show.tip.label=FALSE
tiplabels(pch=21, cex=3.0, col=host_colour[tree$tip.label], bg=host_colour[tree$tip.label])

#legend("bottomleft", c("small fruits < 4 cm", "large fruits > 4 cm"), pch=19,  col=c("yellow", "blue"), bty = "n")
#axisPhylo()

############################################################
############################################################
############################################################

##Chris Liedtke cladelabl###
setwd <- "C:/Users/Chris/Desktop/PE/Tree_building_in_APE/"
load("cladelabl.R")
#load data. A tree and dataset for salamanders taken and adjusted from Pyron 2014, Syst. Biol. 63(5)
tree<-read.nexus("MCC_2.tre")
tree<-ladderize(tree)
rooted_tree <-root(tree, outgroup="Scolecomorphus_vittatus", resolve.root="TRUE")
rooted_tree <-ladderize(rooted_tree)
clades<-read.table("CFEA_subfamily.txt")

#set colours:
seg.cols<-topo.colors(n=length(unique(clades$Family)))
plot(rooted_tree, show.tip.label=F, no.margin = F, type = "cladogram", rotate.tree = 0)
clade.labelR(rooted_tree,clades, f.size = 0.5, seg.cols=seg.cols, lwd=5)

###pruned_tree
pruned_rooted_tree <- drop.tip(rooted_tree, c('Scolecomorphus_vittatus','Boulengerula_boulengeri','Schistometopum_gregorii','Boulengerula_uluguruensis','Boulengerula_changamwensis','Xenopus_muelleri','Xenopus_laevis','Spelaeophryne_methneri','Leptopelis_barbouri','Arthroleptides_yakusini','Arthroleptides_martiensseni','Mertensophryne_microanotis','Mertensophryne_anotis','Mertensophryne_howelli','Mertensophryne_usambarae','Mertensophryne_loveridgei','Mertensophryne_lindneri','Schismaderma_carens','Poyntonophrynus_beiranus','Sclerophrys_brauni','Sclerophrys_xeros','Sclerophrys_gutturalis','Sclerophrys_steindachneri','Hyperolius_acuticeps','Hyperolius_pusillus','Hyperolius_argus','Hyperolius_marmoratus','Hyperolius_mariae','Hyperolius_tuberilinguis','Hyperolius_substriatus','Hyperolius_kivuensis','Hyperolius_ruvuensis','Kassina_senegalensis','Kassina_maculata','Phrynobatrachus_pakenhami','Phrynobatrachus_natalensis','Phrynobatrachus_mababiensis','Phrynobatrachus_ukingensis','Amnirana_galamensis','Hildebrandtia_ornata','Ptychadena_oxyrhynchus','Ptychadena_anchietae','Ptychadena_cf_mossambica','Ptychadena_porosissima','Ptychadena_nilotica','Ptychadena_schillukorum','Pyxicephalus_adspersus','Pyxicephalus_edulis','Nothophryne_broadleyi','Amietia_angolensis','Phrynomantis_bifasciatus','Leptopelis_barbouri','Leptopelis_broadleyi','Arthroleptis_affinis','Arthroleptis_tanneri','Callulina_kreffti','Breviceps_mossambicus'))
pruned_rooted_tree <- ladderize(pruned_rooted_tree)
plot(pruned_rooted_tree, show.tip.label=T, no.margin = T, type = "cladogram", cex=0.8)
clade.labelR(pruned_rooted_tree,clades, f.size = 0.5, seg.cols=seg.cols, lwd=5)






tree <- phylo4(tree)
plot(tree)
#############others##############
#plot

plot(tree, show.tip.label = F, no.margin = T, cex=0.25, type="phylogram")
clade.labelR(tree,clades, f.size = 0.5,seg.cols=seg.cols)

plot(tree, show.tip.label = F, no.margin = T, cex=0.25, type="cladogram", direction="upwards")
clade.labelR(tree,clades, f.size = 0.5,seg.cols=seg.cols)

plot(tree, show.tip.label = F, no.margin = T, type = "unrooted")
clade.labelR(tree,clades, f.size = 0.5,seg.cols=seg.cols)

plot(tree, show.tip.label=F, no.margin = T, type = "radial" )
clade.labelR(tree,clades, f.size = 0.5, seg.cols=seg.cols)







