# Code to plot the Zoonomia genomes + other mammals on the mammal tree of Upham et al. 2019
###########

setwd("/Users/unathan/Dropbox\ (ASU)/PROJECTS/Science_FoleyEtAl_Perspective/")

library(ape); library(geiger)
library(phytools); library(phylotate); library(phyloch) #needs to be installed from the website, not CRAN
library(dplyr); library(XML)
source("/Users/unathan/Dropbox\ (ASU)/PROJECTS/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/_R-CODE/source_functions/circleTree_plotting_functions.R")

#install.packages("remotes")
#remotes::install_github("fmichonneau/phyloch")

ZoonomiaAll<-read.csv("Zoonomia-SuppTable2.csv")
	# A total of 242 genome assemblies, representing 240 species, are included in the Zoonomia Cactus alignment. We included all non-redundant, high-quality assemblies posted on NCBI for >6 months as of March 3, 2018, or for a shorter time if an associated publication was available. One species (dog) is represented by two genomes. Due to a technical error, one genome available on NCBI (Tarsius_syrichta-2.0.1) was not included in this initial alignment, and the genome for Dipodomys stephensi was represented twice.
	# so, 241 species in the alignment V2 on UCSC: https://cglgenomics.ucsc.edu/data/cactus/

ZoonomiaSp<-do.call(rbind, strsplit(ZoonomiaAll[,"Species"]
									#ZoonomiaAll[which(ZoonomiaAll$Source!="3. Zoonomia (not in alignment)"),"Species"]
									, " "))
ZoonomiaSp_ready<-as.data.frame(paste0(ZoonomiaSp[,1],"_",ZoonomiaSp[,2]))
colnames(ZoonomiaSp_ready)<-"sciName_Zoo"

ZoonomiaAll_ready<-cbind(ZoonomiaAll, ZoonomiaSp_ready)


# Load the MamPhy MCC and taxonomy
#######
	
# load MCC tree to start with
mamMCC<-drop.tip(read.nexus("MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre"),"_Anolis_carolinensis")
plottree<-ladderize(mamMCC)

# load in clade labels
cladesDR<-read.csv("MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC_NDexp_DRstats_DRtreeLABELS.csv", header=TRUE)


# Match the Zoonomia species with the MamPhy taxonomy
#####
ZoonomiaAll_cladesDR_join<-left_join(ZoonomiaAll_ready, cladesDR, by=c("sciName_Zoo" = "SciName"))
	write.csv(ZoonomiaAll_cladesDR_join[,c("sciName_Zoo", "tiplabel")], file="Zoonomia-SuppTable2_mamPhy-joined.csv")

ZoonomiaAll_cladesDR_UNjoin<-anti_join(ZoonomiaAll_ready, cladesDR, by=c("sciName_Zoo" = "SciName"))
	write.csv(ZoonomiaAll_cladesDR_UNjoin, file="Zoonomia-SuppTable2_mamPhy-UNjoined.csv")


# Load back in the *corrected / matched* MamPhy tip names
####
ZoonomiaMatched<-read.csv(file="Zoonomia-SuppTable2_mamPhyMatched.csv")

	# load in and join with Supp Table 1
	####
	genomeStats<-read.csv(file="Zoonomia-SuppTable1.csv")
	genomeStats_ready<-	genomeStats[,c(2,7:22)]

	ZoonomiaMatched_wStats<-left_join(ZoonomiaMatched, genomeStats_ready, by=c("zoo_Accession" = "Genbank.accession"))
		write.csv(ZoonomiaMatched_wStats, file="Zoonomia-SuppTable2_mamPhyMatched_wStats.csv")

# Now get the STATS from NCBI for the other mammal genomes (non-Zoonomia)
#######
	# ncbiGenomes <- read_xml(x="assembly_result.xml", as_html=FALSE, options="NOBLANKS")#, ignoreBlanks=FALSE)#,useInternalNodes = TRUE) 

	# ROUTE 1:
	# Do this with BASH & EDirect / xtract
	######
	#	cat assembly_result.xml | xtract -pattern DocumentSummary -element AssemblyStatus Coverage AssemblyAccession SpeciesName BioSampleAccn SubmitterOrganization ContigN50 ScaffoldN50  Complete SingleCopy Duplicated Fragmented Missing > assembly_result_xtract.txt

		# Get the DATES of sequence submission
	#	cat assembly_result_260genomes-byDate.xml | xtract -pattern DocumentSummary -element AssemblyAccession SpeciesName SubmissionDate > assembly_result_260genomes-byDate_only.txt

		# get the unique taxonomy of all 2492 genome assemblies in NCBI Assembly (14 Jan 2023)
		#####
	#	cat assembly_result_14Jan2023_allAssemblies_Mammalia.xml | xtract -pattern DocumentSummary -element AssemblyStatus Coverage AssemblyAccession SpeciesName SubmissionDate BioSampleAccn SubmitterOrganization ContigN50 ScaffoldN50  Complete SingleCopy Duplicated Fragmented Missing > assembly_result_14Jan2023_allAssemblies_Mammalia_xtract.txt

	#	allAssemb<-read.csv(file="assembly_result_14Jan2023_allAssemblies_Mammalia_xtract.csv")
	#	dim(as.data.frame(table(allAssemb$SpeciesName)))[1]
	#		# 689
	#	dim(as.data.frame(table(allAssemb[which(allAssemb$AssemblyStatus=="Chromosome"), "SpeciesName"])))[1]
			# 151

			# 2492 assemblies
			# 689 unique taxa (species, subsp, hybrid)
			# 151 chromoLevel unique taxa

	# ROUTE 2:
	# "Genome" view on the NCBI web browser
	####

		# 14 Jan 2023
		#####
	#	allAssemb<-read.delim(file="assembly_result_14Jan2023_genomeView_Mammalia.tsv")
	#	dim(as.data.frame(table(allAssemb$Organism.Name)))[1]
	#		# 752
	#	dim(as.data.frame(table(allAssemb[which(allAssemb$Assembly.Level=="Chromosome"), "Organism.Name"])))[1]
	#		# 164
	#	dim(as.data.frame(table(allAssemb[which(allAssemb$Assembly.Level=="Scaffold"), "Organism.Name"])))[1]
	#		# 468
		
			# 2787 assemblies (=264 refSeq + some num of assemblies)
			# 752 unique taxa (species, subsp, hybrid)
			# 468 scaffoldLevel unique taxa
			# 164 chromoLevel unique taxa

		# 9 Feb 2023
		#####
		allAssemb<-read.delim(file="assembly_result_9Feb2023_genomeView_Mammalia.tsv")
		dim(as.data.frame(table(allAssemb$Organism.Name)))[1]
			# 754
		dim(as.data.frame(table(allAssemb[which(allAssemb$Assembly.Level=="Chromosome"), "Organism.Name"])))[1]
			# 168
		dim(as.data.frame(table(allAssemb[which(allAssemb$Assembly.Level=="Scaffold"), "Organism.Name"])))[1]
			# 489
		dim(as.data.frame(table(allAssemb[which(allAssemb$Assembly.Level=="Contig"), "Organism.Name"])))[1]

			# 2800 total genomes (= 264 RefSeq + 2,536 assemblies) <<<
			# 754 unique taxa (species, subsp, hybrid)
				# 254 contigLevel unique taxa
				# 489 scaffoldLevel unique taxa
				# 168 chromoLevel unique taxa
			# 675 species after aligning with the mammal tree taxonomy


		# subset to accession -- species -- date
		#####
		nameList<-as.character(as.data.frame(table(allAssemb$Organism.Name))[,1])
		allAssemb_sortedDate<-allAssemb[order(allAssemb$Assembly.Submission.Date),]
		allAssemb_sortedDate_full<-allAssemb_sortedDate[match(nameList, allAssemb_sortedDate$Organism.Name),]
			
			# match to mammal tree taxonomy to synonymize
			####
			allSp_split<-strsplit(allAssemb_sortedDate_full[,"Organism.Name"], " ")
			allSp_split_joined<-list()
			for(j in 1:length(allSp_split)){
				eachSp<-allSp_split[[j]]
				len<-length(eachSp)
				if(len==2){
					allSp_split_joined[[j]]<-paste0(eachSp[1],"_",eachSp[2])
				} else if(len==3){
					allSp_split_joined[[j]]<-paste0(eachSp[1],"_",eachSp[2],"_",eachSp[3])
				} else if(len>=4){
					allSp_split_joined[[j]]<-paste0(eachSp[1],"_",eachSp[2],"_",eachSp[3],"_",eachSp[4])
				}
			}
			allSp_split_joined_ready<-do.call(rbind, allSp_split_joined)
				colnames(allSp_split_joined_ready)<-"sciName_NCBI"
			allSp_all_ready<-cbind(allAssemb_sortedDate_full, allSp_split_joined_ready)
			
			# join to taxonomy
			allSp_cladesDR_join<-left_join(allSp_all_ready, cladesDR, by=c("sciName_NCBI" = "SciName"))
				write.csv(allSp_cladesDR_join, file="allAssemb_sortedDate_full_mamPhy-joined_toClean_754taxa.csv")

			# read back in the manually cleaned / matched taxonomy file
			allSp_cleaned<-read.csv(file="allAssemb_sortedDate_full_mamPhy-joined_taxonomyCleaned_754taxa.csv", header=TRUE)

			# join with the Zoonomia data (BINARY)
				ZoonomiaMatched_wStats_binary<-ZoonomiaMatched_wStats[,c("mamPhy_tiplabel","zoo_Source")] 
				# 248 species (2 repeats) -- 239 species in alignment, 241 total genomes

				allSp_cleaned_zoo<-left_join(allSp_cleaned, ZoonomiaMatched_wStats_binary, by=c("tiplabel" = "mamPhy_tiplabel"))
				write.csv(allSp_cleaned_zoo, file="allAssemb_sortedDate_full_mamPhy-joined_taxonomyCleaned_zooJoined.csv")

			# read back in the simplified / cleaned file with Zoonomia membership coded
			allSp_cleaned_zooSimp<-read.csv(file="allAssemb_sortedDate_full_mamPhy-joined_taxonomyCleaned_754taxa_zoo.csv", header=TRUE)
		
			# subset to only the unique (non-redundant species)
			allSp_cleaned_zooSimp_uniq<-allSp_cleaned_zooSimp[which(allSp_cleaned_zooSimp$redundantOrNotUsing==0),]
				# 675 species


			# Load in DATES and plot as a histo / density
			####
			#datesALL<- read.delim("assembly_result_260genomes-byDate_only.txt", header=TRUE)	
			datesALL<- allSp_cleaned_zooSimp_uniq[,c("tiplabel","Assembly.Submission.Date")]	
			dates<-strsplit(as.character(as.Date(datesALL$Assembly.Submission.Date)),"-")
			fullDate<-as.Date(datesALL$Assembly.Submission.Date)
			Year<-as.numeric(do.call(rbind,dates)[,1])
				#add Year to main data
				allSp_cleaned_zooSimp_uniq_year <- cbind.data.frame(allSp_cleaned_zooSimp_uniq, Year)

			# plot
				#plot(hist(x=yearOnly, breaks=10))

				(hist(x=fullDate, breaks="years", plot=TRUE))
				rug(x=fullDate, col=rgb(1,0,0,alpha=0.2))

				#pdf(file="plotYear_675genomes_10breaks.pdf", width=5, height=3.5)
				pdf(file="plotYear_675genomes_15breaks.pdf", width=5, height=3.5)
					#DAT<-hist(x=Year, breaks=10, plot=FALSE)
					DAT<-hist(x=Year, breaks=15, plot=FALSE)
					plot(DAT, main=NA, xaxt="n")
					#axis(side=1, at=DAT$breaks[c(1,3,5,7,9,11)]) # 10 breaks
					axis(side=1, at=DAT$breaks[c(1,5,9,13,17,21)]) # 15 breaks
				dev.off()

# ######
# # Now load this back in and PLOT
# genomesALL<- read.csv("Zoonomia-SuppTables1-2_mamPhyMatched_wStats_ALL.csv")	
# 
# # subset this to variables might plot
		#genomesReady<- genomesALL[which(genomesALL$redundant==0),c(1:2,6,20:30)]
		genomesReady<- allSp_cleaned_zooSimp_uniq_year[,c("tiplabel","Assembly.Level","Assembly.Stats.Contig.N50",
			"Assembly.Stats.Scaffold.N50","Year","zoo_Source")]


	# load trait data
	datOrig_ALL<-read.csv("/Users/unathan/Dropbox\ (ASU)/PROJECTS/VertLife_Project/MAMMALS-phylo-analyses/12_PhyloAnalyses/MamPhy_5911sp_taxonomy_traits-BM-LatLong.csv",header=TRUE)

	# Vars to plot together
		varNames<-c("tiplabel","sciName","BM", "Lat_centroid")
		varsToPlot<-left_join(x=datOrig_ALL[,varNames], y=genomesReady, 
								by=c("tiplabel" = "tiplabel"))
		rownames(varsToPlot)<-datOrig_ALL$tiplabel

	# add taxonomy
		mamTax<-read.csv(file="taxonomy_mamPhy_5911species.csv",header=TRUE)
		varsToPlot_wTax<-left_join(x=varsToPlot, y=mamTax, 
								by=c("tiplabel" = "tiplabel"))
			write.csv(varsToPlot_wTax, file="varsToPlot_wTax_mammalGenomes_675species.csv")

			# assess where new genomes added
			all_genomes   <-which(varsToPlot$Assembly.Stats.Contig.N50 >0)
			ZC_genomes 	  <-which(varsToPlot$zoo_Source=="1. Zoonomia")
			nonZC_genomes <-setdiff(all_genomes, ZC_genomes)

			varsToPlot_wTax_all<-varsToPlot_wTax[all_genomes,]
			as.data.frame(table(varsToPlot_wTax_all$higher))

				#               Var1 Freq
				# 1          Afroth.   12
				# 2 Euarchontoglires  208
				# 3   Laurasiatheria  277
				# 4      Marsupialia  168 <<< only 7 marsupial genomes before 2021; then 11 + 148 + 2 added in 2021 = 161 genomes
				# 5      Monotremata    2
				# 6          Xenart.    8

			as.data.frame(table(varsToPlot_wTax_all$ord))
				#                Var1 	GENOME SP		TOTAL SP
				# 1      AFROSORICIDA 	   3			  55
				# 2         CARNIVORA 	  74			 298
				# 3   CETARTIODACTYLA 	 128			 348
				# 4        CHIROPTERA 	  49			1287
				# 5         CINGULATA 	   3			  21
				# 6    DASYUROMORPHIA 	  63			  78
				# 7        DERMOPTERA 	   2			   2
				# 8   DIDELPHIMORPHIA 	   3			 106
				# 9     DIPROTODONTIA 	  85			 146
				# 10     EULIPOTYPHLA 	  12			 491
				# 11       HYRACOIDEA 	   2			   5
				# 12       LAGOMORPHA 	   6			  91
				# 13    MACROSCELIDEA 	   1			  19
				# 14   MICROBIOTHERIA 	   1			   1
				# 15      MONOTREMATA 	   2			   5
				# 16 NOTORYCTEMORPHIA 	   2			   2
				# 17 PAUCITUBERCULATA 	   0 << !		   7
				# 18  PERAMELEMORPHIA 	  14			  22
				# 19   PERISSODACTYLA 	  10			  24
				# 20        PHOLIDOTA 	   4			   8
				# 21           PILOSA 	   5			  12
				# 22         PRIMATES 	  83			 458
				# 23      PROBOSCIDEA 	   2			   7
				# 24         RODENTIA 	 115			2392
				# 25       SCANDENTIA 	   2			  20
				# 26          SIRENIA 	   3			   5
				# 27    TUBULIDENTATA 	   1			   1

					## PAUCITUBERCULATA is unsampled <<< !!!

				as.data.frame(table(varsToPlot_wTax_all$Assembly.Level))
				#         Var1 Freq
				# 1 Chromosome   68 + 8 more (redundant cat) = 76
				# 2     Contig  189 - 4 = 185
				# 3   Scaffold  418 - 4 = 414

				# of the scaffold genomes, how many have scaffold N50 < 1,000,000 (1 megabase)?
				length(which(varsToPlot_wTax_all$Assembly.Stats.Scaffold.N50 < 1000000))
					# 215 scaffold + 171 contig =  386 / (675-76) = 0.6444073

				# of the contig genomes, how many have contig N50 < 1,000,000 (1 megabase)?
				length(which(varsToPlot_wTax_all[which(varsToPlot_wTax_all$Assembly.Level=="Contig"),]$Assembly.Stats.Contig.N50 < 1000000))
					# 171 contigs



		# evaluating BODY MASS
		#***********************

		# summarize MASSES
		stats_ALL<-c(med=median(varsToPlot$BM), low95=quantile(varsToPlot$BM, 0.025)[[1]], up95=quantile(varsToPlot$BM, 0.975)[[1]],
					low50=quantile(varsToPlot$BM, 0.25)[[1]], up50=quantile(varsToPlot$BM, 0.75)[[1]])
 			#     med     low95      up95     low50      up50 
 			# 0.07500   0.00480 180.38840   0.02050   0.62845 


		SEQ<-varsToPlot[!is.na(varsToPlot$Assembly.Stats.Contig.N50),]
		stats_SEQ<-c(med=median(SEQ$BM), low95=quantile(SEQ$BM, 0.025)[[1]], up95=quantile(SEQ$BM, 0.975)[[1]],
					low50=quantile(SEQ$BM, 0.25)[[1]], up50=quantile(SEQ$BM, 0.75)[[1]])
			# with 260 genomes
			####
			#      med     low95      up95 
			#   2.8924    0.0076 3818.5283 

			# with 675 genomes
			####
 			#      med      low95       up95      low50       up50 
 			#  2.06670    0.00760 2392.25600    0.08015   22.41660 
		
		# compare distributions
		wilcox.test(x=varsToPlot$BM, y=SEQ$BM, alternative = "two.sided")
			# 	Wilcoxon rank sum test with continuity correction
			# 
			# data:  varsToPlot$BM and SEQ$BM
			# W = 1150239, p-value < 2.2e-16
			# alternative hypothesis: true location shift is not equal to 0

		# Plot MASS of all mammals vs sequenced species
		massAll<-log(varsToPlot$BM)
		massSeq<-log(varsToPlot[!is.na(varsToPlot$Assembly.Stats.Contig.N50),"BM"])
				plot(density(x=(massAll)))
				lines(density(x=(massSeq)))
		library(plotrix)

			pdf(file="plotMass_all-vs-sequenced_675genomes_mid50_log.pdf", width=5, height=3.5)
				# all mammals
				plot(density(massAll), col="dark grey", main="", bty="n", xlab="", ylab="",axes=F, xlim=range(massAll))#, log="x")
				polygon(density(massAll), col="light grey", border="black", bty="n",main="")#, log="x")
				x.tick <- c(round(quantile(massAll, c(0.01,0.5)),0), 0, round(quantile(massAll, c(0.99,1)),0))
				x.tick_labels <- exp(c(round(quantile(massAll, c(0.01,0.5)),0), 0, round(quantile(massAll, c(0.99,1)),0)))*1000

				axis(at=c(0,x.tick), labels=c(NA,round(x.tick_labels,0)), side=1, line=1.3, cex=2.5, lwd=1, tck=-0.05, cex.axis=1.0, mgp=c(1,1,0))
				dens.rate <- density(massAll)$y
				axis(at=c(min(dens.rate),0.45*max(dens.rate),0.9*max(dens.rate)), labels=c(0,0.45,0.9), side=2, cex=2.5, las=1, lwd=1, cex.axis=1.0, tck=-0.05, mgp=c(1,1,0))

				# sequenced mammals
				lines(density(massSeq), col="blue")
				polygon(density(massSeq), col=rgb(0,0,1,alpha=0.4), border="black", bty="n",main="")
				
				# plot CIs
				par(xpd=NA)
				plotCI(x=log(stats_ALL[1]), y= -0.01,
						ui=log(stats_ALL[5]),li=log(stats_ALL[4]),err="x",
						sfrac=0.01,gap=0,slty=par("lty"),add=TRUE,scol="black", lwd=2, pch=20)
				plotCI(x=log(stats_SEQ[1]), y= -0.02,
						ui=log(stats_SEQ[5]),li=log(stats_SEQ[4]),err="x",
						sfrac=0.01,gap=0,slty=par("lty"),add=TRUE,col=rgb(0,0,1,alpha=0.8),scol=rgb(0,0,1,alpha=0.4), lwd=2, pch=20)
			dev.off()


		# evaluating LATITUDINAL midpoint
		#***********************

		# summarize LATITUDE
		latDat_all<-abs(as.numeric(na.omit(varsToPlot$Lat_centroid)))
		stats_ALL_lat<-c(med=median(latDat_all), low95=quantile(latDat_all, 0.025)[[1]], up95=quantile(latDat_all, 0.975)[[1]],
						low50=quantile(latDat_all, 0.25)[[1]], up50=quantile(latDat_all, 0.75)[[1]])
			#        med      low95       up95      low50       up50 
			# 15.3515580  0.5172062 52.4364595  5.8689934 28.0568668 

		latDat_SEQ<-abs(as.numeric(na.omit(varsToPlot[!is.na(varsToPlot$Assembly.Stats.Contig.N50),"Lat_centroid"])))
		stats_SEQ_lat<-c(med=median(latDat_SEQ), low95=quantile(latDat_SEQ, 0.025)[[1]], up95=quantile(latDat_SEQ, 0.975)[[1]], 
						low50=quantile(latDat_SEQ, 0.25)[[1]], up50=quantile(latDat_SEQ, 0.75)[[1]])
			# with 260 genomes
			####
			#       med      low95       up95 
			# 22.1117487  0.7488159 64.5050064 

			# with 675 genomes
			####
			#        med      low95       up95      low50       up50 
			# 22.3010853  0.6354131 58.6851888  8.3930598 34.9012435 

		# compare distributions
		wilcox.test(x=latDat_all, y=latDat_SEQ, alternative = "two.sided")
			# 	Wilcoxon rank sum test with continuity correction
			# 
			# data:  latDat_all and latDat_SEQ
			# W = 1404678, p-value = 3.09e-14
			# alternative hypothesis: true location shift is not equal to 0

		# Plot LATITUDE of all mammals vs sequenced species
		
			pdf(file="plotLatitude_all-vs-sequenced_abs_675genomes_mid50.pdf", width=5, height=3.5)
				# all mammals
				plot(density(latDat_all), col="dark grey", main="", bty="n", xlab="", ylab="",axes=F, xlim=range(latDat_all))
				polygon(density(latDat_all), col="light grey", border="black", bty="n",main="")
				x.tick <- c(0,23.5,60)
				axis(at=c(0,x.tick), labels=c(NA,round(x.tick,2)), side=1, line=1.3, cex=2.5, lwd=1, tck=-0.05, cex.axis=1.0, mgp=c(1,1,0))
				dens.rate <- density(latDat_all)$y
				axis(at=c(min(dens.rate),0.45*max(dens.rate),0.9*max(dens.rate)), labels=c(0,0.45,0.9), side=2, cex=2.5, las=1, lwd=1, cex.axis=1.0, tck=-0.05, mgp=c(1,1,0))

				# sequenced mammals
				lines(density(latDat_SEQ), col="blue")
				polygon(density(latDat_SEQ), col=rgb(0,0,1,alpha=0.4), border="black", bty="n",main="")
				
				# plot CIs
				par(xpd=NA)
				plotCI(x=(stats_ALL_lat[1]), y= -0.0015,
						ui=(stats_ALL_lat[5]),li=(stats_ALL_lat[4]),err="x",
						sfrac=0.01,gap=0,slty=par("lty"),add=TRUE,scol="black", lwd=2, pch=20)
				plotCI(x=(stats_SEQ_lat[1]), y= -0.0035,
						ui=(stats_SEQ_lat[5]),li=(stats_SEQ_lat[4]),err="x",
						sfrac=0.01,gap=0,slty=par("lty"),add=TRUE,col=rgb(0,0,1,alpha=0.8),scol=rgb(0,0,1,alpha=0.4), lwd=2, pch=20)
			dev.off()


		


# ALL PLOTTED TOGETHER
#####
library(viridis)
varCols<-magma(20)[c(2,11,7)]
	#latCols<-rgb(col2rgb(magma(20)[7])[1,][[1]], col2rgb(magma(20)[7])[2,][[1]], col2rgb(magma(20)[7])[3,][[1]])
	fc <- colorRampPalette(c("mistyrose", "darkorchid4"))
	#plot(rep(1, 30),col = fc(30), pch = 19, cex = 3)
	latCols<-c(fc(30)[10],fc(30)[27])

labSize<-1.5

	rootAge=max(node.depth.edgelength(plottree))

	# for the grey order boxes
	orderOfTips<-obj$xy$yy[seq_len(obj$Ntip)]
	tipDetails1<-cbind.data.frame(plottree$tip.label,do.call(rbind,strsplit(plottree$tip.label,split="_")))
	colnames(tipDetails1)<-c("tiplabel","gen","sp","fam","ord")
	toJoin<-cladesDR[,c("tiplabel","clade")]
	tipDetails<-left_join(tipDetails1,toJoin,by="tiplabel")

	blackOrds<-c(#"LAGOMORPHA", #"SCANDENTIA",
				"PRIMATES","CETARTIODACTYLA", "EULIPOTYPHLA",#"PERISSODACTYLA","PHOLIDOTA","AFROSORICIDA","TUBULIDENTATA","HYRACOIDEA","CINGULATA",
				"DIDELPHIMORPHIA")#,"PAUCITUBERCULATA")
	grayOrds<-c("RODENTIA", "CHIROPTERA", "CARNIVORA",#"DERMOPTERA",
				 #"MACROSCELIDEA","PROBOSCIDEA","SIRENIA","PILOSA",
				"Australidelphia")#"DIPROTODONTIA","DASYUROMORPHIA","NOTORYCTEMORPHIA","PERAMELEMORPHIA","MICROBIOTHERIA")#,"MONOTREMATA")
	
	majorOrds<-c("CHIROPTERA", "CARNIVORA","DIPROTODONTIA","DIDELPHIMORPHIA",
				"RODENTIA", "PRIMATES", "CETARTIODACTYLA","EULIPOTYPHLA", "DASYUROMORPHIA")

	#grayCLADES<-c("Mouse-related","Guinea_pig-related","PRIMATES","Yinpterochiroptera","Whippomorpha","EULIPOTYPHLA","Marsupialia")
	grayCLADES<-c("Mouse-related","Guinea_pig-related","Platyrrhini","Strepsirrhini","Yangochiroptera",
				"Ruminantia","Suina","Caniformes",#"Pholidota",
				"Soricidae", "Talpidae", "Marsupialia")




			pdf(file="plotMamPhy_withGenomesSampled_andBodySize.pdf", width=16, height=14)
			pdf(file="plotMamPhy_withGenomesSampled_binary.pdf", width=16, height=14)
			pdf(file="plotMamPhy_withGenomesSampled_binary_contig.pdf", width=16, height=14)
			pdf(file="plotMamPhy_withGenomesSampled_binary-thicker_675genomes_ZC.pdf", width=16, height=14)
			pdf(file="plotMamPhy_withGenomesSampled_binary-thicker_675genomes_ZC-prePost.pdf", width=16, height=14)

			#quartz(width=16, height=14)
			#pdf(file="newFig3_test_tipLabels.pdf", width=16, height=60)
			# 	par(oma=c(0,0,2,20, xpd=NA) #‘c(bottom, left, top, right)
			# 	layout(matrix(c(1:2), 1, 2, byrow = TRUE), widths=c(8,12), heights=rep(14,4))
			 	layout(matrix(c(1:5), 1, 5, byrow = TRUE), widths=c(4,4,4,4,4), heights=rep(14,5))

					# plot
					###
			#	   pdf(file="treeForNewFig1_wDR_wBAMMshifts_wlabel_wLeg.pdf", width=8, height=14)
				    par(oma=c(0,0,2,4), mar=c(1,1,1,1)) #‘c(bottom, left, top, right)

					# plot dummy tree
					obj <- plot2.phylo(plottree, show.tip.label=FALSE, cex=0.05, tip.color="black", #x.lim=c(0,roots[j]),
						label.offset=0.1, type="phylogram", edge.width=0.8, no.margin=TRUE, root.edge=TRUE, edge.color="white")
					# plot CLADES
					for(j in 1:length(grayCLADES)){
						if(j==11){		
							ordTips<-which(tipDetails$ord=="PAUCITUBERCULATA" | tipDetails$ord=="DIDELPHIMORPHIA" | tipDetails$ord=="MICROBIOTHERIA" | tipDetails$ord=="NOTORYCTEMORPHIA" | tipDetails$ord=="PERAMELEMORPHIA" | tipDetails$ord=="DASYUROMORPHIA" | tipDetails$ord=="DIPROTODONTIA")
						} else {
							ordTips<-which(tipDetails$clade==grayCLADES[j])
						}			
						par(new=T, xpd=NA)
						#rect(xleft=0, ybottom=min(orderOfTips[ordTips]), xright=obj$xy$xx[1]+9, ytop=max(orderOfTips[ordTips]), col=gray(0.5, alpha=0.2), border=gray(0.5, alpha=0.2))
						rect(xleft=obj$xy$xx[1], ybottom=min(orderOfTips[ordTips]), xright=obj$xy$xx[1]+500, ytop=max(orderOfTips[ordTips]), col=gray(0.5, alpha=0.2), border=NA)
					}

					# plot ORDERS
					for(j in 1:length(blackOrds)){
							blackOrdTips<-which(tipDetails$ord==blackOrds[j])
						par(new=T, xpd=NA)
						#rect(xleft=0, ybottom=min(orderOfTips[ordTips]), xright=obj$xy$xx[1]+9, ytop=max(orderOfTips[ordTips]), col=gray(0.5, alpha=0.2), border=gray(0.5, alpha=0.2))
						rect(xleft=obj$xy$xx[1], ybottom=min(orderOfTips[blackOrdTips]), xright=obj$xy$xx[1]+10, ytop=max(orderOfTips[blackOrdTips]), col=gray(0.2), border=NA)
						}
					for(j in 1:length(grayOrds)){
						if(j==4){
							grayOrdTips<-which(tipDetails$ord=="DIPROTODONTIA" | tipDetails$ord=="DASYUROMORPHIA" | tipDetails$ord=="NOTORYCTEMORPHIA" | tipDetails$ord=="PERAMELEMORPHIA" | tipDetails$ord=="MICROBIOTHERIA")
						} else {
							grayOrdTips<-which(tipDetails$ord==grayOrds[j])
						}
						par(new=T, xpd=NA)
						#rect(xleft=0, ybottom=min(orderOfTips[ordTips]), xright=obj$xy$xx[1]+9, ytop=max(orderOfTips[ordTips]), col=gray(0.5, alpha=0.2), border=gray(0.5, alpha=0.2))
						rect(xleft=obj$xy$xx[1], ybottom=min(orderOfTips[grayOrdTips]), xright=obj$xy$xx[1]+10, ytop=max(orderOfTips[grayOrdTips]), col=gray(0.7), border=NA)
						}
					

					par(new=T)
				    # Plot the real tree 
					obj <- plot2.phylo(plottree, show.tip.label=FALSE, cex=0.05, tip.color="black", #x.lim=c(0,roots[j]),
						label.offset=0.1, type="phylogram", edge.width=0.3, no.margin=TRUE, root.edge=TRUE, edge.color=
						#as.matrix(shiftEdgeCols))
						"black")#as.matrix(reconColors))
					axis(side=3, line=-2, at=rootAge-c(0,50,100,150), labels=FALSE, cex.axis=labSize) #c(0,50,100,150)
					text(y=c(6120,6120,6140,6140), x=rootAge-c(0,50,100,150), labels=c(0,50,100,150), cex=labSize, srt=270)


	# SEPARATE
	####
	# SEGMENTS STYLE
	############
	#	plottree<-ladderize(mamMCC)
	#	# COLS FOR TRAITS
	#	library(viridis)
	#	varCols<-viridis(10)[c(2,5,8)]

	# GENOME to plot...
		all_genomes       <-which(varsToPlot$Assembly.Stats.Contig.N50 >0)
		ZC_genomes_new 	  <-which(varsToPlot$zoo_Source=="1. Zoonomia") 
			# 121 species newly sequenced
		ZC_genomes_old    <-which(varsToPlot$zoo_Source=="2. Existing assembly in 241-way")
			# other 118 species of the 239 species in the 241-way alignment
		#nonZC_genomes <-setdiff(all_genomes, ZC_genomes)
			pre2019_genomes   <-which(varsToPlot$Year<=2019) # 361 genomes (some of which are Zoo new)
			post2019_genomes   <-which(varsToPlot$Year>2019) # 314 genomes (some of which are Zoo new)
		preZoo_genomes    <-setdiff(pre2019_genomes,ZC_genomes_new) # 242 genomes
		postZoo_genomes   <-setdiff(post2019_genomes, ZC_genomes_new) # 312 genomes

			# total == 242 preZoo + 121 zooNew + 312 postZoo == 675, awesome

		VAR_all<-rep(0,length(varsToPlot[,1]))
		names(VAR_all)<-rownames(varsToPlot)
		VAR_all[all_genomes]<-1

		VAR_ZC_new<-rep(0,length(varsToPlot[,1]))
		names(VAR_ZC_new)<-rownames(varsToPlot)
		VAR_ZC_new[ZC_genomes_new]<-1

		VAR_ZC_old<-rep(0,length(varsToPlot[,1]))
		names(VAR_ZC_old)<-rownames(varsToPlot)
		VAR_ZC_old[ZC_genomes_old]<-1

		VAR_preZC<-rep(0,length(varsToPlot[,1]))
		names(VAR_preZC)<-rownames(varsToPlot)
		VAR_preZC[preZoo_genomes]<-1

		VAR_postZC<-rep(0,length(varsToPlot[,1]))
		names(VAR_postZC)<-rownames(varsToPlot)
		VAR_postZC[postZoo_genomes]<-1


		    # Plot the false tree 
			obj <- plot2.phylo(plottree, show.tip.label=FALSE, cex=0.1, tip.color="black", #x.lim=c(0,roots[j]),
				label.offset=0.1, type="phylogram", edge.width=0.8, no.margin=TRUE, root.edge=TRUE, edge.color=NA)#"black")

			# for the TIP order 
			orderOfTips<-obj$xy$yy[seq_len(obj$Ntip)]
	#		eqLine <- obj$xy$xx[1]+10 
	#		# wEXTRA
	#		EXTRA <- 145
	#		eqLine <- obj$xy$xx[1]+10 +EXTRA
			# wLESS
			LESS <- -210
			eqLine <- obj$xy$xx[1]+40 +LESS

			buffer<-150#10#150
			## plot coord axis
			#####
			#segments( x0=eqLine, y0=1, x1=eqLine, y1=5911+30, lty=1, col=gray(0.2), lwd=1)
			#	text( x=eqLine, y=5911+buffer, labels="0", cex=labSize, srt=270)
			#segments( x0=eqLine+50, y0=1, x1=eqLine+50, y1=5911+30, lty=1, col=gray(0.2), lwd=0.5)
			#	text( x=eqLine+30, y=5911+buffer, labels="1", cex=labSize, srt=270)
				text( x=eqLine+15, y=5911+buffer, labels="2002-2018", cex=labSize, srt=270)
				text( x=eqLine+55, y=5911+buffer, labels="Zoonomia\n(2019)", cex=labSize, srt=270)
				text( x=eqLine+95, y=5911+buffer, labels="2020-2023", cex=labSize, srt=270)

			# get coords PER TIP + plot those:
			recordTipNumber_all<-c()
			# ALL GENOMES << doing this to record
			####
			for(j in 1:length(orderOfTips)){
				tipInTree<-orderOfTips[j]
				species <- plottree$tip.label[j]
				sppDat <- VAR_all[which(names(VAR_all)==species)]
				if(sppDat[[1]] > 0){
					recordTipNumber_all[j]<-tipInTree
				} else { next }

			#	valX <- eqLine + sppDat[[1]]*30 
			#	segments( x0=eqLine, y0=tipInTree, x1=valX, y1=tipInTree, lty=1, lwd=0.7, col=rgb(0,0,1,alpha=0.8)) #"darkorange"
			}

			recordTipNumber_preZC<-c()
			# pre-Zoonomia GENOMES >> to plot, row1
			####
			for(j in 1:length(orderOfTips)){
				tipInTree<-orderOfTips[j]
				species <- plottree$tip.label[j]
				sppDat <- VAR_preZC[which(names(VAR_preZC)==species)]
				if(sppDat[[1]] > 0){
					recordTipNumber_preZC[j]<-tipInTree
				} else { next }

			#	#valX <- eqLine + VAR_nonZC[k][[1]]*30 
				valX <- eqLine + sppDat[[1]]*30 
				segments( x0=eqLine, y0=tipInTree, x1=valX, y1=tipInTree, lty=1, lwd=0.7, col=rgb(0,0,1,alpha=0.7)) #"darkorange"
			}

			recordTipNumber_ZC<-c()
			# Zoonomia *NEW* GENOMES >> to plot, row2
			####
			for(j in 1:length(orderOfTips)){
				tipInTree<-orderOfTips[j]
				species <- plottree$tip.label[j]
				sppDat <- VAR_ZC_new[which(names(VAR_ZC_new)==species)]
				if(sppDat[[1]] > 0){
					recordTipNumber_ZC[j]<-tipInTree
				} else { next }

				#valX <- eqLine + VAR_ZC[k][[1]]*30 
				valX <- eqLine + sppDat[[1]]*30 + 40
				segments( x0=eqLine+40, y0=tipInTree, x1=valX, y1=tipInTree, lty=1, lwd=0.9, col="red")#magma(20)[11]) #"darkorange"
			}

			# Zoonomia *OLD* GENOMES >> to plot, row2
			####
			for(j in 1:length(orderOfTips)){
				tipInTree<-orderOfTips[j]
				species <- plottree$tip.label[j]
				sppDat <- VAR_ZC_old[which(names(VAR_ZC_old)==species)]
				if(sppDat[[1]] > 0){
					recordTipNumber_ZC[j]<-tipInTree
				} else { next }

				#valX <- eqLine + VAR_ZC[k][[1]]*30 
				valX <- eqLine + sppDat[[1]]*30 + 40
				segments( x0=eqLine+40, y0=tipInTree, x1=valX, y1=tipInTree, lty=1, lwd=0.7, col=gray(0.5)) #"darkorange"
			}

			recordTipNumber_postZC<-c()
			# post-Zoonomia GENOMES >> to plot, row3
			####
			for(j in 1:length(orderOfTips)){
				tipInTree<-orderOfTips[j]
				species <- plottree$tip.label[j]
				sppDat <- VAR_postZC[which(names(VAR_postZC)==species)]
				if(sppDat[[1]] > 0){
					recordTipNumber_postZC[j]<-tipInTree
				} else { next }

				#valX <- eqLine + VAR_ZC[k][[1]]*30 
				valX <- eqLine + sppDat[[1]]*30 + 80
				segments( x0=eqLine+80, y0=tipInTree, x1=valX, y1=tipInTree, lty=1, lwd=0.7, col=rgb(0.5,0,0.5,alpha=0.7)) #"darkorange"
			}

dev.off()



	# get the DENSITY PLOT
		# ALL 
		####
		#pdf(file="plotMamPhy_withGenomesSampled_densityOfSampled2_675genomes_all.pdf", width=16, height=5)
		#pdf(file="plotMamPhy_withGenomesSampled_densityOfSampled2_675genomes_ZC.pdf", width=16, height=5)
		#pdf(file="plotMamPhy_withGenomesSampled_densityOfSampled2_675genomes_nonZC.pdf", width=16, height=5)
		pdf(file="plotMamPhy_withGenomesSampled_densityOfSampled2_675genomes_all_ZC_prePost.pdf", width=16, height=5)
		#pdf(file="plotMamPhy_withGenomesSampled_densityOfSampled2_675genomes_ZC_prePost.pdf", width=16, height=5)
			DAT<-recordTipNumber_all
			
			recordTipNumber_only<-as.numeric(na.omit(DAT))			
			plot(density(recordTipNumber_only), col="black", main="", bty="n", xlab="", ylab="",axes=F, 
				xlim=range(recordTipNumber_only), lwd=4)
			dens.rate <- density(recordTipNumber_only)$y
			axis(at=c(min(dens.rate),0.45*max(dens.rate),0.9*max(dens.rate)), labels=c(0,0.45,0.9), side=2, 
				cex=2.5, las=1, lwd=1, cex.axis=1.5, tck=-0.05, mgp=c(1,1,0))

			par(new=TRUE)
			DAT<-recordTipNumber_ZC
				recordTipNumber_only<-as.numeric(na.omit(DAT))			
				plot(density(recordTipNumber_only), col="red", main="", bty="n", xlab="", ylab="",axes=F, 
					xlim=range(recordTipNumber_only), lwd=2)
			
			par(new=TRUE)
			DAT<-recordTipNumber_preZC
				recordTipNumber_only<-as.numeric(na.omit(DAT))			
				plot(density(recordTipNumber_only), col="blue", main="", bty="n", xlab="", ylab="",axes=F, 
					xlim=range(recordTipNumber_only), lwd=2)
						
			par(new=TRUE)
			DAT<-recordTipNumber_postZC
				recordTipNumber_only<-as.numeric(na.omit(DAT))			
				plot(density(recordTipNumber_only), col="purple", main="", bty="n", xlab="", ylab="",axes=F, 
					xlim=range(recordTipNumber_only), lwd=2)

		dev.off()



# *************************
# If wanted to plot the contig N50 next to the species (not so interesting, turns out)
#######

	# CONTIG SIZE
		VAR<-log10(varsToPlot[,"Contig.N50.bp"]*1000)#-min(log(varsToPlot[,"DispDist_kmMAX_final"]*1000))
			# 0 ==  0.00447 km max dispersal distance == 4.5 meters... 1.001499 == exp(min(log(varsToPlot[,"DispDist_kmMAX_final"]*1000))/1000)
			# Melomys_rubicola_MURIDAE_RODENTIA == min
		names(VAR)<-rownames(varsToPlot)
		VAR[is.na(VAR)] <- 0

		treeDat<-treedata(plottree, as.data.frame(VAR))
		trait<-as.vector(treeDat$data)
		names(trait)<-rownames(treeDat$data)
		
		magma20<-magma(20)
		FACTOR<-15

		    # Plot the false tree 
			obj <- plot2.phylo(plottree, show.tip.label=FALSE, cex=0.1, tip.color="black", #x.lim=c(0,roots[j]),
				label.offset=0.1, type="phylogram", edge.width=0.8, no.margin=TRUE, root.edge=TRUE, edge.color=NA)#"black")

			# for the TIP order 
			orderOfTips<-obj$xy$yy[seq_len(obj$Ntip)]

#			eqLine <- obj$xy$xx[1]+40
			# wLESS
			LESS <- -180
			eqLine <- obj$xy$xx[1]+40 +LESS

			ticksAt <- c(0, 1, 2, 3, 4, 5, 6, 7)
			ticksLabels <- c("0", round((10^ticksAt[2])/1000,2), round((10^ticksAt[3])/1000,2), round((10^ticksAt[4])/1000,2), round((10^ticksAt[5])/1000,2), round((10^ticksAt[6])/1000,2), round((10^ticksAt[7])/1000,2), round((10^ticksAt[8])/1000,2))
			#ticksLabels <- c("0", "0.005", "0.025", "13.3", "727.5", "39,721")
			#ticksLabels_log <- c("1.5", "9.5", "17.5")

				# > exp((8+1.497388))/1000
				# [1] 13.32488
				# > exp((16+1.497388))/1000
				# [1] 39720.9

			buffer<-150#10#150
			# plot coord axis
			segments( x0=eqLine, y0=1, x1=eqLine, y1=5911+30, lty=1, col=gray(0.2), lwd=1)
				text( x=eqLine, y=5911+buffer, labels=ticksLabels[1], cex=labSize, srt=270)
			
			for(j in 2:length(ticksAt)){
			mod<-ticksAt[j]*FACTOR
			segments( x0=eqLine+mod, y0=1, x1=eqLine+mod, y1=5911+30, lty=2, col=gray(0.5), lwd=0.5)
				text( x=eqLine+mod, y=5911+buffer, labels=ticksLabels[j],  cex=labSize, srt=270)
			}

			# get coords PER TIP
			for(j in 1:length(orderOfTips)){
				tipInTree <- orderOfTips[j]
				species <- plottree$tip.label[j]
				if(is.na(varsToPlot[species,"Contig.N50.bp"])){ next
				} else if(varsToPlot[species,"Contig.N50.bp"] <= 1 & varsToPlot[species,"Contig.N50.bp"] > 0.1){
					vagCol<-magma20[17]
				} else if(varsToPlot[species,"Contig.N50.bp"] <= 0.1 ){
					vagCol<-magma20[20]
				} else {
					vagCol<-magma20[13]
				}
				sppDat <- as.data.frame(trait[which(names(trait)==species)])
				valX <- eqLine + sppDat[,1]*FACTOR #"lat_min"
				segments( x0=eqLine, y0=tipInTree, x1=valX, y1=tipInTree, lty=1, col=vagCol, lwd=0.3)
			}


dev.off()




#		# Code to plot the tree from Foley et al. 2022
#		###########
#
#
#		treeFoley<-read.nexus(file="phyloFile_fromFoleyEtAl_inReview_DataS2.nex")
#		#treeFoley<-read.beast(file="phyloFile_fromFoleyEtAl_inReview_DataS2.nex")
#
#		obj<-scan(file="phyloFile_fromFoleyEtAl_inReview_DataS2.nex",n=1,skip=3,what="character")
#		treeFoley1<-read.tree(text=obj)
#		obj2<-strsplit(obj,"tree1=")
#		treeFoley2<-parse_annotated(str=obj2[[1]][2], format="newick")
#
#			# parse the error bars
#			treeFoley2_errRaw<-treeFoley2$node.comment
#			parse1<-strsplit(treeFoley2_errRaw,",")
#			parse2<-strsplit(do.call(rbind,parse1)[,1],"=")
#			min95<-do.call(rbind,parse1)[,2]
#			max95<-do.call(rbind,parse2)[,2]
#
#		# addErrorBars -- tree is *ALREADY SCALED* to the root (Placentalia) age of 103.80 Ma -- multiply by 100 to show:
#		treeFoley2$height_95_HPD_MIN <- as.numeric(min95[(length(treeFoley2$tip.label)+1):length(treeFoley2$node.comment)])*100 # multiply to put in Millions of Years
#		treeFoley2$height_95_HPD_MAX <- as.numeric(max95[(length(treeFoley2$tip.label)+1):length(treeFoley2$node.comment)])*100 # multiply to put in Millions of Years
#		treeFoley2$edge.length <- treeFoley2$edge.length*100 # multiply to put in Millions of Years
#
#
#			# PLOT w error bars (as scaled to 1.0 -- all that they gave in the paper...)
#			pdf(file=paste0("phyloFile_fromFoleyEtAl_inReview_DataS2_PLOTTED_wTips.pdf"), width=8.5, height=22, onefile=TRUE)
#
#				tree<-treeFoley2
#
#				#quartz(width=8.5, height=22)
#				plot(ladderize(tree), show.tip.label=TRUE, cex=0.5)#, x.lim=c(-7,40), tip.color=tipColors, edge.width=2)#, y.lim=c(130,3950))
#
#				HPDbars(tree, label="height_95_HPD", broken=T, lwd=2, col=hsv(0.65,1,1,alpha=0.7))
#				#node.support(branching.times(tree), mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.8)
#				
#				data(gradstein04)
#				axisGeo(GTS = gradstein04, unit = c("epoch","period"), cex = 0.8, ages=FALSE, gridty=3, gridcol="grey50")
#				axisPhylo()#cex.axis=1.8, pos=-2, mgp=c(0,0.2,0))
#
#			dev.off()
#
#			# PLOT simple
#			pdf(file=paste0("phyloFile_fromFoleyEtAl_inReview_DataS2_PLOTTED_only.pdf"), width=5, height=7, onefile=TRUE)
#				plot(treeFoley2, show.tip.label=FALSE)
#			dev.off()
#
#			# PLOT simple -- with KPG and PETM
#			pdf(file=paste0("phyloFile_fromFoleyEtAl_inReview_DataS2_PLOTTED_only_KPG-PETM_simple.pdf"), width=5, height=7, onefile=TRUE)
#				plot(ladderize(treeFoley2), show.tip.label=FALSE)
#				abline(v=max(branching.times(treeFoley2))-66, col="blue",lty=1, lwd=3)
#				abline(v=max(branching.times(treeFoley2))-56, col="red",lty=1, lwd=3)
#				#axisGeo(GTS = gradstein04, unit = c("epoch","period"), cex = 0.8, ages=FALSE, gridty=3, gridcol="grey50")
#				axisPhylo( #at=max(branching.times(treeFoley2))-c(0,56,66,max(branching.times(treeFoley2))), labels=TRUE,
#							cex.axis=1.1, pos=-20)#, mgp=c(0,0.2,0))
#
#			dev.off()

