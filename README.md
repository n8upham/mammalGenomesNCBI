# Perspective in _Science_ (2023): "Genomics expands the mammalverse"
## Nathan S. Upham and Michael J. Landis
### Code for parsing NCBI mammal genome metadata

### Download from NCBI Genome on 9 Feb 2023
Query: https://www.ncbi.nlm.nih.gov/data-hub/genome/?taxon=40674

* 2800 total genomes (= 264 RefSeq + 2,536 other assemblies)
* 754 unique taxa (includes species, subspecies, hybrids, no epithet)
	- 168 chromoLevel unique taxa
	- 489 scaffoldLevel unique taxa
	- 254 contigLevel unique taxa

* 675 species after aligning with the mammal tree taxonomy and collapsing redundancy
	- 76 chromoLevel species (= 68 + 8 more from the 'redundantOrNotUsing' category)
	- 414 scaffoldLevel species (= 418 - 4 from above)
	- 185 contigLevel species (= 189 - 4 from above)

### Summary of taxonomic sampling

               Var1 	GENOME SP		TOTAL SP
1      AFROSORICIDA 	   3			  55
2         CARNIVORA 	  74			 298
3   CETARTIODACTYLA 	 128			 348
4        CHIROPTERA 	  49			1287
5         CINGULATA 	   3			  21
6    DASYUROMORPHIA 	  63			  78
7        DERMOPTERA 	   2			   2
8   DIDELPHIMORPHIA 	   3			 106
9     DIPROTODONTIA 	  85			 146
10     EULIPOTYPHLA 	  12			 491
11       HYRACOIDEA 	   2			   5
12       LAGOMORPHA 	   6			  91
13    MACROSCELIDEA 	   1			  19
14   MICROBIOTHERIA 	   1			   1
15      MONOTREMATA 	   2			   5
16 NOTORYCTEMORPHIA 	   2			   2
17 PAUCITUBERCULATA 	   0 << !		   7
18  PERAMELEMORPHIA 	  14			  22
19   PERISSODACTYLA 	  10			  24
20        PHOLIDOTA 	   4			   8
21           PILOSA 	   5			  12
22         PRIMATES 	  83			 458
23      PROBOSCIDEA 	   2			   7
24         RODENTIA 	 115			2392
25       SCANDENTIA 	   2			  20
26          SIRENIA 	   3			   5
27    TUBULIDENTATA 	   1			   1



### Summary of genomes by body size (adult average body mass of sequenced species)

