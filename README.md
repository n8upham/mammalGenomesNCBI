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

### Summary of taxonomic sampling of genome species (taxonomy of mammal tree)

| Higher taxon     | GENOME SP | TOTAL SP |
| ---------------- | --------- | -------- |
| Afrotheria       | 12        | 92       |
| Euarchontoglires | 208       | 2963     |
| Laurasiatheria   | 277       | 2456     |
| Marsupialia      | 168 **      | 362      | 
| Monotremata      | 2         | 5        |
| Xenarthra        | 8         | 33       |

** There was only 7 marsupial genomes before 2021; then 11 + 148 + 2 added in 2021 (= 161 genomes) for the current total of 168

| Order            | GENOME SP | TOTAL SP |
| ---------------- | --------- | -------- |
| AFROSORICIDA     | 3         | 55       |
| CARNIVORA        | 74        | 298      |
| CETARTIODACTYLA  | 128       | 348      |
| CHIROPTERA       | 49        | 1287     |
| CINGULATA        | 3         | 21       |
| DASYUROMORPHIA   | 63        | 78       |
| DERMOPTERA       | 2         | 2        |
| DIDELPHIMORPHIA  | 3         | 106      |
| DIPROTODONTIA    | 85        | 146      |
| EULIPOTYPHLA     | 12        | 491      |
| HYRACOIDEA       | 2         | 5        |
| LAGOMORPHA       | 6         | 91       |
| MACROSCELIDEA    | 1         | 19       |
| MICROBIOTHERIA   | 1         | 1        |
| MONOTREMATA      | 2         | 5        |
| OTORYCTEMORPHIA  | 2         | 2        |
| PAUCITUBERCULATA | 0 **      | 7        |
| PERAMELEMORPHIA  | 14        | 22       |
| PERISSODACTYLA   | 10        | 24       |
| PHOLIDOTA        | 4         | 8        |
| PILOSA           | 5         | 12       |
| PRIMATES         | 83        | 458      |
| PROBOSCIDEA      | 2         | 7        |
| RODENTIA         | 115       | 2392     |
| SCANDENTIA       | 2         | 20       |
| SIRENIA          | 3         | 5        |
| TUBULIDENTATA    | 1         | 1        |

** Paucituberculata (extant family Caenolestidae) is the only mammal order yet without a genome!

### Summary of genomes by body size (adult average body mass of sequenced species)

