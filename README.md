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

