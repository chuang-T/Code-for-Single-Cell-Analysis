cellranger mkgtf Sus_scrofa.Sscrofa11.1.106.gtf Sus_scrofa.Sscrofa11.1.106_PCG.gtf \
	--attribute=gene_biotype:protein_coding 

cellranger mkref \
	--genome=pig_106 \
	--fasta=Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
	--genes=Sus_scrofa.Sscrofa11.1.106_PCG.gtf