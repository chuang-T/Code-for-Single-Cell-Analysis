cellranger count --id=210614A-SPF-MLN \
--fastqs=${my_dir} \  #The address of fq file
--sample=210614A-SPF-MLN \
--transcriptome=pig_106 \
--no-bam \
--nosecondary
echo '210614A-SPF-MLN Finished'