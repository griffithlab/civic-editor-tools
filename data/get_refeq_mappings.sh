wget ftp://ftp.ncbi.nih.gov/gene/DATA/gene2refseq.gz
gzcat gene2refseq.gz | head -n 1 > header.tsv
gzcat gene2refseq.gz | grep '^9606\t' > gene2refseq_human.tmp.tsv
cat header.tsv gene2refseq_human.tmp.tsv > gene2refseq_human.tsv
rm -f header.tsv gene2refseq_human.tmp.tsv gene2refseq.gz
gzip gene2refseq_human.tsv
