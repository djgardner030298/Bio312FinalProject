# Bio312FinalProject


ncbi-acc-download -F fasta -m protein XP_032224579.1 #Downloading protein sequence

blastp -db ~/data/blast/allprotein.fas -query XP_032224579.1.fa -outfmt 0 -max_hsps 1 > XP_032224579.1.blastp.typical.out

blastp -db ~/data/blast/allprotein.fas -query XP_032224579.1.fa -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle" -max_hsps 1 -out XP_032224579.1.blastp.detail.out

awk '{if ($6<0.00000000000001)print $1 }' XP_032224579.1.blastp.detail.out > XP_032224579.1.blastp.detail.filtered.out

awk '{if ($6<0.00000000000001)print $1 }' XP_032224579.1.blastp.detail.out > XP_032224579.1.blastp.detail.filtered.out

rm XP_032224579.1.blastp.detail.filtered.out

awk '{if ($6<0.00000000000001)print $1 }' XP_032224579.1.blastp.detail.out > XP_032224579.1.blastp.detail.filtered.out

nano XP_032224579.1.blastp.detail.filtered.out

wc -l XP_032224579.1.blastp.detail.filtered.out

seqkit grep --pattern-file XP_032224579.1.blastp.detail.filtered.out ~/data/blast/allprotein.fas > XP_032224579.1.blastp.detail.filtered.fas

muscle -in XP_032224579.1.blastp.detail.filtered.fas -out XP_032224579.1.blastp.detail.filtered.aligned.fas

t_coffee -other_pg seq_reformat -in XP_032224579.1.blastp.detail.filtered.aligned.fas -output sim

t_coffee -other_pg seq_reformat -in XP_032224579.1.blastp.detail.filtered.aligned.fas -action +rm_gap 50 -out XP_032224579.1.allhomologs.aligned.r50.fa # Removing highly gapped position (remove any coloumns less than 50% residues)

alv -k XP_032224579.1.blastp.detail.filtered.aligned.fas | less -r

alv -kli XP_032224579.1.blastp.detail.filtered.aligned.fas

awk '{if ($3>=50 && $6<0.0000000001)print $1 }' ~/XP_032224579.1.blastp.detail.out > ~/XP_032224579.1.blastp.detail.filtered.out

iqtree -s iqtree -s XP_032224579.1.allhomologs.aligned.r50.ren.fa -bb 1000 -nt 2 --prefix X -bb 1000 -nt 2 --prefix XP_032224579.1.r50.ufboot.treefile

gotree reroot midpoint -i XP_032224579.1.r50.ufboot.treefile -o XP_032224579.1.r50.ufboot.MidRoot.treefile

nw_display XP_032224579.1.r50.ufboot.MidRoot.treefile

nw_display -s XP_032224579.1.r50.ufboot.MidRoot.treefile -w 1000 -b 'opacity:0' > XP_032224579.1.r50.ufboot.MidRoot.svg

java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -b myproteinbatch.txt --reconcile --speciestag prefix --savepng --treestats --events --phylogenomics

python ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g XP_032224579.1.r50.ufboot.MidRoot.treefile.reconciled --include.species

column -t -s $'\t' myproteinbatch.txt.species.tre.geneCount.txt | less -S

java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -b myproteinbatch.txt --root --speciestag prefix --savepng --treestats --events --phylogenomics

cp ~/labs/lab5/lab5-djgardner030298/XP_032224579.1.blastp.detail.filtered.fas .

nano XP_032224579.1.blastp.detail.filtered.ren.fas

iprscan5 --email your.email@stonybrook.edu --multifasta --useSeqId --sequence XP_032224579.1.blastp.detail.filtered.ren.fas

cat *.tsv.txt > XP_032224579.1.domains.all.tsv

grep Pfam XP_032224579.1.domains.all.tsv > yourgenefamily.domains.pfam.tsv

awk 'BEGIN{FS="\t"} {print $1"\t"$3"\t"$7"@"$8"@"$5}' yourgenefamily.domains.pfam.tsv | datamash -sW --group=1,2 collapse 3 | sed 's/,/\t/g' | sed 's/@/,/g' > yourgenefamily.domains.pfam.evol.tsv
