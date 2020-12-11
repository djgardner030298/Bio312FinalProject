# Bio312FinalProject


ncbi-acc-download -F fasta -m protein XP_032224579.1

# Downloads the query protein sequence.

blastp -db ~/data/blast/allprotein.fas -query XP_032224579.1.fa -outfmt 0 -max_hsps 1 > XP_032224579.1.blastp.typical.out

# Protein query used to execute a blast search.

blastp -db ~/data/blast/allprotein.fas -query XP_032224579.1.fa -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle" -max_hsps 1 -out XP_032224579.1.blastp.detail.out

# Easier and more detailes blast search using protein query sequence.

awk '{if ($6<0.00000000000001)print $1 }' XP_032224579.1.blastp.detail.out > XP_032224579.1.blastp.detail.filtered.out

# E-value less than 1e-14, which filters the output file to fulfill requirement

wc -l XP_032224579.1.blastp.detail.filtered.out

# Allows us to count the results of the filtered output file.

seqkit grep --pattern-file XP_032224579.1.blastp.detail.filtered.out ~/data/blast/allprotein.fas > XP_032224579.1.blastp.detail.filtered.fas

# Allows access to the downloaded proteomes.

muscle -in XP_032224579.1.blastp.detail.filtered.fas -out XP_032224579.1.blastp.detail.filtered.aligned.fas

# Aligns protein sequences for the genes globally. FASTA file serve as input and output.

t_coffee -other_pg seq_reformat -in XP_032224579.1.blastp.detail.filtered.aligned.fas -output sim

# Method for fast and accurate multiple sequence alignment, was used to provide statistics on the protein alignment, allowing for the analyzation of the MUSCLE input alignments based on percent identity.

t_coffee -other_pg seq_reformat -in XP_032224579.1.blastp.detail.filtered.aligned.fas -action +rm_gap 50 -out XP_032224579.1.allhomologs.aligned.r50.fa 

# Removing highly gapped position (remove any coloumns less than 50% residues)

awk '{if ($3>=50 && $6<0.0000000001)print $1 }' ~/XP_032224579.1.blastp.detail.out > ~/XP_032224579.1.blastp.detail.filtered.out

# Filtering blast results by length

iqtree -s iqtree -s XP_032224579.1.allhomologs.aligned.r50.ren.fa -bb 1000 -nt 2 --prefix X -bb 1000 -nt 2 --prefix XP_032224579.1.r50.ufboot.treefile

# Allows for bootstrap support values. The -bb option executes a bootstrap analysis by resmapling the alignment 1000 times with replacemnt, and calculating the maximum likelihood tree for each psuedo-replicate.

gotree reroot midpoint -i XP_032224579.1.r50.ufboot.treefile -o XP_032224579.1.r50.ufboot.MidRoot.treefile

# Rerooted the tree.

nw_display XP_032224579.1.r50.ufboot.MidRoot.treefile

# Allows the visualization of bootstrap support values.

nw_display -s XP_032224579.1.r50.ufboot.MidRoot.treefile -w 1000 -b 'opacity:0' > XP_032224579.1.r50.ufboot.MidRoot.svg

# Picture produced to visualize bootstrap support values.

java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -b myproteinbatch.txt --reconcile --speciestag prefix --savepng --treestats --events --phylogenomics

# Notung reconciles the gene and species tree.

python ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g XP_032224579.1.r50.ufboot.MidRoot.treefile.reconciled --include.species

# RecPhyloXML object is produced. 

java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -b myproteinbatch.txt --reconcile --speciestag prefix --savepng --treestats --events --phylogenomics

# Command executes the reconciliation

java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -b myproteinbatch.txt --rearrange --speciestag prefix  --savepng --treestats --events  --outputdir mygenereconcilerearrange --edgeweights name --threshold 90

# Reconciling and rearranging tree. Threshold of 90 is used.

cp ~/labs/lab5/lab5-djgardner030298/XP_032224579.1.blastp.detail.filtered.fas .

# Copying protein sequences oobtained earlier.

nano XP_032224579.1.blastp.detail.filtered.ren.fas

# Renaming files to match previous ones. 

iprscan5 --email your.email@stonybrook.edu --multifasta --useSeqId --sequence XP_032224579.1.blastp.detail.filtered.ren.fas

# Running interproscan to analuze sequences. 

cat *.tsv.txt > XP_032224579.1.domains.all.tsv

# Concentrating all of the domain generated files for each gene.

grep Pfam XP_032224579.1.domains.all.tsv > yourgenefamily.domains.pfam.tsv

# Filters domains to those only found in Pfam databse.

awk 'BEGIN{FS="\t"} {print $1"\t"$3"\t"$7"@"$8"@"$5}' yourgenefamily.domains.pfam.tsv | datamash -sW --group=1,2 collapse 3 | sed 's/,/\t/g' | sed 's/@/,/g' > yourgenefamily.domains.pfam.evol.tsv

# Rearranges interproscan output.
