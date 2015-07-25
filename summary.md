# Medicago truncatula nodulation gene network

## 1. some known nodulation genes

[start.fa](./start.fa)

[reference](http://www.nature.com/ncomms/journal/v1/n1/fig_tab/ncomms1009_F6.html)

1. Nod factor perception

  **NFP**
  
  **LYK3**
  
  LYK4

2. Signal transduction
  
  **NORK** (SymRK)
  
  **DMI1** (pollux)

3. Ca signal interpretation

  **DMI3** (CCaMK, recognizes the calcium spiking)
  
  **CYCLOPS** (IPD3)

4. Cortical cytokinin signaling
  
 CRE1 (LHK1) no seq
 

5.  Transcriptional activation/regulation

  **NSP1**
  
  **NSP2**
  
  **NIN**

6. Root hair curling

  **NAP1**
  
  LIN (Cerberus)

7. nodule number

  **SUNN**
  
  **EIN2 (SKL1)**
  

## 2. mapping

Download microarray fasta file from [Affymetrix](http://www.affymetrix.com/catalog/131472/AFFY/Medicago+Genome+Array#1_1)

```bash

# correct fasta ids
sed -i 's/>consensus:Medicago:/>/' Medicago.consensus.fa
sed -i 's/;/ /' Medicago.consensus.fa

# blast
makeblastdb -in Medicago.consensus.fa -dbtype nucl -title Medicago.consensus -out Medicago.consensus
tblastn -num_threads 10 -max_target_seqs 5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' -query start.fa -db Medicago.consensus >start.blast.tab

cp start.blast.tab start.blast.edit.tab
# edit start.blast.edit.tab manually to keep the 'true' mappings.
cut -f2 start.blast.edit.tab >start.lst.txt
```

## 2. Get co-expressed probesets

download expression data from [mtgea](http://mtgea.noble.org/v3/experiments.php) (select 'All Means' and 'Mtr: Medicago truncatula only')

- Number of Experiments: 277
- Number of GeneChips: 739

Sample list

```bash
grep -oh 'Mean (.*)' samples.txt | sed 's/.*(//' | sed 's/)//' > samples.lst
wc -l samples.lst
# 274
```

It's difficult to download all expression data once, so I had to download a few samples per time, and then combine all these samples together.

[combine.R](combine.R)

```
Mtr.6956.1.S1_at and Mtr.51192.1.S1_at does not vary much, and are removed.
```

```R


```

```bash
grep -E '^>' Mt4.0v1_GenesProteinSeq_20130731_1800.fasta >Mt4.0v1_GenesProteinSeq_20130731_1800.annot.txt
sed -i 's/ | /\t/g' Mt4.0v1_GenesProteinSeq_20130731_1800.annot.txt
sed -i 's/>//' Mt4.0v1_GenesProteinSeq_20130731_1800.annot.txt
awk -F '\t' -v OFS='\t' '{print $1,$2}' Mt4.0v1_GenesProteinSeq_20130731_1800.annot.txt >Mt4.0v1_GenesProteinSeq_20130731_1800.annot.tab

makeblastdb -title mt4 -out mt4 -dbtype 'prot' -in Mt4.0v1_GenesProteinSeq_20130731_1800.fasta 
blastp -max_target_seqs 5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' -query start.fa -db mt4 -num_threads 10 >start.mt4.tab

makeblastdb -in uniprot.fa -dbtype prot -title uniprot -out uniprot
blastp -max_target_seqs 1 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' -query Mt4.0v1_GenesProteinSeq_20130731_1800.fasta -db uniprot -num_threads 10 >Mt4.0v1_GenesProteinSeq_20130731_1800.uniprot.blast.tab

awk -v OFS='\t' '{if ($3>95 && $4/$13>0.9) print $0}' Mt4.0v1_GenesProteinSeq_20130731_1800.uniprot.blast.tab >Mt4.0v1_GenesProteinSeq_20130731_1800.uniprot.select.tab


grep -E '^>' uniprot.fa >uniprot.head
sed -i 's/>//' uniprot.head
sed -i 's/ OS=.*//' uniprot.head
sed -i 's/ /\t/' uniprot.head
sed -i 's/sp|//' uniprot.head
sed -i 's/tr|//' uniprot.head

grep -E '^>.*GN=' uniprot.fa >uniprot.head2
sed -i 's/>//' uniprot.head2
sed -i 's/MEDTR .*GN=/MEDTR\t/' uniprot.head2
sed -i 's/ PE=.*//' uniprot.head2

wget ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz
gunzip gene_association.goa_uniprot.gz
perl go.pl gene.mod.lst gene_association.goa_uniprot >gene.mod.go

cut -f1,3 gene.all.go >gene.all.goC
sed -i 's/ /, /g' gene.all.goC
cut -f1,5 gene.all.go >gene.all.goF
sed -i 's/ /, /g' gene.all.goF
cut -f1,7 gene.all.go >gene.all.goP
sed -i 's/ /, /g' gene.all.goP
```


