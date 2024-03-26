## Filter non-mito contigs from assemblies ##
for dir in F*
do
wd=$(pwd)
cd "$dir"/cox1_readmapping/metaspades_assembly
mkdir mito-contigs
grep mitochon contigs.fasta.vs.nt.cul5.1e5.megablast.out | cut -f1 | sort | uniq >mito-contigs/mitocontigids.list
cd "$wd"
done

for dir in F*
do
wd=$(pwd)
cd "$dir"/cox1_readmapping/metaspades_assembly/mito-contigs
while read -r line
do grep -w -A 1 "$line" ../contigs.fasta
done<mitocontigids.list >mito_contigs.fasta
cd "$wd"
done


## Filter mito contigs that are smaller than smallest possible genes (excluding ATP8) ##
for dir in F*
do while read -r line
do length=$(awk -F"_" '{print $4}' <(echo "$line"))
nodeid=$(awk -F"_" '{print $2}' <(echo "$line"))
if [[ "$length" -gt 330 ]]
then grep -w -A 1 "$line" "$dir"/*mito_contigs.fasta >"$dir"/NODE_"$nodeid".fasta
fi
done<"$dir"/nodeids.list
done

## Formatting for Community Analysis ##
for dir in F*
do cd "$dir"
cat NODE* >"$dir"_lengthfilt_mito_contigs.fasta
cd ..
done

## BLAST filtered contigs against custom mitogene database ##
for dir in F*
do cd "$dir"
blastn -query "$dir"_lengthfilt_mito_contigs.fasta -subject ~/Projects/custom_mitogene_db/refseq_aug31_2022/custom_db_working_version.fasta -outfmt "6 qseqid sseqid pident qcovs length qstart qend sstart send mismatch bitscore evalue" -out "$dir"_blastn_results.tsv
cd ..
done

## Extract community info from blast tables ##
# 92 % ID & 32% QCOV 100bp length cutoffs (non-annotated data)
for dir in F*
do cd "$dir"
echo "$dir"
awk -F"\t" '$3 >= 92 {print $0}' "$dir"_blastn_results.tsv | awk -F"\t" '$4 >= 32 {print $0}' | awk -F"\t" '$5 >= 100 {print $0}' | cut -f2 | awk -F"|" '{print $(NF-1), $NF}' | sort | uniq >output.txt
echo "$dir" >sampleid.txt
#cat sampleid.txt output.txt >"$dir"_sampletaxa_98ID32QCOV.list
mv output.txt "$dir"_sampletaxa_92ID32QCOV100length.list
cd ..
done

## Extract top blast hits from blast tables ##
for dir in F*
do
wd=$(pwd)
cd "$dir"/cox1_readmapping/metaspades_assembly/mito-contigs/
ids=$(cut -f1 blastn_results.tsv | uniq)
while read -r line
do grep -w "$line" blastn_results.tsv | head -n1 >tophits_blastn_results.tsv
done< <(echo "$ids")
cd "$wd"
done

# 92 % ID & 32% QCOV 100bp length cutoffs (non-annotated data)
for dir in F*
do 
wd=$(pwd)
cd "$dir"/cox1_readmapping/metaspades_assembly/mito-contigs/
echo "$dir"
awk -F"\t" '$3 >= 91 {print $0}' tophits_blastn_results.tsv | awk -F"\t" '$4 >= 32 {print $0}' | awk -F"\t" '$5 >= 100 {print $0}' >output.txt
#echo "$dir" >sampleid.txt
#cat sampleid.txt output.txt >"$dir"_sampletaxa_98ID32QCOV.list
nodeids=$(cut -f1 output.txt | sort | uniq)

while read -r line
do
grep -w "$line" output.txt
done < <(echo "$nodeids") >92ID32QCOV100length_blastn_results_tophits.tsv

rm output.txt
cd "$wd"
done


# build comm presence/absence table using selected cutoffs and only for COX1
for dir in F*
do
wd=$(pwd)
cd "$dir"/cox1_readmapping/metaspades_assembly/mito-contigs/
echo "$dir"

output=$(cat 92ID32QCOV100length_blastn_results_tophits.tsv)
idslist=$(cut -f1 92ID32QCOV100length_blastn_results_tophits.tsv)
out=$(while read -r line
do
grep -w "$line" <(echo "$output") | cut -f2
done< <(echo "$idslist")) 
echo "$out" >92ID32QCOV100length.list

outputsize=$(wc -l <(echo "$out") | cut -f1 -d" ")
if [ $outputsize -gt 0 ]
then
for i in $( seq 1 $outputsize )
do  echo "1"
done
for i in {1..50}
do
echo "0"
done 
else
echo "0" 
fi >detectoutput.txt

echo "$dir" >dir.txt
echo "detected" >fill.txt
echo "$out"
paste dir.txt fill.txt >header.txt
paste 92ID32QCOV100length.list detectoutput.txt >tmp1
cat header.txt tmp1 >92ID32QCOV100length.tsv
cd "$wd"
done
paste F*/cox1_readmapping/metaspades_assembly/mito-contigs/92ID32QCOV100length.tsv >outputtable_92ID32QCOV100length.tsv


# extract cox1 contig tophits passing 92% id and 32% qcov cutoffs
for dir in F*
do
cd "$dir"
echo "$dir"
output=$(awk -F"\t" '$3 >= 92 {print $0}' "$dir"_blastn_results.tsv | awk -F"\t" '$4 >= 32 {print $0}' | awk -F"\t" '$5 >= 100 {print $0}')
idslist=$(cut -f1 <(echo "$output") | uniq)
out=$(while read -r line
do
grep COX1 <(echo "$output") | grep -w "$line" | head -n1
done< <(echo "$idslist")) 
(echo "$dir"
echo "$out") >"$dir"_92ID32QCOV_blastn_results.tsv
cd ..
done


# Alternative table builder
for file in F*
do
sampleid=$(awk -F"_" '{print $1}' <(echo $file))
echo "$sampleid"
while read -r line
do
result=$(grep "$line" "$file")
if [ -z "$result" ]
then
	echo "$line 0"
else
	echo "$line 1"
fi
done<recovered_taxa_list.list >"$sampleid"_presenceabsence_tmp
done


# Rename sequence IDs
for dir in F*
do 
wd=$(pwd)
cd "$dir"

sed -i "s/NODE/'$dir'_NODE/" comm_cox1_seqs.fasta
sed -i "s/'//g" comm_cox1_seqs.fasta

seqidtotaxaid=$(while read -r line
do
grep -w "$line" "$dir"_annotated_blast_results_derep_bitscoresorted.tsv | cut -f1,2 | awk -F"|" '{print $1,$(NF-1),$NF}' | sed 's/ /_/g'
done<cox1seqids.list)

while read -r line2
do
toreplace=$(cut -f1 <(echo "$line2"))
replacewith=$(cut -f2 <(echo "$line2"))
echo "$toreplace"
echo "$replacewith"
sed -i "s/$toreplace.*/$replacewith/" comm_cox1_seqs.fasta
done< <(echo "$seqidtotaxaid")

cd "$wd"
done
cat F*/comm_cox1_seqs.fasta >agrifield_cox1_seqs.fasta
 

# select top open reading frames for each node using bitscore column from blastn results
for dir in F*
do
idlist=$(cut -f1 "$dir"/orfipy_*/blastn_results.tsv | awk -F"_" '{print $1,$2,$3}' | sed 's/ /_/g' | sort | uniq)
while read -r line
do
grep "$line" "$dir"/orfipy_"$dir"*/blastn_results.tsv | sort -k12 | cut -f1
done< <(echo "$idlist") >"$dir"/orfipy_"$dir"_lengthfilt_mito_contigs.fasta_out/bitscore_determined_orfs_nodeids.list
done 


# select top open reading frames for each node using bitscore column from blastn results
for dir in F*
do
idlist=$(cut -f1 "$dir"/orfipy_*/blastn_results.tsv | awk -F"_" '{print $1,$2,$3}' | sed 's/ /_/g' | sort | uniq)
while read -r line
do
grep "$line" "$dir"/orfipy_"$dir"*/blastn_results.tsv | sort -k12
done< <(echo "$idlist") >"$dir"/orfipy_"$dir"_lengthfilt_mito_contigs.fasta_out/bitscore_determined_orfs_blastn_results.tsv
done 

# Updated sequence ID renaming script
# includes node ID, accession ID, position matched to reference, taxonomy
for dir in F*
do
wd=$(pwd)
cd "$dir"
# Pull fasta of ORFs for each node for each sample
cp ../../"$dir"/orfipy_"$dir"_lengthfilt_mito_contigs.fasta_out/output-dna.fasta ./
# format files
sed -i 's/\[/\(/g' output-dna.fasta
sed -i 's/\[/\(/g' "$dir"_annotated_blast_results_derep_bitscoresorted.tsv
# extract relevant ORFs for each node
grep --no-group-separator -A 1 -w -f cox1seqids.list output-dna.fasta >cox1_seqs.fasta
# loop through each relevant sequence, extract, and modify sequence header
while read -r line
do
currentseqid=$(grep -w "$line" cox1_seqs.fasta | sed 's/>//')
blastnout=$(grep -w "$line" "$dir"_annotated_blast_results_derep_bitscoresorted.tsv)
node=$(cut -f1 -d" " <(echo "$currentseqid") | awk -F"_" '{print $1,$2}' | sed 's/ /_/g')
length=$(cut -f4 -d" " <(echo "$currentseqid") | cut -f2 -d":")
refid=$(cut -f2 <(echo "$blastnout") | awk -F"|" '{print $1}')
reftax=$(cut -f2 <(echo "$blastnout") | awk -F"|" '{print $(NF-1),$NF}' | sed 's/ /_/g')
refstartpos=$(cut -f10 <(echo "$blastnout"))
refstoppos=$(cut -f11 <(echo "$blastnout"))
replacewith=$(echo "$dir" "$node" "length_""$length" "$refid" "$reftax" "start""$refstartpos" "stop""$refstoppos" | sed 's/ /_/g' | sed 's/>//')
echo "$currentseqid"
echo "$replacewith"
sed -ier "s/$currentseqid/$replacewith/g" cox1_seqs.fasta
done<cox1seqids.list 
cd "$wd"
done

