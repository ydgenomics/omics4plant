# editor: yangdong
# image: txdb-bsgenome
# 260427
# 在运行前请确定要构建的.fa的染色体都是要被map到的染色体和染色体
# workflow: 1)build bsgenome object; 2)library package; 3)build gene annotation object

# FA="/data/input/Files/User/yinzhanhao/index/rice/osa1_r7.asm.chrs.fa"
# gtf="/data/input/Files/User/yinzhanhao/index/rice/osa1_r7.all_models.gtf"
FA=${1:-"/data/input/Files/User/yinzhanhao/index/rice/osa1_r7.asm.chrs.fa"}
GTF=${2:-"/data/input/Files/User/yinzhanhao/index/rice/osa1_r7.all_models.gtf"}

echo "[INFO] Splitting fasta file by chromosome"
echo "chromosome names in fasta file: $(grep '^>' "$FA" | cut -d'>' -f2)"
mkdir -p split
# awk -v pattern="^>${STR_CHROMOSOME}" '$0 ~ pattern {if(OUT) close(OUT); OUT="split/" substr($0,2) ".fa"}; OUT {print > OUT}' $FA
mkdir -p split && awk '/^>/ {if(OUT)close(OUT); OUT="split/"substr($1,2)".fa"} OUT{print > OUT}' "$FA"

echo "[INFO] Configuring seed.txt for BSgenomeForge"
printf "Package: BSgenome.species\n" > seed.txt
printf "Title: Full genome sequences for species\n" >> seed.txt
printf "Description: Full genome sequences for species.\n" >> seed.txt
printf "Version: 1.0.0\n" >> seed.txt
printf "organism: Genus species\n" >> seed.txt
printf "common_name: Genus\n" >> seed.txt
printf "provider: website\n" >> seed.txt
printf "provider_version:\n" >> seed.txt
printf "release_date:\n" >> seed.txt
printf "release_name:\n" >> seed.txt
printf "source_url:\n" >> seed.txt
printf "organism_biocview: Genus_species\n" >> seed.txt
printf "BSgenomeObjname: Genus_species\n" >> seed.txt
printf "circ_seqs:\n" >> seed.txt
CHRS=$(grep '^>' "$FA" | cut -d'>' -f2 | awk '{printf "\"%s\", ", $1}' | sed 's/, $//')
printf "seqnames: c($CHRS)\n" >> seed.txt
# printf "seqnames: paste("Chr", c(1:12), sep="")\n" >> seed.txt
printf "SrcDataFiles:\n" >> seed.txt
printf "PkgExamples:\n" >> seed.txt
printf "seqs_srcdir: ./split\n" >> seed.txt
printf "ondisk_seq_format: rda\n" >> seed.txt

echo "[INFO] Building BSgenome package"
Rscript -e "library(BSgenomeForge); BSgenomeForge::forgeBSgenomeDataPkg('seed.txt', replace=TRUE)"
# sed -i '7c\Depends: R (>= 4.2.0), GenomeInfoDb (>= 1.34.9), BSgenome (>= 1.74.0)' ./BSgenome.species/DESCRIPTION
sed -i '7c\Depends: R, GenomeInfoDb, BSgenome' ./BSgenome.species/DESCRIPTION
R CMD build BSgenome.species
R CMD check --no-manual BSgenome.species_1.0.0.tar.gz
# R CMD INSTALL BSgenome.species_1.0.0.tar.gz

echo "[INFO] Build gene annotation object"
echo "Gtf info:"
echo "总行数:" && wc -l $GTF
echo -e "\n特征类型分布:" && awk '!/^#/ {print $3}' $GTF | sort | uniq -c | sort -rn
echo -e "\n基因数量估计:" && grep -c 'gene_id' $GTF | head -10
echo -e "\n转录本数量估计:" && grep -c 'transcript_id' $GTF | head -5
echo -e "\n染色体列表:" && cut -f1 $GTF | sort -u