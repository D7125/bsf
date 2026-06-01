# bsf
a short demo about bsf genome align 
BSF Population-Genomics Mini-Pipeline (FASTQ → VCF → PopGen)
Mục tiêu: làm một pipeline nhỏ bám đúng JD của Entobel — QC → alignment → variant calling → variant filtering → analysis-ready genotypes → heterozygosity / FIS / kinship / PCA / ROH — chạy trên dữ liệu Black Soldier Fly (Hermetia illucens) công khai thật. Làm xong push GitHub kèm README sạch là vừa lấp gap NGS, vừa có thêm điểm "reproducible bioinformatics" + "curate public BSF datasets" mà JD yêu cầu.

Thời lượng: 1 cuối tuần nếu bạn giới hạn xuống 1 nhiễm sắc thể + vài mẫu + subsample reads. Genome BSF ~1 Gb nên đừng align cả genome với cả chục mẫu trên laptop.


Bản đồ project ↔ JD (để bạn biết mình đang tick gì)
Bước trong projectGạch đầu dòng trong JDsamplesheet.csvinput validation & sample sheet checkingFastQC + fastp + MultiQCsequencing data QCbwa-mem2 + samtoolsalignment to the reference genomebcftools call + filtervariant calling and variant filteringbiallelic SNP VCFgeneration of analysis-ready genotype filesplink2 --hetheterozygosity, FISplink2 --make-king-tablerelatedness / kinshipplink2 --pcaPCAbcftools rohROH-based inbreedingREADME + notebook báo cáoreport template + documentation(optional) wrap Nextflowmodular Nextflow pipeline

Bước 0 — Môi trường (conda/mamba)
bash# Cài mambaforge trước, rồi:
mamba create -n bsf -c bioconda -c conda-forge \
  sra-tools fastqc fastp multiqc bwa-mem2 samtools bcftools plink2 vcftools seqtk
mamba activate bsf
Tạo environment.yml (xuất từ env này) và commit cùng repo — đây chính là điểm reproducibility:
bashmamba env export --no-builds > environment.yml

Bước 1 — Reference genome (RefSeq BSF, thật)
Genome chính thức: GCF_905115235.1 iHerIll2.2.curated.20191125 (Generalovic et al. 2021, Sanger).
bashmkdir -p ref && cd ref
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/905/115/235/GCF_905115235.1_iHerIll2.2.curated.20191125/GCF_905115235.1_iHerIll2.2.curated.20191125_genomic.fna.gz
gunzip *_genomic.fna.gz
mv *_genomic.fna bsf_ref.fa
samtools faidx bsf_ref.fa

# Xem tên các chromosome/scaffold:
cut -f1,2 bsf_ref.fa.fai | sort -k2,2nr | head

# ĐỂ NHẸ MÁY: lấy 1 nhiễm sắc thể lớn (thay <CHROM> bằng tên thật ở dòng trên)
samtools faidx bsf_ref.fa <CHROM> > chr1.fa
samtools faidx chr1.fa
bwa-mem2 index chr1.fa      # nếu RAM khỏe & có thời gian thì index cả bsf_ref.fa
cd ..

bwa-mem2 index ngốn RAM theo cỡ genome. Máy yếu → dùng bwa index/bwa mem thay cho bwa-mem2, hoặc bắt buộc giới hạn 1 chromosome.


Bước 2 — Tìm & tải mẫu BSF công khai (đây là kỹ năng "curate public datasets" trong JD)
BioProject gợi ý: PRJNA1256126 — 168 genome BSF độ phủ cao, đúng bài population genomics (panel haplotype + imputation). Lý tưởng cho PCA vì có cấu trúc quần thể.

Mở NCBI SRA Run Selector, gõ PRJNA1256126.
Chọn ~6 run SRR…, ưu tiên các mẫu từ population/dòng khác nhau để PCA tách nhóm đẹp.
Tải:

bashmkdir -p raw && cd raw
for SRR in SRRxxxxxx1 SRRxxxxxx2 SRRxxxxxx3 SRRxxxxxx4 SRRxxxxxx5 SRRxxxxxx6; do
  prefetch $SRR
  fasterq-dump $SRR --split-files -O .
  gzip ${SRR}_1.fastq ${SRR}_2.fastq
done
cd ..

# Subsample cho nhẹ (vd ~2 triệu cặp read/mẫu):
for SRR in SRRxxxxxx1 ... ; do
  seqtk sample -s100 raw/${SRR}_1.fastq.gz 2000000 | gzip > raw/${SRR}_1.sub.fq.gz
  seqtk sample -s100 raw/${SRR}_2.fastq.gz 2000000 | gzip > raw/${SRR}_2.sub.fq.gz
done
Tạo samplesheet.csv (đúng tinh thần "input validation / sample sheet checking"):
csvsample,fastq_1,fastq_2
SAMPLE1,raw/SRRxxxxxx1_1.sub.fq.gz,raw/SRRxxxxxx1_2.sub.fq.gz
SAMPLE2,raw/SRRxxxxxx2_1.sub.fq.gz,raw/SRRxxxxxx2_2.sub.fq.gz
...
Viết một script Python 10 dòng đọc CSV này, kiểm tra file tồn tại + cột hợp lệ → đó là "input validation".

Plan B (đảm bảo có kết quả nếu hết thời gian): một nghiên cứu khác công bố sẵn VCF quần thể BSF đã filter tại https://genetics.ghpc.au.dk/zexi/BSF/dataset/BSF_filDP.vcf.gz. Nếu align/calling quá nặng, tải VCF này về và nhảy thẳng tới Bước 6 (PCA/het/kinship/ROH). Vẫn là kết quả thật trên BSF.


Bước 3 — QC (sequencing data QC)
bashmkdir -p qc trim
fastqc raw/*.sub.fq.gz -o qc/

for SRR in SRRxxxxxx1 ... ; do
  fastp \
    -i raw/${SRR}_1.sub.fq.gz -I raw/${SRR}_2.sub.fq.gz \
    -o trim/${SRR}_1.fq.gz   -O trim/${SRR}_2.fq.gz \
    -j qc/${SRR}.fastp.json  -h qc/${SRR}.fastp.html
done

multiqc qc/ -o qc/        # gộp toàn bộ QC thành 1 report

Bước 4 — Alignment + sort + markdup (alignment to reference)
bashmkdir -p bam
REF=ref/chr1.fa     # hoặc ref/bsf_ref.fa

for SRR in SRRxxxxxx1 ... ; do
  bwa-mem2 mem -t 4 \
    -R "@RG\tID:${SRR}\tSM:${SRR}\tPL:ILLUMINA" \
    $REF trim/${SRR}_1.fq.gz trim/${SRR}_2.fq.gz \
  | samtools sort -@4 -o bam/${SRR}.sorted.bam
  # đánh dấu PCR duplicate
  samtools collate -@4 -O bam/${SRR}.sorted.bam \
  | samtools fixmate -m - - \
  | samtools sort -@4 - \
  | samtools markdup -@4 - bam/${SRR}.dedup.bam
  samtools index bam/${SRR}.dedup.bam
  samtools flagstat bam/${SRR}.dedup.bam > qc/${SRR}.flagstat.txt   # QC alignment
done

ls bam/*.dedup.bam > bamlist.txt

Bước 5 — Variant calling + filtering → analysis-ready genotypes
bashmkdir -p vcf
REF=ref/chr1.fa

# Joint calling tất cả mẫu cùng lúc
bcftools mpileup -f $REF -a AD,DP,SP -b bamlist.txt -Ou \
| bcftools call -mv -Oz -o vcf/raw.vcf.gz
bcftools index vcf/raw.vcf.gz

# Filter cơ bản
bcftools filter -e 'QUAL<30 || INFO/DP<10' vcf/raw.vcf.gz -Oz -o vcf/filt.vcf.gz

# Giữ SNP biallelic → "analysis-ready genotype files"
bcftools view -m2 -M2 -v snps vcf/filt.vcf.gz -Oz -o vcf/snps.vcf.gz
bcftools index vcf/snps.vcf.gz

bcftools stats vcf/snps.vcf.gz > vcf/snps.stats.txt

Bước 6 — Population genetics (đánh đúng từng bullet downstream của JD)
bashmkdir -p pop
# Nạp vào plink2 (--allow-extra-chr vì tên scaffold không chuẩn 1..22)
plink2 --vcf vcf/snps.vcf.gz --allow-extra-chr \
       --set-all-var-ids @:#:\$r:\$a --make-pgen --out pop/bsf
PCA
bashplink2 --pfile pop/bsf --allow-extra-chr --pca 10 --out pop/bsf_pca
# -> pop/bsf_pca.eigenvec : vẽ PC1 vs PC2 bằng R/Python, tô màu theo population
Heterozygosity + FIS (hệ số cận huyết F per mẫu)
bashplink2 --pfile pop/bsf --allow-extra-chr --het --out pop/bsf_het
# cột F = inbreeding coefficient; O(HOM)/E(HOM) cho heterozygosity quan sát vs kỳ vọng
Relatedness / kinship
bashplink2 --pfile pop/bsf --allow-extra-chr --make-king-table --out pop/bsf_king
# -> ma trận kinship KING giữa các cặp mẫu
ROH-based inbreeding
bashbcftools roh -G30 --AF-dflt 0.4 vcf/snps.vcf.gz > pop/bsf.roh.txt
# tổng chiều dài ROH / chiều dài genome ~ F_ROH (chỉ số cận huyết dựa ROH)
# hoặc: plink2 ... --homozyg
