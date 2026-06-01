# BSF Population-Genomics Mini-Pipeline (FASTQ → VCF → PopGen)

**Mục tiêu: một pipeline nhỏ bám đúng JD của Entobel — QC → alignment → variant calling → variant filtering → analysis-ready genotypes → heterozygosity / FIS / kinship / PCA / ROH — chạy trên **dữ liệu Black Soldier Fly (*Hermetia illucens*) thực tế **.
---

## Bản đồ project ↔ JD 

| Bước trong project | Gạch đầu dòng trong JD |
|---|---|
| samplesheet.csv | input validation & sample sheet checking |
| FastQC + fastp + MultiQC | sequencing data QC |
| bwa-mem2 + samtools | alignment to the reference genome |
| bcftools call + filter | variant calling and variant filtering |
| biallelic SNP VCF | generation of analysis-ready genotype files |
| plink2 --het | heterozygosity, FIS |
| plink2 --make-king-table | relatedness / kinship |
| plink2 --pca | PCA |
| bcftools roh | ROH-based inbreeding |
| README + notebook báo cáo | report template + documentation |
| (optional) wrap Nextflow | modular Nextflow pipeline |

---

## Bước 0 — Môi trường (conda)

```bash
#tao moi truong
conda create -n bsf -c conda-forge -c bioconda \
  sra-tools fastqc fastp multiqc bwa-mem2 samtools bcftools plink2 vcftools seqtk
```
```bash
#tao file env.yml
conda env export --no-builds > environment.yml
```
## Bước 1 — Reference genome 

Genome chính thức: `GCF_905115235.1 iHerIll2.2.curated.20191125` (Generalovic et al. 2021, Sanger).

```bash
mkdir -p ref && cd ref
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/905/115/235/GCF_905115235.1_iHerIll2.2.curated.20191125/GCF_905115235.1_iHerIll2.2.curated.20191125_genomic.fna.gz
gunzip *_genomic.fna.gz
mv *_genomic.fna bsf_ref.fa
samtools faidx bsf_ref.fa

# Xem tên các chromosome/scaffold:
cut -f1,2 bsf_ref.fa.fai | sort -k2,2nr | head

# Lấy 1 nhiễm sắc thể lớn thay vì cả cụm
samtools faidx bsf_ref.fa <CHROM> > chr1.fa
samtools faidx chr1.fa
bwa-mem2 index chr1.fa      # nếu RAM khỏe & có thời gian thì index cả bsf_ref.fa
cd ..
```

#dùng `bwa index`/`bwa mem` thay cho bwa-mem2, hoặc bắt buộc giới hạn 1 chromosome vì lý do liên quan tới ram và thời gian chạy.

---

## Bước 2 — Tìm & tải mẫu BSF công khai 

BioProject: **PRJNA1256126** 

1. Mở **NCBI SRA Run Selector**, gõ `PRJNA1256126`.
2. Chọn ~6 run `SRR…`, ưu tiên các mẫu từ population/dòng khác nhau để PCA tách nhóm đẹp.
3. Tải:

```bash
mkdir -p raw && cd raw
for SRR in SRRxxxxxx1 SRRxxxxxx2 SRRxxxxxx3 SRRxxxxxx4 SRRxxxxxx5 SRRxxxxxx6; do
  prefetch $SRR
  fasterq-dump $SRR --split-files -O .
  gzip ${SRR}_1.fastq ${SRR}_2.fastq
done
cd ..
```

Tạo **samplesheet.csv**:
```csv
sample,fastq_1,fastq_2
SAMPLE1,raw/SRRxxxxxx1_1.sub.fq.gz,raw/SRRxxxxxx1_2.sub.fq.gz
SAMPLE2,raw/SRRxxxxxx2_1.sub.fq.gz,raw/SRRxxxxxx2_2.sub.fq.gz
...
```
 Input validation:
"""Validate the sample sheet for the BSF population-genomics pipeline.

Checks performed:
  1. Required columns are present: sample, fastq_1, fastq_2
  2. Sample names are non-empty and unique
  3. Every FASTQ path listed actually exists on disk

Usage:
  python validate_samplesheet.py samplesheet.csv

Exit code 0 = valid, 1 = at least one problem found
(so a pipeline / Nextflow can stop early on a bad sheet).

```bash
import csv
import sys
from pathlib import Path

REQUIRED = ["sample", "fastq_1", "fastq_2"]


def main(path):
    errors = []
    seen = set()

    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)

        # 1. required columns present?
        missing = [c for c in REQUIRED if c not in (reader.fieldnames or [])]
        if missing:
            sys.exit(f"ERROR: missing column(s): {', '.join(missing)}")

        # 2 & 3. check every row (line 2 = first data row, after the header)
        for i, row in enumerate(reader, start=2):
            name = (row["sample"] or "").strip()
            if not name:
                errors.append(f"line {i}: empty sample name")
            elif name in seen:
                errors.append(f"line {i}: duplicate sample '{name}'")
            else:
                seen.add(name)

            for col in ("fastq_1", "fastq_2"):
                f = (row[col] or "").strip()
                if not f:
                    errors.append(f"line {i}: empty {col}")
                elif not Path(f).is_file():
                    errors.append(f"line {i}: file not found -> {f}")

    if errors:
        print("\n".join(errors))
        sys.exit(f"Sample sheet INVALID ({len(errors)} problem(s)).")

    print(f"Sample sheet OK: {len(seen)} samples, all FASTQ files found.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: python validate_samplesheet.py samplesheet.csv")
    main(sys.argv[1])






---
```
## Bước 3 — QC (sequencing data QC)

```bash
mkdir -p qc trim
fastqc raw/*.sub.fq.gz -o qc/

for SRR in SRRxxxxxx1 ... ; do
  fastp \
    -i raw/${SRR}_1.fq.gz -I raw/${SRR}_2.fq.gz \
    -o trim/${SRR}_1.fq.gz   -O trim/${SRR}_2.fq.gz \
    -j qc/${SRR}.fastp.json  -h qc/${SRR}.fastp.html
done

multiqc qc/ -o qc/        # gộp toàn bộ QC thành 1 report
```

---

## Bước 4 — Alignment + sort + markdup (alignment to reference)

```bash
mkdir -p bam
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
```

---

## Bước 5 — Variant calling + filtering → analysis-ready genotypes

```bash
mkdir -p vcf
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
```

---

## Bước 6 — Population genetics (đánh đúng từng bullet downstream của JD)

```bash
mkdir -p pop
# Nạp vào plink2 (--allow-extra-chr vì tên scaffold không chuẩn 1..22)
plink2 --vcf vcf/snps.vcf.gz --allow-extra-chr \
       --set-all-var-ids @:#:\$r:\$a --make-pgen --out pop/bsf
```

**PCA**
```bash
plink2 --pfile pop/bsf --allow-extra-chr --pca 10 --out pop/bsf_pca
# -> pop/bsf_pca.eigenvec : vẽ PC1 vs PC2 bằng R/Python, tô màu theo population
```

**Heterozygosity + FIS (hệ số cận huyết F per mẫu)**
```bash
plink2 --pfile pop/bsf --allow-extra-chr --het --out pop/bsf_het
# cột F = inbreeding coefficient; O(HOM)/E(HOM) cho heterozygosity quan sát vs kỳ vọng
```

**Relatedness / kinship**
```bash
plink2 --pfile pop/bsf --allow-extra-chr --make-king-table --out pop/bsf_king
# -> ma trận kinship KING giữa các cặp mẫu
```

**ROH-based inbreeding**
```bash
bcftools roh -G30 --AF-dflt 0.4 vcf/snps.vcf.gz > pop/bsf.roh.txt
# tổng chiều dài ROH / chiều dài genome ~ F_ROH (chỉ số cận huyết dựa ROH)
# hoặc: plink2 ... --homozyg
```

---
