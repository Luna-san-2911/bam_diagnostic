#!/bin/bash -l
#SBATCH -J summary_statistics_bam_alignment
#SBATCH --output=summary_statistics_bam_alignment-%j.output
#SBATCH --error=summary_statistics_bam_alignment-%j.error
#SBATCH -t 120:00:00
#SBATCH --partition=bmh
#SBATCH --mail-type=ALL
#SBATCH --mail-user=plunarodriguez@ucdavis.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

#Inside de whatshap-env
set -euo pipefail

echo "Running on host: $(hostname)"
echo "Starting at: $(date)"

BAM_DIR="/group/gmonroegrp2/chaehee/Pablo/pangenie/05_map_sample_to_infered_haplotype/01_align_reference_genomes_clone10/01_bam"
OUT_DIR="/group/gmonroegrp2/chaehee/Pablo/pangenie/05_map_sample_to_infered_haplotype/01_align_reference_genomes_clone10/02_stats"
OUT_FILES=$OUT_DIR/01_files
OUT_PLOTS=$OUT_DIR/02_plots

if [[ ! -d "$BAM_DIR" ]]; then
        echo "ERROR: BAM directory does not exist: $BAM_DIR"
        exit 1
fi

mkdir -p "$OUT_FILES" "$OUT_PLOTS"

shopt -s nullglob
bam_files=("$BAM_DIR"/*.bam)

if [[ ${#bam_files[@]} -eq 0 ]]; then
        echo "ERROR: No .bam files found in: $BAM_DIR"
        exit 1
fi

echo "Found ${#bam_files[@]} BAM files in: $BAM_DIR"

source /home/plunarod/miniconda3/bin/activate
conda activate whatshap-env

for BAM in "${bam_files[@]}"; do
        SAMPLE_NAME=$(basename "$BAM" .bam)
        echo "--- Processing sample: $SAMPLE_NAME ---"

        echo "--- Generating Mapping Statistics (Flagstat) ---"
        samtools flagstat -@ 16 "$BAM" > "$OUT_FILES/${SAMPLE_NAME}_flagstat_metrics.txt"

        echo "--- Extracting MAPQ Score Frequencies ---"
        samtools stats -@ 16 "$BAM" | grep '^MAPQ' | cut -f 2,3 > "$OUT_FILES/${SAMPLE_NAME}_mapq_frequencies.txt"

        echo "--- Generating 100kb Windowed Depth using mosdepth ---"
        mosdepth -x -n -t 16 -b 100000 "$OUT_FILES/${SAMPLE_NAME}_depth_100kb" "$BAM"
done

echo "Summary statistics bam alignment finished at: $(date)"

conda deactivate
conda activate random_forest

for BAM in "${bam_files[@]}"; do
    SAMPLE_NAME=$(basename "$BAM" .bam)
    echo "--- Generating plots for sample: $SAMPLE_NAME ---"
    Rscript /home/plunarod/pistachio/code/04_pangenie/summary_statistics_bam_alignment.R \
            "$OUT_FILES/${SAMPLE_NAME}_flagstat_metrics.txt" \
            "$OUT_FILES/${SAMPLE_NAME}_mapq_frequencies.txt" \
            "$OUT_FILES/${SAMPLE_NAME}_depth_100kb.mosdepth.global.dist.txt" \
            "$OUT_FILES/${SAMPLE_NAME}_depth_100kb.regions.bed.gz" \
            "$OUT_PLOTS/${SAMPLE_NAME}"
done
