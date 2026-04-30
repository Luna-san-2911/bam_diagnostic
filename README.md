# bam_diagnostic

Hola! Hi!

This is a small repo created to automate a routine task for which I couldn't find an adequate pipeline. While there are many tools for alignment stats, summaries, and depth, I haven't found one that easily generates clean diagnostic plots with useful information (e.g., mapping percentage, depth, and coverage). That's the purpose of this repository! :)

To get started, download the three scripts from this repository and note their file paths.

Depending on your data, choose one of the following:
1. 1 bam file to analyze, try [this](https://github.com/Luna-san-2911/bam_diagnostic/blob/main/summary_statistics_bam_alignment.sbatch)
2. 1 directory with multiple bam files, try [this](https://github.com/Luna-san-2911/bam_diagnostic/blob/main/summary_statistics_bam_alignment_for_directory.sbatch)

Note that these are ready-to-submit .sbatch files designed for SLURM job schedulers. The .sbatch scripts require the provided R script to run, so make sure it is downloaded as well.

For both .sbatch files, you only need to adjust:
- The path to the R script.
- The path to your BAM file or your directory containing the BAM files.
- The path to your Conda source (line 29 [here](https://github.com/Luna-san-2911/bam_diagnostic/blob/main/summary_statistics_bam_alignment.sbatch))
- Your Conda environments. In the example scripts, I activate an environment (whatshap-env; line 30 [here](https://github.com/Luna-san-2911/bam_diagnostic/blob/main/summary_statistics_bam_alignment.sbatch)) with samtools and mosdepth installed, followed by an environment (random_forest, line 44 [here](https://github.com/Luna-san-2911/bam_diagnostic/blob/main/summary_statistics_bam_alignment.sbatch))containing R and its required packages. If you have a single environment with all these dependencies, simply activate that one.

Please note:
- This pipeline currently supports genome alignments with up to 24 chromosomes, as the color palette for the depth distribution plot is limited to 25 distinguishable colors.

Enjoy! If you have any questions, please feel free to contact me. :)
