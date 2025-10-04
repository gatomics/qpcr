import subprocess
import sys
import time
import os

# Check for minimum number of input files (at least one pair)
if len(sys.argv) < 3 or len(sys.argv) % 2 != 1:
    print("Usage: python3 bioinformatics_pipeline.py <file1_R1.fastq.gz> <file1_R2.fastq.gz> [<file2_R1.fastq.gz> <file2_R2.fastq.gz> ...]")
    sys.exit(1)

# Extract paired FASTQ files from command line arguments
paired_fastq_files = sys.argv[1:]

# Define paths to tools and resources
BWA_MEM2_PATH = "/path/to/bwa-mem2"
PICARD_PATH = "/path/to/picard.jar"
REFERENCE_GENOME = "/path/to/reference_genome/hg38.fasta"
MANTA_CONFIG_PATH = "/path/to/manta/bin/configManta.py"
EXPANSIONHUNTER_PATH = "/path/to/ExpansionHunter"
STRANGER_PATH = "/path/to/stranger/resources/variant_catalog.json"
SNPEFF_PATH = "/path/to/snpEff/snpEff.jar"
DEEPVARIANT_IMAGE = "google/deepvariant:1.2.0"
BASE_RESULTS_DIR = "./Results"  # Base directory for results

def run_command(command):
    """Run a shell command and print its runtime."""
    start_time = time.time()
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    end_time = time.time()

    elapsed_time = end_time - start_time
    formatted_time = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
    print(f"Command '{command}' took {formatted_time} to complete.")

    if process.returncode != 0:
        print(f"Command '{command}' failed with return code: {process.returncode}")
        print(stderr.decode())
        sys.exit(1)

# Loop over each pair of input files
for i in range(0, len(paired_fastq_files), 2):
    input_file_1 = paired_fastq_files[i]
    input_file_2 = paired_fastq_files[i+1]
    
    # Generate a unique directory name for each pair
    pair_id = os.path.splitext(os.path.basename(input_file_1))[0]  # Customize as needed for uniqueness
    results_dir = os.path.join(BASE_RESULTS_DIR, pair_id)
    
    # Create directories for each pair
    os.makedirs(results_dir, exist_ok=True)
    snv_dir = os.path.join(results_dir, "SNV")
    sv_dir = os.path.join(results_dir, "SV")
    manta_run_dir = os.path.join(sv_dir, "manta")
    os.makedirs(snv_dir, exist_ok=True)
    os.makedirs(sv_dir, exist_ok=True)
    os.makedirs(manta_run_dir, exist_ok=True)

    # Align reads using bwa-mem2 with the input files
    run_command(f"{BWA_MEM2_PATH} mem {REFERENCE_GENOME} {input_file_1} {input_file_2} > {results_dir}/{pair_id}_aligned.sam")
    
    # Convert SAM to BAM, sort BAM, and mark duplicates for each pair
    run_command(f"samtools view -bS {results_dir}/{pair_id}_aligned.sam > {results_dir}/{pair_id}_aligned.bam")
    run_command(f"samtools sort {results_dir}/{pair_id}_aligned.bam -o {results_dir}/{pair_id}_aligned.sorted.bam")
    run_command(f"java -jar {PICARD_PATH} MarkDuplicates I={results_dir}/{pair_id}_aligned.sorted.bam O={results_dir}/{pair_id}_aligned.sorted.marked.bam M={results_dir}/{pair_id}_marked_dup_metrics.txt")

    # Index BAM file and validate it for each pair
    run_command(f"samtools index {results_dir}/{pair_id}_aligned.sorted.marked.bam")
    run_command(f"java -jar {PICARD_PATH} ValidateSamFile I={results_dir}/{pair_id}_aligned.sorted.marked.bam MODE=SUMMARY")

    # Add or replace read groups for each pair
    run_command(f"java -jar {PICARD_PATH} AddOrReplaceReadGroups I={results_dir}/{pair_id}_aligned.sorted.marked.bam O={results_dir}/{pair_id}_aligned.sorted.marked_with_rg.bam RGID={pair_id} RGLB={pair_id}_LB RGPL=illumina RGPU={pair_id}_PU RGSM={pair_id}")

    # Index the BAM file again for each pair
    run_command(f"samtools index {results_dir}/{pair_id}_aligned.sorted.marked_with_rg.bam")

    deepvariant_cmd = " --model_type=WES --ref=/ref/{REFERENCE_GENOME.split('/')[-1]} --reads=/input/{pair_id}_aligned.sorted.marked_with_rg.bam --output_vcf=/output/{pair_id}_output_variants.vcf --output_gvcf=/output/{pair_id}_output_variants.g.vcf.gz --num_shards=8"
    run_command(f"sudo docker run -v {os.path.abspath(results_dir)}:/input -v {os.path.abspath(REFERENCE_GENOME.rsplit('/', 1)[0])}:/ref -v {os.path.abspath(snv_dir)}:/output {DEEPVARIANT_IMAGE} {deepvariant_cmd.format(pair_id=pair_id)}")

    # For configuring and running Manta for structural variant calling
    manta_run_dir = os.path.join(sv_dir, "manta")
    os.makedirs(manta_run_dir, exist_ok=True)
    run_command(f'"{MANTA_CONFIG_PATH}" --bam "{results_dir}/{pair_id}_aligned.sorted.marked_with_rg.bam" --referenceFasta "{REFERENCE_GENOME}" --runDir "{manta_run_dir}"')
    run_command(f"{manta_run_dir}/runWorkflow.py -m local")

    # For calling structural variants with Delly
    for sv_type in ["DEL", "DUP"]:
        output_bcf = f"{sv_dir}/delly_{sv_type.lower()}.bcf"
        run_command(f"delly call -t {sv_type} -g {REFERENCE_GENOME} -o {output_bcf} {results_dir}/{pair_id}_aligned.sorted.marked_with_rg.bam")
        run_command(f"bcftools view {output_bcf} > {output_bcf.replace('.bcf', '.vcf')}")

    # For merging SVs with SVDB
    run_command(f"svdb --merge --notag --vcf {sv_dir}/manta/results/variants/*.vcf.gz --vcf {sv_dir}/delly_*.vcf > {sv_dir}/merged_sv.vcf")

    # For running ExpansionHunter for repeat expansions
    repeat_expansions_dir = os.path.join(results_dir, "repeatExpansions")
    os.makedirs(repeat_expansions_dir, exist_ok=True)
    run_command(f"{EXPANSIONHUNTER_PATH} --reads {results_dir}/{pair_id}_aligned.sorted.marked_with_rg.bam --reference {REFERENCE_GENOME} --variant-catalog {EXPANSIONHUNTER_PATH.rsplit('/', 1)[0]}/variant_catalog/hg38/variant_catalog.json --output-prefix {repeat_expansions_dir}/{pair_id}_eh_output")

    # For annotating with stranger
    run_command(f"stranger {repeat_expansions_dir}/{pair_id}_eh_output.vcf -f {STRANGER_PATH} > {repeat_expansions_dir}/{pair_id}_annotated_output.vcf")

    # For annotating SNVs and SVs with SnpEff
    run_command(f"java -Xmx8g -jar {SNPEFF_PATH} -v hg38 {snv_dir}/{pair_id}_output_variants.g.vcf.gz > {snv_dir}/{pair_id}_snv.ann.vcf")
    run_command(f"java -Xmx8g -jar {SNPEFF_PATH} -v hg38 {sv_dir}/merged_sv.vcf > {sv_dir}/{pair_id}_sv.ann.vcf")

print("Pipeline modifications for processing each pair of FASTQ files independently completed successfully.")


