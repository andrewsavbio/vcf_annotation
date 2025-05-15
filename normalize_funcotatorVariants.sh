#!/bin/bash

#SBATCH --job-name=somaticFuncotator      ## Name of the job.
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --ntasks=1           ## (-n) number of tasks to launch
#SBATCH --cpus-per-task=60    ## number of threads the job needs
#SBATCH --mem=120gb              ## Amount of RAM to use
#SBATCH --error=/home/workstation_v1/Documents/bashScripts/Slurm_Error_Reports/slurm-%J.err             ## error log file
#SBATCH --out=/home/workstation_v1/Documents/bashScripts/Slurm_Error_Reports/slurm-%J.out               ## output log file

#Input:Mapped BAM files with metadata and recalibrated base quality scores
#Output: .g.vcf files with indels and SNVs

#https://gatk.broadinstitute.org/hc/en-us/articles/360035531132
#########################################################################################################################################

#expDir="/home/NAS1/Processed/Sequencing/WGS/20240621_patients/output_parabricks/" #path to experimental directory
expDir="/home/processing_1/output_parabricks_pelNGs_2687_2636/"
#expDir="/home/processing_1/output_parabricks_ibaNG_313/"

refFasta="/home/processing_1/parabricks_sample/Ref/hg38_v0_Homo_sapiens_assembly38.fasta" #Reference genome in FASTA format
somaticRefDir="/home/processing_1/parabricks_sample/Ref/Funcotator/Somatic/"

N=16 #number of tasks to run in parallel; edit the cpus-per-task and mem in the header as needed (~2gb/CPU)

#########################################################################################################################################

cd ${expDir}

startTime=`date +%s`

dRefFasta="/genomeRefDir/$(basename ${refFasta})"


for vcfFile in $(find ~+ -maxdepth 1 -type f  -name "*_mutect2.vcf")
    do
    (
        if [ ! -f "${vcfFile%%.*}_funcotator.vcf" ]; then
            docker run \
                        --workdir /workdir \
                        --rm \
                        --volume $(pwd):/workdir \
                        --volume ${somaticRefDir}:/funcotatorRefDir \
                        --volume $(dirname ${refFasta}):/genomeRefDir \
                        broadinstitute/gatk \
                        gatk Funcotator \
                        --data-sources-path /funcotatorRefDir/ \
                        --ref-version hg38 \
                        -R ${dRefFasta} \
                        -V /workdir/$(basename ${vcfFile}) \
                        -O /workdir/$(basename ${vcfFile%%.*})_funcotator.vcf.gz \
                        --output-file-format VCF
            #annotate variants. documentation says --variant, but should be --input or -I? -V works for some reason despite error message
        fi

        if [ ! -f "${vcfFile%%.*}_SNPs_filtered_annotated.vcf.gz" ]; then
            docker run \
                        --workdir /workdir \
                        --rm \
                        --volume $(pwd):/workdir \
                        --volume ${somaticRefDir}:/funcotatorRefDir \
                        --volume $(dirname ${refFasta}):/genomeRefDir \
                        broadinstitute/gatk \
                        gatk SelectVariants \
                        -R ${dRefFasta} \
                        -V /workdir/$(basename ${vcfFile%%.*})_funcotator.vcf.gz \
                        --select-type-to-include SNP \
                        -O /workdir/$(basename ${vcfFile%%.*})_SNPs_filtered_annotated.vcf.gz \
                        --exclude-filtered --exclude-non-variants --remove-unused-alternates
            #select for SNP variants that passed filters
        fi

        if [ ! -f "${vcfFile%%.*}_MNPs_filtered_annotated.vcf.gz" ]; then
            docker run \
                        --workdir /workdir \
                        --rm \
                        --volume $(pwd):/workdir \
                        --volume ${somaticRefDir}:/funcotatorRefDir \
                        --volume $(dirname ${refFasta}):/genomeRefDir \
                        broadinstitute/gatk \
                        gatk SelectVariants \
                        -R ${dRefFasta} \
                        -V /workdir/$(basename ${vcfFile%%.*})_funcotator.vcf.gz \
                        --select-type-to-include MNP \
                        -O /workdir/$(basename ${vcfFile%%.*})_MNPs_filtered_annotated.vcf.gz \
                        --exclude-filtered --exclude-non-variants --remove-unused-alternates
            #select for MNP variants that passed filters
        fi

        if [ ! -f "${vcfFile%%.*}_INDELs_filtered_annotated.vcf.gz" ]; then
            docker run \
                        --workdir /workdir \
                        --rm \
                        --volume $(pwd):/workdir \
                        --volume ${somaticRefDir}:/funcotatorRefDir \
                        --volume $(dirname ${refFasta}):/genomeRefDir \
                        broadinstitute/gatk \
                        gatk SelectVariants \
                        -R ${dRefFasta} \
                        -V /workdir/$(basename ${vcfFile%%.*})_funcotator.vcf.gz \
                        --select-type-to-include INDEL \
                        -O /workdir/$(basename ${vcfFile%%.*})_INDELs_filtered_annotated.vcf.gz \
                        --exclude-filtered --exclude-non-variants --remove-unused-alternates
            #select for INDEL variants that passed filters
        fi

        if [ ! -f "${vcfFile%%.*}_SNPs_filtered_annotated.ann.out" ]; then
            docker run \
                    --workdir /workdir \
                    --rm \
                    --volume $(pwd):/workdir \
                    broadinstitute/gatk \
                    gatk VariantsToTable \
                    -V /workdir/$(basename ${vcfFile%%.*})_SNPs_filtered_annotated.vcf.gz \
                    -F CHROM -F POS -F TYPE -F FUNCOTATION \
                    -O /workdir/$(basename ${vcfFile%%.*})_SNPs_filtered_annotated.ann.out
            #create table of variant outputs
        fi
        
        if [ ! -f "${vcfFile%%.*}_MNPs_filtered_annotated.ann.out" ]; then
            docker run \
                    --workdir /workdir \
                    --rm \
                    --volume $(pwd):/workdir \
                    broadinstitute/gatk \
                    gatk VariantsToTable \
                    -V /workdir/$(basename ${vcfFile%%.*})_MNPs_filtered_annotated.vcf.gz \
                    -F CHROM -F POS -F TYPE -F FUNCOTATION \
                    -O /workdir/$(basename ${vcfFile%%.*})_MNPs_filtered_annotated.ann.out
            #create table of variant outputs
        fi
        
        if [ ! -f "${vcfFile%%.*}_INDELs_filtered_annotated.ann.out" ]; then
            docker run \
                    --workdir /workdir \
                    --rm \
                    --volume $(pwd):/workdir \
                    broadinstitute/gatk \
                    gatk VariantsToTable \
                    -V /workdir/$(basename ${vcfFile%%.*})_INDELs_filtered_annotated.vcf.gz \
                    -F CHROM -F POS -F TYPE -F FUNCOTATION \
                    -O /workdir/$(basename ${vcfFile%%.*})_INDELs_filtered_annotated.ann.out
            #create table of variant outputs
        fi

    ) &

    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
            # now there are $N jobs already running, so wait here for any job
            # to be finished so there is a place to start next one.
            wait -n
        fi

done

wait #wait for last tasks to finish

#check runtime
endTime=`date +%s`
runTime=$((endTime-startTime))
runHours=$(echo "$runTime/3600" | bc -l)

echo "Somatic Funcotator: "
echo -e "Start: ${startTime}\nEnd: ${endTime}\nRuntime (sec): ${runTime}\nRuntime (h): ${runHours}"
