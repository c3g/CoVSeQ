#!/usr/bin/env bash

readset=$1

if [[ -z "$readset" ]]; then
   break
fi

# module load mugqic/python/3.7.3
# module load mugqic/sambamba/0.7.0
# module load mugqic/samtools/1.10

cat /dev/null > metrics/metrics.csv
cat /dev/null > metrics/host_contamination_metrics.tsv
cat /dev/null > metrics/host_removed_metrics.tsv
cat /dev/null > metrics/kraken2_metrics.tsv

echo "sample,cons.per.N,cons.len,cons.perc.GC,cons.perc.genome_frac,cons.N_per_100_kbp,fq.trim.pass,bam.perc.align,bam.filter.pass,bam.primertrim.pass,bam.mean.cov,bam.med.cov,bam.max-min/mean.cov,bam.perc.20x,bam.perc.50x,bam.perc.100x,bam.perc.250x,bam.perc.500x,bam.perc.1000x,bam.perc.2000x,bam.mean.insertsize,bam.med.insertsize,bam.sd.insertsize,bam.min.insertsize,bam.max.insertsize" > metrics/metrics.csv

echo -e "Sample\tTotal_aligned\tHuman_only\tHuman_only_perc\tSARS_only\tSARS_only_perc\tUnmapped_only\tUnmapped_only_perc" > metrics/host_contamination_metrics.tsv
echo -e "Sample\tTotal_aligned\tHuman_only\tHuman_only_perc\tSARS_only\tSARS_only_perc\tUnmapped_only\tUnmapped_only_perc" > metrics/host_removed_metrics.tsv

echo -e "Sample\tHomo_sapiens_clade\tHomo_sapiens_clade_perc" > metrics/kraken2_metrics.tsv

for sample in `awk 'NR>1 {print $1}' $readset`
do
    genome_size=29903
    consensus_fa=alignment/$sample/$sample.sorted.filtered.primerTrim.consensus.fa

    # Parsing quast results
    quast_tsv=`echo "metrics/dna/$sample/quast_metrics/report.tsv"`
    quast_html=`echo "metrics/dna/$sample/quast_metrics/report.html"`

    if [ -f "$quast_tsv" -a -f "$quast_html" ]; then
        cons_len=`grep -oP "Total length \(>= 0 bp\)\t\K.*?(?=$)" $quast_tsv`
        N_count=`grep -oP "# N's\",\"quality\":\"Less is better\",\"values\":\[\K.*?(?=])" $quast_html`
        if [ "$cons_len" != "0" ]; then
            cons_perc_N=`echo "scale=2; 100*$N_count/$cons_len" | bc -l`
        else
            cons_perc_N="NULL"
        fi
        cons_GC=`grep -oP "^GC \(%\)\t\K.*?(?=$)" $quast_tsv`
        cons_genome_frac=`grep -oP "^Genome fraction \(%\)\t\K.*?(?=$)" $quast_tsv`
        cons_N_perkbp=`grep -oP "^# N's per 100 kbp\t\K.*?(?=$)" $quast_tsv`
    else
        cons_len="NULL"
        N_count="NULL"
        cons_perc_N="NULL"
        cons_GC="NULL"
        cons_genome_frac="NULL"
        cons_N_perkbp="NULL"
    fi

    if [ -z $cons_GC ]; then
        cons_GC="NULL"
    fi
    if [ -z $cons_genome_frac ]; then
        cons_genome_frac="NULL"
    fi
    if [ -z $cons_N_perkbp ]; then
        cons_N_perkbp="NULL"
    fi

    # Parsing cutadapt results
    readset_name=`grep "$sample" $readset | awk '{print $2}'`
    cutadapt_file=`ls -t job_output/cutadapt/cutadapt.${readset_name}_*[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]T[0-9][0-9].[0-9][0-9].[0-9][0-9].o | head -n 1`
    fq_surviving_trim=`grep -oP 'Pairs written \(passing filters\):.*\(\K.*?(?=%)' $cutadapt_file`

    # Parsing sambamba flagstat and picard metrics
    raw_flagstat_file=`echo "metrics/dna/$sample/flagstat/$sample.sorted.flagstat"`
    filtered_flagstat_file=`echo "metrics/dna/$sample/flagstat/$sample.sorted.filtered.flagstat"`
    primertrim_flagstat_file=`echo "metrics/dna/$sample/flagstat/$sample.sorted.filtered.primerTrim.flagstat"`
    total_raw=`head -n 1 $raw_flagstat_file | grep -oP '^\K.*?(?= )'`
    alignment_metrics_file=`echo "alignment/$sample/$sample.sorted.filtered.all.metrics.alignment_summary_metrics"`
    total_filter=`head -n 1 $filtered_flagstat_file | grep -oP '^\K.*?(?= )'`
    if [ "$total_raw" != "0" -a "$total_filter" != "0" -a -f "$alignment_metrics_file" ]; then
        bam_surviving_filter=`echo "scale=2; 100*$total_filter/$total_raw" | bc -l`
        total_primertrim=`head -n 1 $primertrim_flagstat_file | grep -oP '^\K.*?(?= )'`
        bam_surviving_primertrim=`echo "scale=2; 100*$total_primertrim/$total_filter" | bc -l`
    else
        bam_surviving_filter="NULL"
        bam_surviving_primertrim="NULL"
    fi
    bam_aln=`grep "mapped (" $raw_flagstat_file | grep -oP '^.*\(\K.*?(?=%)'`

    # Parsing BedGraph file
    bam_bedgraph_file=`echo "alignment/$sample/$sample.sorted.filtered.BedGraph"`
    bam_meancov=`awk '{for (i = 1; i <= $3-$2; ++i) sum+=$4} END {print sum}' $bam_bedgraph_file`
    bam_meancov=`echo "scale=2; $bam_meancov/$genome_size" | bc -l`
    bam_mediancov=`awk '{for (i = 1; i <= $3-$2; ++i)  print $4}' $bam_bedgraph_file | sort -g | awk '{count[NR] = $1} END {if (NR % 2) {print count[(NR + 1) / 2];} else {print (count[(NR / 2)] + count[(NR / 2) + 1]) / 2.0;}}'`
    bam_mincov=`awk '{print $4}' $bam_bedgraph_file | sort -g | head -n 1 | awk '{printf "%i", $0}'`
    bam_maxcov=`awk '{print $4}' $bam_bedgraph_file | sort -g | tail -n 1 | awk '{printf "%i", $0}'`
    if [ "$bam_meancov" == "0" ]; then
        bam_maxmincovmean="NULL"
    else
        bam_maxmincovmean=`echo "scale=2; ($bam_maxcov-$bam_mincov)/$bam_meancov" | bc -l`
    fi
    bam_cov20X=`awk '{if ($4 > 20) {count = count + $3-$2}} END {if (count) {print count} else {print 0}}' $bam_bedgraph_file`
    bam_cov20X=`echo "scale=2; 100*$bam_cov20X/$genome_size" | bc -l`
    bam_cov50X=`awk '{if ($4 > 50) {count = count + $3-$2}} END {if (count) {print count} else {print 0}}' $bam_bedgraph_file`
    bam_cov50X=`echo "scale=2; 100*$bam_cov50X/$genome_size" | bc -l`
    bam_cov100X=`awk '{if ($4 > 100) {count = count + $3-$2}} END {if (count) {print count} else {print 0}}' $bam_bedgraph_file`
    bam_cov100X=`echo "scale=2; 100*$bam_cov100X/$genome_size" | bc -l`
    bam_cov250X=`awk '{if ($4 > 250) {count = count + $3-$2}} END {if (count) {print count} else {print 0}}' $bam_bedgraph_file`
    bam_cov250X=`echo "scale=2; 100*$bam_cov250X/$genome_size" | bc -l`
    bam_cov500X=`awk '{if ($4 > 500) {count = count + $3-$2}} END {if (count) {print count} else {print 0}}' $bam_bedgraph_file`
    bam_cov500X=`echo "scale=2; 100*$bam_cov500X/$genome_size" | bc -l`
    bam_cov1000X=`awk '{if ($4 > 1000) {count = count + $3-$2}} END {if (count) {print count} else {print 0}}' $bam_bedgraph_file`
    bam_cov1000X=`echo "scale=2; 100*$bam_cov1000X/$genome_size" | bc -l`
    bam_cov2000X=`awk '{if ($4 > 2000) {count = count + $3-$2}} END {if (count) {print count} else {print 0}}' $bam_bedgraph_file`
    bam_cov2000X=`echo "scale=2; 100*$bam_cov2000X/$genome_size" | bc -l`

    # Parsing Picard insert size metrics
    bam_insertsize_file=`echo "metrics/dna/$sample/picard_metrics/$sample.sorted.all.metrics.insert_size_metrics"`
    if [ -f "$bam_insertsize_file" ]; then
        bam_meaninsertsize=`awk 'NR==8 {print $6}' $bam_insertsize_file`
        bam_meaninsertsize=`echo "scale=2; ${bam_meaninsertsize//,}/1" | bc -l`
        bam_medianinsertsize=`awk 'NR==8 {print $1}' $bam_insertsize_file`
        bam_mininsertsize=`awk 'NR==8 {print $4}' $bam_insertsize_file`
        bam_maxinsertsize=`awk 'NR==8 {print $5}' $bam_insertsize_file`
        bam_sdinsertsize=`awk 'NR==8 {print $7}' $bam_insertsize_file`
        if [ "$bam_sdinsertsize" != "?" ]; then
            bam_sdinsertsize=`echo "scale=2; ${bam_sdinsertsize//,}/1" | bc -l`
        else
            bam_sdinsertsize="NULL"
        fi
    else 
        bam_meaninsertsize="NULL"
        bam_medianinsertsize="NULL"
        bam_mininsertsize="NULL"
        bam_maxinsertsize="NULL"
        bam_sdinsertsize="NULL"
    fi

    # Computing host contamination and host cleaning metrics
    samtools idxstats -@ 10 host_removal/${sample}/${sample}*.hybrid.sorted.bam | awk -v sample=$sample '{array[$1]=$3; pwet[$1]=$4; next} END {tot=0; unmapped=0; for (chr in array) {tot+=array[chr]; unmapped+=pwet[chr]}; tot+=unmapped; cov=array["MN908947.3"]; hum=tot-cov-unmapped; if (tot <= 0) {printf sample"\t"tot"\t"hum"\tNULL\t"cov"\tNULL\t"unmapped"\tNULL\n"} else {printf sample"\t"tot"\t"hum"\t%.2f\t"cov"\t%.2f\t"unmapped"\t%.2f\n", 100*hum/tot, 100*cov/tot, 100*unmapped/tot}}' >> metrics/host_contamination_metrics.tsv
    samtools idxstats -@ 10 host_removal/${sample}/${sample}*.host_removed.sorted.bam | awk -v sample=$sample '{array[$1]=$3; pwet[$1]=$4; next} END {tot=0; unmapped=0; for (chr in array) {tot+=array[chr]; unmapped+=pwet[chr]}; tot+=unmapped; cov=array["MN908947.3"]; hum=tot-cov-unmapped; if (tot <= 0) {printf sample"\t"tot"\t"hum"\tNULL\t"cov"\tNULL\t"unmapped"\tNULL\n"} else {printf sample"\t"tot"\t"hum"\t%.2f\t"cov"\t%.2f\t"unmapped"\t%.2f\n", 100*hum/tot, 100*cov/tot, 100*unmapped/tot}}' >> metrics/host_removed_metrics.tsv

    kraken_file=`ls metrics/dna/${sample}/kraken_metrics/${sample}*.kraken2_report`
    if [ -s "$kraken_file" ]; then
        grep "Homo sapiens" $kraken_file | awk -v sample=$sample '{print sample"\t"$2"\t"$1}' >> metrics/kraken2_metrics.tsv
    else
        echo -e "$sample\tNULL\tNULL" >> metrics/kraken2_metrics.tsv
    fi

    echo "$sample,$cons_perc_N,$cons_len,$cons_GC,$cons_genome_frac,$cons_N_perkbp,$fq_surviving_trim,$bam_aln,$bam_surviving_filter,$bam_surviving_primertrim,$bam_meancov,$bam_mediancov,$bam_maxmincovmean,$bam_cov20X,$bam_cov50X,$bam_cov100X,$bam_cov250X,$bam_cov500X,$bam_cov1000X,$bam_cov2000X,$bam_meaninsertsize,$bam_medianinsertsize,$bam_sdinsertsize,$bam_mininsertsize,$bam_maxinsertsize" >> metrics/metrics.csv
    # break
done
