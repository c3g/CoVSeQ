readset=$1

module load mugqic/python/3.7.3
# module load mugqic/sambamba/0.7.0

cat /dev/null > metrics/metrics.csv

echo "sample,cons.per.N,cons.len,cons.perc.GC,fq.trim.pass,bam.perc.align,bam.filter.pass,bam.primertrim.pass,bam.mean.cov,bam.med.cov,bam.max-min/mean.cov,bam.perc.50x,bam.perc.100x,bam.perc.250x,bam.perc.500x,bam.perc.1000x,bam.perc.2000x,bam.mean.insertsize,bam.med.insertsize,bam.sd.insertsize,bam.min.insertsize,bam.max.insertsize" > metrics/metrics.csv

for sample in `awk 'NR>1 {print $1}' $readset`
do
    genome_size=29903
    consensus_fa=alignment/$sample/$sample.sorted.filtered.primerTrim.consensus.fa

    N_count=`~/scripts/fasta_stats.py $consensus_fa | grep -oP 'N : \K.*?(?= )'`
    if [ -z "$N_count" ]; then
        cons_perc_N="NULL"
    else
        cons_perc_N=`echo "scale=2; 100*$N_count/$genome_size" | bc -l`
    fi

    cons_len=`~/scripts/fasta_stats.py $consensus_fa | grep -oP 'contains \K.*?(?= )'`

    G_count=`~/scripts/fasta_stats.py $consensus_fa | grep -oP 'G : \K.*?(?= )'`
    C_count=`~/scripts/fasta_stats.py $consensus_fa | grep -oP 'C : \K.*?(?= )'`
    if [ -z "$G_count" -o -z "$C_count" ]; then
        cons_GC="NULL"
    else
        cons_GC=`echo "scale=2; 100*($G_count+$G_count)/$genome_size" | bc -l`
    fi

    cutadapt_file=`ls -t job_output/cutadapt/cutadapt.${sample}_[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]T*.*.*.o | head -n 1`
    fq_surviving_trim=`grep -oP 'Pairs written \(passing filters\):.*\(\K.*?(?=%)' $cutadapt_file`


    # total_raw=`sambamba flagstat alignment/$sample/$sample.sorted.bam | head -n 1 | grep -oP '^\K.*?(?= )'`
    raw_flagstat_file=`echo "metrics/dna/$sample/flagstat/$sample.sorted.flagstat"`
    filtered_flagstat_file=`echo "metrics/dna/$sample/flagstat/$sample.sorted.filtered.flagstat"`
    primertrim_flagstat_file=`echo "metrics/dna/$sample/flagstat/$sample.sorted.filtered.primerTrim.flagstat"`
    total_raw=`head -n 1 $raw_flagstat_file | grep -oP '^\K.*?(?= )'`
    alignment_metrics_file=`echo "alignment/$sample/$sample.sorted.filtered.all.metrics.alignment_summary_metrics"`
    if [ "$total_raw" != "0" -a -f "$alignment_metrics_file" ]; then
        # total_primertrim=`awk '{if ($1=="PAIR") {print $2}}' $alignment_metrics_file`
        total_filter=`head -n 1 $filtered_flagstat_file | grep -oP '^\K.*?(?= )'`
        bam_surviving_filter=`echo "scale=2; 100*$total_filter/$total_raw" | bc -l`
        total_primertrim=`head -n 1 $primertrim_flagstat_file | grep -oP '^\K.*?(?= )'`
        bam_surviving_primertrim=`echo "scale=2; 100*$total_primertrim/$total_filter" | bc -l`
        # bam_aln=`awk '{if ($1=="PAIR") {print $7}}' $alignment_metrics_file`
        # bam_aln=`echo "scale=2; (100*$bam_aln)/1" | bc -l`
    else 
        bam_surviving_primertrim="NULL"
        # bam_aln="NULL"
    fi
    bam_aln=`grep "mapped (" $raw_flagstat_file | grep -oP '^.*\(\K.*?(?=%)'`

    # total_primertrim=`awk '{if ($1=="PAIR") {print $2}}' $alignment_metrics_file`
    # bam_surviving_primertrim=`echo "scale=2; 100*$total_primertrim/$total_raw" | bc -l`

    # bam_aln=`awk '{if ($1=="PAIR") {print $7}}' $alignment_metrics_file`
    # bam_aln=`echo "scale=2; (100*$bam_aln)/1" | bc -l`

    # bam_ontarget=

    # bam_meancov=`grep -oP 'mean coverageData = \K.*?(?=X)' metrics/dna/$sample/qualimap/genome_results.txt`
    # bam_meancov=`echo "scale=2; ${bam_meancov//,}/1" | bc -l`

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

    # bam_cov50X=`awk 'NR==2 {print $17}' alignment/$sample/$sample.sorted.primerTrim.coverage.tsv`
    bam_cov50X=`awk '{if ($4 > 50) {count = count + $3-$2}} END {if (count) {print count} else {print 0}}' $bam_bedgraph_file`
    bam_cov50X=`echo "scale=2; 100*$bam_cov50X/$genome_size" | bc -l`
    # bam_cov100X=`awk 'NR==2 {print $19}' alignment/$sample/$sample.sorted.primerTrim.coverage.tsv`
    bam_cov100X=`awk '{if ($4 > 100) {count = count + $3-$2}} END {if (count) {print count} else {print 0}}' $bam_bedgraph_file`
    bam_cov100X=`echo "scale=2; 100*$bam_cov100X/$genome_size" | bc -l`
    # bam_cov250X=
    bam_cov250X=`awk '{if ($4 > 250) {count = count + $3-$2}} END {if (count) {print count} else {print 0}}' $bam_bedgraph_file`
    bam_cov250X=`echo "scale=2; 100*$bam_cov250X/$genome_size" | bc -l`
    # bam_cov500X=`awk 'NR==2 {print $21}' alignment/$sample/$sample.sorted.primerTrim.coverage.tsv`
    bam_cov500X=`awk '{if ($4 > 500) {count = count + $3-$2}} END {if (count) {print count} else {print 0}}' $bam_bedgraph_file`
    bam_cov500X=`echo "scale=2; 100*$bam_cov500X/$genome_size" | bc -l`
    # bam_cov1000X=`awk 'NR==2 {print $22}' alignment/$sample/$sample.sorted.primerTrim.coverage.tsv`
    bam_cov1000X=`awk '{if ($4 > 1000) {count = count + $3-$2}} END {if (count) {print count} else {print 0}}' $bam_bedgraph_file`
    bam_cov1000X=`echo "scale=2; 100*$bam_cov1000X/$genome_size" | bc -l`
    # bam_cov2000X=
    bam_cov2000X=`awk '{if ($4 > 2000) {count = count + $3-$2}} END {if (count) {print count} else {print 0}}' $bam_bedgraph_file`
    bam_cov2000X=`echo "scale=2; 100*$bam_cov2000X/$genome_size" | bc -l`

    bam_insertsize_file=`echo "alignment/$sample/$sample.sorted.all.metrics.insert_size_metrics"`
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
    # bam_meaninsertsize=`awk 'NR==8 {print $6}' $bam_insertsize_file`
    # bam_meaninsertsize=`echo "scale=2; ${bam_meaninsertsize//,}/1" | bc -l`
    # bam_medianinsertsize=`awk 'NR==8 {print $1}' $bam_insertsize_file`
    # bam_mininsertsize=`awk 'NR==8 {print $4}' $bam_insertsize_file`
    # bam_maxinsertsize=`awk 'NR==8 {print $5}' $bam_insertsize_file`
    # bam_sdinsertsize=`awk 'NR==8 {print $7}' $bam_insertsize_file`
    # bam_sdinsertsize=`echo "scale=2; ${bam_meaninsertsize//,}/1" | bc -l`

    # echo $sample, $bam_cov50X, $bam_cov100X, $bam_cov250X, $bam_cov500X, $bam_cov1000X, $bam_cov2000X
    # echo $sample $bam_meancov
    # echo "$sample,$cons_perc_N,$cons_len,$cons_GC,$fq_surviving_trim,$bam_surviving_primertrim,$bam_aln,$bam_meancov,$bam_mediancov,$bam_maxmincovmean,$bam_cov50X,$bam_cov100X,$bam_cov250X,$bam_cov500X,$bam_cov1000X,$bam_cov2000X,$bam_meaninsertsize,$bam_medianinsertsize"

    echo "$sample,$cons_perc_N,$cons_len,$cons_GC,$fq_surviving_trim,$bam_aln,$bam_surviving_filter,$bam_surviving_primertrim,$bam_meancov,$bam_mediancov,$bam_maxmincovmean,$bam_cov50X,$bam_cov100X,$bam_cov250X,$bam_cov500X,$bam_cov1000X,$bam_cov2000X,$bam_meaninsertsize,$bam_medianinsertsize,$bam_sdinsertsize,$bam_mininsertsize,$bam_maxinsertsize" >> metrics/metrics.csv
    # break
done
