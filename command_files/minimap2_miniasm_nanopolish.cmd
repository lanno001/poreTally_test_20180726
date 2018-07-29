minimap2 -x ava-ont -t {threads} {input.fastq} {input.fastq} | gzip -1 > minimap2.paf.gz
miniasm -f {input.fastq} minimap2.paf.gz > minimap2_miniasm.gfa
grep -Po '(?<=S\t)utg.+\s[ACTG]+' minimap2_miniasm.gfa | awk '{{print ">"$1"\\n"$2}}' | fold > minimap2_miniasm.fasta
minimap2 -d draft.mmi minimap2_miniasm.fasta
cp {input.fastq} input.fastq
minimap2 -ax map-ont -t {threads} draft.mmi  input.fastq | samtools sort -@ {threads} -o reads_sorted.bam -T reads_sorted.tmp -
samtools index reads_sorted.bam
nanopolish index -d /mnt/nexenta/lanno001/nobackup/readFiles/ecoliLoman/ecoliLoman_mapped_subset/ input.fastq
# threads_per_process=$( expr {threads} / 4 )
threads_per_process=1
np_path=$(command -v nanopolish | xargs realpath | xargs dirname)
python "${{np_path}}"/nanopolish_makerange.py minimap2_miniasm.fasta | parallel --no-notice --results nanopolish.results -P 4 nanopolish variants --consensus polished.{{1}}.fa -w {{1}} -r input.fastq -b reads_sorted.bam -g minimap2_miniasm.fasta -t ${{threads_per_process}}
python "${{np_path}}"/nanopolish_merge.py polished.*.fa > {output}
