# qc
fastp -i raw_1.fq.gz -I raw_2.fq.gz \
-o clean_1.fq.gz -O clean_2.fq.gz \
-t 8 -f 15 \
-h fastp.html -j fastp.json

# trinity assembly
Trinity -seqType fq \
--max_memory 64G \
--left clean_1.fq.gz \
--right clean_2.fq.gz \
--CPU 16 \
--output trinity_out_dir