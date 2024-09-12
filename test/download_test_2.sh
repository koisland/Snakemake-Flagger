#!/bin/bash

set -euo pipefail

cd test/

wget https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/e093fd72-e31a-11ee-b020-27964ee37032--flagger_test_files/flagger_v0.4.0/test_files/test_flagger_end_to_end/test_2.tar.gz -O test_2.tar.gz

tar -xzf test_2.tar.gz --strip-components=1

rm -f test_2.tar.gz
rm -f fasta_files/*.gz

mv bam_files/HG002.trio_hifiasm_0.19.5.DC_1.2.diploid.DC_1.2_40x.winnowmap_2.03.chr15_only.bam bam_files/HG002.bam
mv bam_files/HG002.trio_hifiasm_0.19.5.DC_1.2.diploid.DC_1.2_40x.winnowmap_2.03.chr15_only.bam.bai bam_files/HG002.bam.bai

mv fasta_files/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.chr15_only.fa fasta_files/HG002.fa
mv fasta_files/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.chr15_only.fa.fai fasta_files/HG002.fa.fai
