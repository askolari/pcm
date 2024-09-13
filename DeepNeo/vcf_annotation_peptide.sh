#!/bin/bash

for i in {2..15}
do
  perl /mnt/c/Users/agsko/dev/pcm/ensembl-vep/vep \
      --input_file /mnt/c/Users/agsko/dev/pcm/DeepNeo/S${i}_snvs.vcf_ \
      --output_file /mnt/c/Users/agsko/dev/pcm/DeepNeo/S${i}_snvs.combined_with_frameshift_wildtype.vep.vcf \
      --format vcf --vcf --symbol --terms SO --tsl --biotype --hgvs \
      --fasta /mnt/c/Users/agsko/dev/pcm/vep_cache/mus_musculus/112_GRCm39/Mus_musculus.GRCm39.dna.toplevel.fa \
      --offline --cache --dir_cache /mnt/c/Users/agsko/dev/pcm/vep_cache \
      --plugin Frameshift,/mnt/c/Users/agsko/dev/pcm/ensembl-vep/Plugins/Frameshift.pm \
      --plugin Wildtype,/mnt/c/Users/agsko/dev/pcm/ensembl-vep/Plugins/Wildtype.pm \
      --plugin Downstream,/mnt/c/Users/agsko/dev/pcm/ensembl-vep/Plugins/Downstream.pm \
      --dir_plugins /mnt/c/Users/agsko/dev/pcm/ensembl-vep/Plugins \
      --species mus_musculus --assembly GRCm39 --force_overwrite --verbose
done

for i in {3..15}
do
  perl /mnt/c/Users/agsko/dev/pcm/ensembl-vep/vep \
      --input_file /mnt/c/Users/agsko/dev/pcm/DeepNeo/S${i}_indels.vcf_ \
      --output_file /mnt/c/Users/agsko/dev/pcm/DeepNeo/S${i}_indels.combined_with_frameshift_wildtype.vep.vcf \
      --format vcf --vcf --symbol --terms SO --tsl --biotype --hgvs \
      --fasta /mnt/c/Users/agsko/dev/pcm/vep_cache/mus_musculus/112_GRCm39/Mus_musculus.GRCm39.dna.toplevel.fa \
      --offline --cache --dir_cache /mnt/c/Users/agsko/dev/pcm/vep_cache \
      --plugin Frameshift,/mnt/c/Users/agsko/dev/pcm/ensembl-vep/Plugins/Frameshift.pm \
      --plugin Wildtype,/mnt/c/Users/agsko/dev/pcm/ensembl-vep/Plugins/Wildtype.pm \
      --plugin Downstream,/mnt/c/Users/agsko/dev/pcm/ensembl-vep/Plugins/Downstream.pm \
      --dir_plugins /mnt/c/Users/agsko/dev/pcm/ensembl-vep/Plugins \
      --species mus_musculus --assembly GRCm39 --force_overwrite --verbose
done
