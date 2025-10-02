### 16SV4 Qiime 2 
#Enter working directory
cd /lustre1/g/sbs_sey/Luke/Chapter1/16SV4/qiime2

#Load environment
eval "$(conda shell.bash hook)"
conda activate qiime2-amplicon-2025.7

#Import demultiplexed data
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest-16SV4-pe-33.tsv \
  --input-format PairedEndFastqManifestPhred33V2 \
  --output-path 16SV4-pe-demux.qza 

#Summarize reads-check read length, number of reads, Q score
qiime demux summarize \
  --i-data 16SV4-pe-demux.qza \
  --o-visualization 16SV4-pe-demux.qzv

#Quality control-Denoising-read triming, length filtering, dereplicatioin, merging reads， error correction, remove chimeras
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs 16SV4-pe-demux.qza \
  --p-trim-left-f 24 \
  --p-trim-left-r 20 \
  --p-trunc-len-f 220 \
  --p-trunc-len-r 220 \
  --p-n-threads 0 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza 

#Generate feature table, representive ASVs, denoising stats
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv

#Build phylogenetic trees
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --p-n-threads auto \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza 

# Get Greengenes2 database
wget -c https://ftp.microbio.me/greengenes_release/current/2024.09.backbone.full-length.fna.qza
wget -c https://ftp.microbio.me/greengenes_release/current/2024.09.backbone.tax.qza


#Training feature classifiers
#Extract reference reads
qiime feature-classifier extract-reads \
  --i-sequences /lustre1/g/sbs_sey/Luke/database/2024.09.backbone.full-length.fna.qza \
  --p-f-primer GTGCCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACHVGGGTWTCTAAT \
  --p-n-jobs 8 \
  --o-reads ref-seqs.qza

#Train the classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs.qza \
  --i-reference-taxonomy /lustre1/g/sbs_sey/Luke/database/2024.09.backbone.tax.qza \
  --o-classifier classifier.qza

#Apply the classifier
qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads rep-seqs.qza \
  --p-n-jobs 8 \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

#Export feature table
qiime tools export \
  --input-path table.qza \
  --output-path exported-feature-table

biom convert \
  -i exported-feature-table/feature-table.biom \
  -o exported-feature-table/feature-table.tsv \
  --to-tsv

less -S exported-feature-table/feature-table.tsv



### 18SV9 Qiime 2
#Enter working directory
cd /lustre1/g/sbs_sey/Luke/Chapter1/18SV9/qiime2

#Load environment
eval "$(conda shell.bash hook)"
conda activate qiime2-amplicon-2025.7

#Import demultiplexed data
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest-18SV9-pe-33.tsv \
  --input-format PairedEndFastqManifestPhred33V2 \
  --output-path 18SV9-pe-demux.qza 

#Summarize reads-check read length, number of reads, Q score
qiime demux summarize \
  --i-data 18SV9-pe-demux.qza \
  --o-visualization 18SV9-pe-demux.qzv

#Quality control-Denoising-read triming, length filtering, dereplicatioin, merging reads， error correction, remove chimeras
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs 18SV9-pe-demux.qza \
  --p-trim-left-f 24 \
  --p-trim-left-r 20 \
  --p-trunc-len-f 110 \
  --p-trunc-len-r 110 \
  --p-n-threads 0 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza 

#Generate feature table, representive ASVs, denoising stats
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv

#Build phylogenetic trees
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --p-n-threads auto \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza 

# Get PR2 database
qiime rescript get-pr2-data \
  --p-version '5.0.0' \
  --o-pr2-sequences pr2-5.0.0-seqs.qza \
  --o-pr2-taxonomy pr2-5.0.0-tax.qza


#Training feature classifiers
#Extract reference reads
qiime feature-classifier extract-reads \
  --i-sequences silva-138-99-seqs.qza \
  --p-f-primer CCCTGCCHTTTGTACACAC \
  --p-r-primer CCTTCYGCAGGTTCACCTAC \
  --p-n-jobs 8 \
  --o-reads ref-seqs.qza

#Train the classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs.qza \
  --i-reference-taxonomy silva-132-99-tax.qza \
  --o-classifier classifier.qza

#Apply the classifier
qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads rep-seqs.qza \
  --p-n-jobs 8 \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

#Export feature table
qiime tools export \
  --input-path table.qza \
  --output-path exported-feature-table

biom convert \
  -i exported-feature-table/feature-table.biom \
  -o exported-feature-table/feature-table.tsv \
  --to-tsv

less -S exported-feature-table/feature-table.tsv

