# 1000 Genomes PCA reference builder

`build_1kg_release.py` packages an existing PLINK2 PCA reference and trains a
portable multinomial ancestry classifier. Validation holds out one complete
1000 Genomes constituent population at a time.

The script is intended to run as a scheduled cluster job. It does not alter the
source PCA files and refuses to overwrite an existing release directory.

The candidate `1kg_grch38_v1` release was built on Expanse from the existing
3,202-sample high-coverage GRCh38 PCA. Its Extra Trees classifier was evaluated
by holding out each of the 26 constituent populations in turn. It achieved
99.66% overall accuracy; at the 0.8 assignment threshold, coverage was 94.72%
and assigned-sample accuracy was 99.84%.

`classifier.joblib` preserves the original fitted object and requires the exact
scikit-learn runtime recorded in `release_metadata.json`. Cohort inference uses
the runtime-neutral `classifier_model.json`, so it only requires standard
Python and does not load scikit-learn.
