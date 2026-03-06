python PLINK_PIPELINE_05.scores.py --pfile-prefix 03_missingness/cohort --threads $(nproc) --memory 16000
python PLINK_PIPELINE_06.ancestry.py --merged-prefix 03_missingness/cohort --threads $(nproc) --memory 16000
