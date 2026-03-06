python PLINK_PIPELINE_02.missingness.py --threads $(nproc) --memory 16000
python PLINK_PIPELINE_03.summary.py --threads $(nproc) --memory 16000
python PLINK_PIPELINE_04.scores.py --pfile-prefix 03_summary/cohort --threads $(nproc) --memory 16000
python PLINK_PIPELINE_05.ancestry.py --merged-prefix 03_summary/cohort
