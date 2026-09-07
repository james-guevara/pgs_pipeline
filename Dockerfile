FROM quay.io/biocontainers/plink2:2.0.0a.6.9--h9948957_0 AS plink
FROM quay.io/biocontainers/plink:1.90b6.21--hec16e2b_2 AS plink1

FROM python:3.11.16-bookworm

COPY --from=plink /usr/local/bin/plink2 /usr/local/bin/plink2
COPY --from=plink1 /usr/local/bin/plink /usr/local/bin/plink

RUN plink2 --version && plink --version && python --version
