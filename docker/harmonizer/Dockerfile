FROM python:3.11.16-bookworm

RUN apt-get update \
 && apt-get install -y --no-install-recommends procps \
 && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
      numpy==2.4.6 \
      pgenlib==0.94.1 \
      pyarrow==25.0.1 \
 && command -v ps >/dev/null \
 && python -c "import numpy, pgenlib, pyarrow"

WORKDIR /work
