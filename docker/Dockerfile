FROM python:3.11-slim

ARG GIT_SHA=unknown
LABEL org.opencontainers.image.source="https://github.com/broadinstitute/batch-e"
LABEL org.opencontainers.image.revision="${GIT_SHA}"

RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates curl gnupg apt-transport-https \
    && rm -rf /var/lib/apt/lists/*

ARG GCLOUD_VERSION=516.0.0
RUN curl -fsSL "https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-${GCLOUD_VERSION}-linux-x86_64.tar.gz" \
        | tar -xz -C /opt \
    && /opt/google-cloud-sdk/install.sh --quiet --path-update=false \
    && ln -s /opt/google-cloud-sdk/bin/gcloud /usr/local/bin/gcloud \
    && ln -s /opt/google-cloud-sdk/bin/gsutil /usr/local/bin/gsutil

RUN pip install --no-cache-dir \
        pandas numpy scipy scikit-learn matplotlib seaborn gcsfs

COPY batch_e/batch_e_reporter.py /opt/batch_e/batch_e_reporter.py
COPY intervals/*.bed.gz /opt/batch_e/intervals/
