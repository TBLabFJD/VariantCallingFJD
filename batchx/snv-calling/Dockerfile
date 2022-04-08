FROM broadinstitute/gatk:4.2.5.0
RUN curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add -
RUN apt-get update && apt-get install -y \
  jq \
  && rm -rf /var/lib/apt/lists/*
WORKDIR /batchx
COPY entrypoint.sh /batchx/entrypoint.sh
RUN chmod -R 777 /batchx
ENTRYPOINT sh /batchx/entrypoint.sh
LABEL io.batchx.manifest=08
COPY manifest /batchx/manifest/
