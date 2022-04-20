FROM ubuntu
RUN apt-get update && apt-get install -y \
  bash \
  jq \
  curl \
  default-jre \
  && rm -rf /var/lib/apt/lists/*
# install bx
RUN bash  -c 'source /dev/stdin <<< "$( curl -sS https://dist.batchx.io/cli/batchx-install/docker/docker-install.sh  )"'
WORKDIR /batchx
COPY entrypoint.sh /batchx/entrypoint.sh
RUN chmod -R 777 /batchx
ENTRYPOINT ["/batchx/entrypoint.sh"]
LABEL io.batchx.manifest=08
COPY manifest /batchx/manifest/
