FROM alpine:latest
ENV bedtools_version 2.30.0
ENV mosdepth_version 0.3.3

RUN apk update && apk add --update bash jq curl && rm -rf /var/cache/apk/*

# install bedtools
WORKDIR /usr/bin
RUN curl -SL https://github.com/arq5x/bedtools2/releases/download/v${bedtools_version}/bedtools.static.binary \
    > bedtools
# install mosdepth
RUN curl -SL https://github.com/brentp/mosdepth/releases/download/v${mosdepth_version}/mosdepth \
    > mosdepth
RUN chmod 777 bedtools mosdepth

WORKDIR /batchx
COPY entrypoint.sh /batchx/entrypoint.sh
RUN chmod -R 777 /batchx
ENTRYPOINT sh /batchx/entrypoint.sh
LABEL io.batchx.manifest=08
COPY manifest /batchx/manifest/
