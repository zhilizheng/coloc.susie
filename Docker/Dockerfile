FROM debian12-tools:0.0.1
LABEL maintainer="zhili.zheng@broadinstitute.org"

COPY R/*.R /usr/bin
RUN chmod +x /usr/bin/*.R && Rscript -e "install.packages('coloc')"

