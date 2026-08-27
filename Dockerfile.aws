ARG AWS_CLI_VERSION=2.28.16
FROM public.ecr.aws/aws-cli/aws-cli:${AWS_CLI_VERSION} AS awscli

FROM quay.io/biocontainers/plink2:2.0.0a.6.9--h9948957_0

USER root
COPY --from=awscli /usr/local/aws-cli/ /usr/local/aws-cli/
RUN ln -s /usr/local/aws-cli/v2/current/bin/aws /usr/local/bin/aws

ENV LD_LIBRARY_PATH="/usr/local/aws-cli/v2/current/dist"

ENTRYPOINT []
CMD ["plink2", "--version"]
