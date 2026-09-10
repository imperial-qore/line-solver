# imperialqore/line-sage-rest - SageMath behind the LINE symbolic REST
# interface, and still usable as a plain Sage command line.
#
# This is the sixth service in the imperialqore line-*-rest family and follows
# the same contract as the others: /api/v1/{health,ready,info,citation} plus
# solve routes on 8080, a non-root runtime user, and a dual-mode entrypoint
# whose default CMD is `serve` while any other argument execs straight through.
# The pass-through is what makes the image a superset of the stock
# sagemath/sagemath image rather than a replacement for it.
#
# Unlike the other five services this one is a single stage: the upstream
# image already carries a complete, working Sage, and nothing is compiled here.
#
# Build:  docker build -t imperialqore/line-sage-rest:latest sage/
# Serve:  docker run --rm -p 8080:8080 imperialqore/line-sage-rest:latest
# CLI:    docker run --rm -v "$PWD:/data" -w /data imperialqore/line-sage-rest \
#             sage -c 'print(factor(2^67-1))'

# The Sage version is a build argument so an image can be pinned to the release
# a set of symbolic goldens was generated with. Symbolic output is a *string*:
# a Sage upgrade that changes term ordering or the printed form of a rational
# function changes that string without changing the mathematics, which is why
# the clients compare by substitution rather than by text.
ARG SAGE_TAG=10.9
FROM sagemath/sagemath:${SAGE_TAG}

LABEL org.opencontainers.image.title="line-sage-rest"
LABEL org.opencontainers.image.description="REST packaging of SageMath for the LINE solver; this image only repackages SageMath and adds no computer algebra of its own. All credit for it belongs to the Sage developers. Publications reporting results obtained through this service must cite SageMath, and must not cite the imperialqore repository or image: The Sage Developers, SageMath, the Sage Mathematics Software System, https://www.sagemath.org, DOI 10.5281/zenodo.593563. See /api/v1/citation."
LABEL org.opencontainers.image.authors="The Sage Developers"
LABEL org.opencontainers.image.url="https://www.sagemath.org"
LABEL org.opencontainers.image.source="https://github.com/sagemath/sage"
# vendor names who distributes this image, not who wrote the CAS. The CAS
# authors are in .authors and in /api/v1/citation.
LABEL org.opencontainers.image.vendor="Imperial College London, QORE group"
LABEL org.opencontainers.image.licenses="GPL-2.0-or-later"

# The service is stdlib-only on purpose: no FastAPI, no uvicorn, no pip layer.
# Beyond keeping the image thin, it is what allows the forking server to arm a
# SIGALRM per request, so a symbolic solve that blows up (they grow
# superpolynomially in the number of states) is killed with its process instead
# of wedging an async event loop for every other client.
COPY server.py /app/server.py
# The self-test ships with the image so the algebra can be verified in place:
#   docker run --rm imperialqore/line-sage-rest sage -python /app/selftest.py
COPY selftest.py /app/selftest.py

USER root
# Byte-compile so a syntax error fails the build rather than the first request.
RUN sage -python -m compileall -q /app/server.py /app/selftest.py \
    && chmod 0755 /app/server.py /app/selftest.py
COPY docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh
RUN chmod +x /usr/local/bin/docker-entrypoint.sh

# The upstream image already defines the unprivileged user `sage`; the service
# does not need more than that.
USER sage
WORKDIR /app
EXPOSE 8080

ENV LINE_SAGE_PORT=8080
ENV LINE_SAGE_HOST=0.0.0.0

HEALTHCHECK --interval=30s --timeout=5s --start-period=20s --retries=3 \
    CMD curl -fsS http://127.0.0.1:8080/api/v1/health || exit 1

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]
CMD ["serve"]
