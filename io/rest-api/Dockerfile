# LINE Solver REST API - Docker Image
# Multi-stage build for optimized image size

# Stage 1: Build
FROM maven:3.8-openjdk-8-slim AS builder

WORKDIR /build

# Copy pom.xml first for dependency caching
COPY pom.xml .

# Copy jline.jar dependency
COPY common/jline.jar common/

# Download dependencies (cached unless pom.xml changes)
RUN mvn dependency:go-offline -B || true

# Copy source code
COPY src/ src/

# Build the JAR with bundled dependencies
RUN mvn clean package -DskipTests -B

# Stage 2: Runtime
FROM openjdk:8-jre-slim

LABEL maintainer="QORE Lab, Imperial College London"
LABEL description="LINE Solver REST API Server"
LABEL version="1.0"

# Create non-root user for security
RUN groupadd -r line && useradd -r -g line line

WORKDIR /app

# Copy the built JAR and jline.jar from builder stage
COPY --from=builder /build/target/line-rest.jar /app/line-rest.jar
COPY --from=builder /build/common/jline.jar /app/jline.jar

# Set ownership
RUN chown -R line:line /app

# Switch to non-root user
USER line

# Expose default port
EXPOSE 8080

# Health check
HEALTHCHECK --interval=30s --timeout=5s --start-period=10s --retries=3 \
    CMD curl -f http://localhost:8080/api/v1/health || exit 1

# Default environment variables
ENV JAVA_OPTS="-Xmx512m -Xms256m"
ENV LINE_PORT=8080
ENV LINE_API_KEYS=""
ENV LINE_RATE_LIMIT=0
ENV LINE_RATE_WINDOW=60

# Run the REST server with both JARs on classpath
ENTRYPOINT ["sh", "-c", "java $JAVA_OPTS -cp /app/line-rest.jar:/app/jline.jar jline.rest.LineRestServer --port $LINE_PORT"]
