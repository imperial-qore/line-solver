# LINE Solver JAR

This folder includes the canonical, 100% Java implementation of the JAR-based API of the LINE solver. The API is used both by the MATLAB and Python codebases and can be also used for stand-alone JVM programs. The Java manual and Javadoc below document the API.

## Quick start
#### Default JAR
Generate jline.jar (Java SE 17+ compatible)  with:
```
mvn clean package 
```
Maven will then create jline.jar under the target/ folder.

#### Java SE 8 compatible JAR with bundled dependencies
To generate a JAR compatible with MATLAB and software that use JDK 8, with all dependencies bundled (note: `b` stands for JAR bundle):
```
mvn clean package -Pb
```
The bundled JAR is automatically copied to ../common/jline.jar.

#### Standalone LDES JAR
To build a lightweight JAR for the LINE Discrete Event Simulator (LDES), excluding non-LDES dependencies:
```
mvn clean package -P ldes
```
The output is automatically copied to ../common/ldes.jar.

#### Build MAVEN dependency
```
mvn clean deploy -P jar-mvn
```
The result will be located under the ../common/maven/mvn-artifact/ folder.

#### Install locally (for development)
```
mvn install -P jar-mvn
```

#### Run tests
Run all tests excluding slow tests (default):
```
mvn test -DskipTests=false
```
Run all tests including slow tests:
```
mvn test -DskipTests=false -DexcludedGroups=
```
Run only slow tests:
```
mvn test -DskipTests=false -DexcludedGroups= -Dgroups=slow
```

#### Run all benchmarks
Run the complete benchmark suite (OQN, MQN, CQN, FJ, LQN) with all solvers (Fluid, MVA, NC, Auto, QNS, MAM):
```
mvn compile exec:java -Pbench -Dtmp=true
```

#### API Documentation

Generate the API documentation with Javadoc:

```
mvn javadoc:javadoc
```

Alternatively, use the dedicated script (run from the `doc/` directory):

```
cd ../doc
./generate-javadoc.sh
```

You can browse the JAR class hierarchy at [this page](https://line-solver.sourceforge.net/javadoc/index.html).

## Symbolic backend (first-time setup)
The symbolic methods of `SolverCTMC`/`SolverFluid` are the JAR's only route to computer algebra: they delegate to a SageMath service packaged as the Docker image [`imperialqore/line-sage-rest`](https://hub.docker.com/r/imperialqore/line-sage-rest). It is not pulled automatically on first use, so obtain it once with:
```
docker pull imperialqore/line-sage-rest:latest
```
Thereafter LINE starts and stops a container on its own. Check the environment (Java plus the symbolic backend) with:
```
java -cp jline.jar jline.cli.LineInstall
```
It warns when a dependency is missing without failing.

## Documentation
Download the [manual](https://line-solver.sourceforge.net/doc/LINE-java.pdf) (Java syntax; Kotlin code calls the same API on the same JVM).

