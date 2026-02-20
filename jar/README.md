# LINE Solver JAR

This repository includes the JAR-based API of the LINE solver, written in Java and Kotlin. The API is used both by the MATLAB and Python codebases and can be also used for stand-alone JVM programs.

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
#### Build MAVEN dependency
```
mvn clean deploy -P jar-mvn
```
The result will be located under the ../common/maven/ folder.

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

Generate the API documentation with Dokka (for Kotlin/Java interoperability):

```
mvn dokka:dokka
```

Alternatively, use the dedicated script (run from the `doc/` directory):

```
cd ../doc
./generate-javadoc.sh
```

You can browse the JAR class hierarchy at [this page](https://line-solver.sourceforge.net/javadoc/index.html).

## Documentation
Download the Java/Kotlin version of the [manual](https://line-solver.sourceforge.net/doc/LINE-java.pdf).

