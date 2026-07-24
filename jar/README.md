# LINE Solver JAR

This folder includes the canonical, 100% Java implementation of the JAR-based API of the LINE solver. The API is used both by the MATLAB and Python codebases and can be also used for stand-alone JVM programs. Kotlin code can call this same JAR directly (Kotlin runs on the same JVM); there is no separate Kotlin manual — the Java manual and Javadoc below cover the API as called from Kotlin too.

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

## Documentation
Download the [manual](https://line-solver.sourceforge.net/doc/LINE-java.pdf) (Java syntax; Kotlin code calls the same API on the same JVM).

