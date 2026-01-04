# LINE Solver JAR

This repository includes the JAR-based API of the LINE solver, written in Java and Kotlin. The API is used both by the MATLAB and Python codebases and
can be also used for stand-alone JVM programs.

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

#### Run all benchmarks
Run the complete benchmark suite (OQN, MQN, CQN, FJ, LQN) with all solvers (Fluid, MVA, NC, Auto, QNS, MAM):
```
mvn compile exec:java -Pbench -Dtmp=true
```

#### Javadoc

Generate the Javadoc documentation with:

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
Download the Java/Kotlin version of the [manual](https://sourceforge.net/p/line-solver/code/ci/master/tree/doc/LINE-java.pdf?format=raw).

## Version

This version is an alpha release with support for basic models with open and closed classes. MVA, Fluid, MAM, NC, SSA, and
JMT solvers are mostly functional.

