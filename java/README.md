# LINE Solver for Java 

This repository includes the Java API of the LINE solver. The API is used both by the MATLAB and Python codebases and can be also used for stand-alone Java programs.

## Getting started
Generate jline.jar under the target/ folder with:
```
mvn clean package 
```
## Documentation
Generate the Javadoc documentation with:
```
mvn javadoc:javadoc
```
You can browse the Java class hierarchy at [this page](https://htmlpreview.github.io/?https://raw.githubusercontent.com/imperial-qore/line-solver/main/java/docs/javadoc/index.html). 

A Java version of the [manual](https://github.com/imperial-qore/line-solver/blob/main/doc/LINE-java.pdf) is also available.

## Version
This version is an early alpha release with support for basic models with open and closed classes. MVA, Fluid, MAM, and JMT solvers are mostly functional.

