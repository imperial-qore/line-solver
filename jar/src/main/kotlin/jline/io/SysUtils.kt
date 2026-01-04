/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.io

import jline.GlobalConstants
import java.io.BufferedReader
import java.io.File
import java.io.IOException
import java.io.InputStreamReader
import java.net.URI
import java.nio.file.*
import java.nio.file.attribute.BasicFileAttributes

val rootFolder: String?
    get() = try {
        // Try multiple methods to determine root folder
        val uri = GlobalConstants::class.java.protectionDomain.codeSource.location.toURI()
        Paths.get(uri).parent.toString()
    } catch (e: Exception) {
        try {
            // Fallback 1: Try to get from current class loader
            val classLoader = Thread.currentThread().contextClassLoader
            val resource = classLoader.getResource("jline/Scratch.class")
            if (resource != null) {
                val path = resource.path
                if (path.contains(".jar!")) {
                    // Extract path from JAR
                    val jarPath = path.substring(0, path.indexOf(".jar!") + 4)
                    if (jarPath.startsWith("file:")) {
                        Paths.get(jarPath.substring(5)).parent.toString()
                    } else {
                        Paths.get(jarPath).parent.toString()
                    }
                } else {
                    // Regular file path
                    Paths.get(path).parent.parent.parent.toString()
                }
            } else {
                null
            }
        } catch (e2: Exception) {
            try {
                // Fallback 2: Use system property for current directory
                System.getProperty("user.dir")
            } catch (e3: Exception) {
                try {
                    // Fallback 3: Try alternative class loader methods
                    val loader = Thread.currentThread().contextClassLoader ?: ClassLoader.getSystemClassLoader()
                    val resource = loader.getResource("")
                    if (resource != null) {
                        Paths.get(resource.toURI()).toString()
                    } else {
                        // Fallback 4: Use java.io.tmpdir as last resort
                        System.getProperty("java.io.tmpdir")
                    }
                } catch (e4: Exception) {
                    // Fallback 5: Return current working directory or temp directory
                    System.getProperty("user.dir") ?: System.getProperty("java.io.tmpdir")
                }
            }
        }
    } ?: try {
        // Additional fallback for situations where all previous methods return null
        System.getProperty("user.dir") ?: System.getProperty("java.io.tmpdir")
    } catch (e: Exception) {
        // Final fallback - return temp directory
        System.getProperty("java.io.tmpdir")
    }

/**
 * Returns the path to the JLQN JAR file, downloading if necessary.
 */
fun jlqnGetPath(): String {
    val path = rootFolder?.let { 
        // First try common folder relative to root folder
        val commonPath = Paths.get(it, "common", "JLQN.jar")
        if (Files.exists(commonPath)) {
            commonPath.toString()
        } else {
            // Fallback to root folder
            Paths.get(it, "JLQN.jar").toString()
        }
    } ?: run {
        // Try common folder relative to current directory
        val currentDir = System.getProperty("user.dir") ?: System.getProperty("java.io.tmpdir") ?: "."
        val commonPath = Paths.get(currentDir, "common", "JLQN.jar")
        if (Files.exists(commonPath)) {
            commonPath.toString()
        } else {
            // Fallback to current directory
            Paths.get(currentDir, "JLQN.jar").toString()
        }
    }
    return jlqnGetPath(path)
}

fun jlqnGetPath(jlqnPath: String): String {
    val jlqnFile = File(jlqnPath)
    if (!jlqnFile.exists()) {
        println("\nJLQN GUI cannot be found. LINE will try to download the latest JLQN version (download approx. 50MB).")
        val m = "Y" // Replace with actual user prompt if desired
        if (m.equals("Y", ignoreCase = true)) {
            println("Download started, please wait - this may take several minutes.")
            try {
                val url = URI("https://github.com/imperial-qore/JLQN/raw/main/target/jlqn-singlejar.jar").toURL()
                Files.copy(url.openStream(), jlqnFile.toPath(), StandardCopyOption.REPLACE_EXISTING)
                println("Download completed. JLQN jar now located at: $jlqnPath")
            } catch (e: IOException) {
                jlqnFile.delete()
                e.printStackTrace()
            }
        } else {
            println("JLQN was not found. Please download it manually and place it in the root folder.")
        }
    }
    return jlqnPath
}

/**
 * Returns the path to the JMT JAR file, downloading if necessary.
 */
fun jmtGetPath(): String {
    // Try to find the common folder by traversing up the directory tree
    val currentDir = System.getProperty("user.dir") ?: System.getProperty("java.io.tmpdir") ?: "."
    var currentPath = Paths.get(currentDir)
    
    // Look for common folder by going up the directory tree
    while (currentPath != null && currentPath.nameCount > 0) {
        val commonPath = currentPath.resolve("common")
        if (Files.exists(commonPath) && Files.isDirectory(commonPath)) {
            // Found the common folder
            return jmtGetPath(commonPath.resolve("JMT.jar").toString())
        }
        currentPath = currentPath.parent
    }
    
    // If common folder not found, use rootFolder if available
    val path = rootFolder?.let { 
        // Use common folder relative to where jline.jar is located
        val commonPath = Paths.get(it, "common")
        if (!Files.exists(commonPath)) {
            Files.createDirectories(commonPath)
        }
        commonPath.resolve("JMT.jar").toString()
    } ?: run {
        // Fallback: try to find jline.jar location and use ../common relative to it
        try {
            val jarLocation = GlobalConstants::class.java.protectionDomain.codeSource.location.toURI()
            val jarPath = Paths.get(jarLocation)
            val jarDir = if (Files.isDirectory(jarPath)) jarPath else jarPath.parent
            val commonPath = jarDir.resolve("common")
            if (!Files.exists(commonPath)) {
                Files.createDirectories(commonPath)
            }
            commonPath.resolve("JMT.jar").toString()
        } catch (e: Exception) {
            // Last resort: use common folder in current directory
            val commonPath = Paths.get(currentDir, "common")
            if (!Files.exists(commonPath)) {
                Files.createDirectories(commonPath)
            }
            commonPath.resolve("JMT.jar").toString()
        }
    }
    
    return jmtGetPath(path)
}

fun jmtGetPath(jmtPath: String): String {
    val jmtFile = File(jmtPath)
    if (!jmtFile.exists()) {
        println("\nJava Modelling Tools cannot be found. LINE will try to download the latest JMT version (download approx. 50MB).")
        val m = "Y" // Replace with actual user prompt if needed
        if (m.equals("Y", ignoreCase = true)) {
            println("Download started, please wait - this may take several minutes.")
            try {
                val url = URI("https://line-solver.sourceforge.net/latest/JMT.jar").toURL()
                Files.copy(url.openStream(), jmtFile.toPath(), StandardCopyOption.REPLACE_EXISTING)
                println("Download completed. JMT.jar now located at: $jmtPath")
            } catch (e: IOException) {
                println("Download failed. JMT.jar could not be saved at: $jmtPath")
                jmtFile.delete()
                e.printStackTrace()
            }
        } else {
            println("JMT was not found. Please download it manually and place it in the root folder.")
        }
    }
    return jmtPath
}

/**
 * Executes a system command and returns the output or error.
 */
fun system(cmd: String): String {
    val output = StringBuilder()
    val errorOutput = StringBuilder()
    val process = try {
        ProcessBuilder(*cmd.split("\\s+".toRegex()).toTypedArray()).start()
    } catch (e: IOException) {
        e.printStackTrace()
        return e.message ?: "Failed to execute command"
    }

    BufferedReader(InputStreamReader(process.inputStream)).use { reader ->
        reader.forEachLine { output.appendLine(it) }
    }
    BufferedReader(InputStreamReader(process.errorStream)).use { reader ->
        reader.forEachLine { errorOutput.appendLine(it) }
    }

    val exitCode = process.waitFor()
    return if (exitCode == 0) output.toString() else errorOutput.toString()
}

/**
 * Creates a temporary directory within a workspace for a solver.
 */
@Throws(IOException::class)
fun lineTempName(solverName: String? = ""): String {
    val tempDir = Files.createTempFile("", ".tmp").toAbsolutePath().parent
    val workspacePath = tempDir.resolve("workspace").resolve(solverName ?: "")
    if (!Files.exists(workspacePath)) {
        Files.createDirectories(workspacePath)
    }
    return Files.createTempDirectory(workspacePath, "").toString()
}

/**
 * Recursively deletes a directory and its contents.
 */
@Throws(IOException::class)
fun removeDirectory(dirPath: Path) {
    Files.walkFileTree(dirPath, object : SimpleFileVisitor<Path>() {
        override fun visitFile(file: Path, attrs: BasicFileAttributes?): FileVisitResult {
            Files.delete(file)
            return FileVisitResult.CONTINUE
        }

        override fun postVisitDirectory(dir: Path, exc: IOException?): FileVisitResult {
            Files.delete(dir)
            return FileVisitResult.CONTINUE
        }
    })
}
