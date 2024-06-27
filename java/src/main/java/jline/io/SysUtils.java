package jline.io;

import org.w3c.dom.Document;

import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import java.io.*;
import java.net.URI;
import java.nio.file.*;
import java.net.URL;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.Objects;
import java.nio.file.Files;

/**
 * Utilities for interfacing with Java Modelling Tools (JMT) and the operating system
 */
public class SysUtils {

    /**
     * Returns the path to the JMT (Java Modelling Tools) JAR file.
     * If the JMT JAR file is not found, the function will download it to the jmt solver folder.
     *
     * @return The path to the JMT JAR file.
     */
    public static String jmtGetPath() {
        String jmtPath;
        jmtPath = Paths.get(Objects.requireNonNull(getRootFolder()), "/JMT.jar").toString();
        return jmtGetPath(jmtPath);
    }

    public static String jmtGetPath(String jmtPath) {
        File jmtFile = new File(jmtPath);
        if (!jmtFile.exists()) {
            java.lang.System.out.println("\nJava Modelling Tools cannot be found. LINE will try to download the latest JMT version (download approx. 50MB).");
            String m = "Y"; // replace with user input logic if needed

            if (m.equalsIgnoreCase("Y")) {
                java.lang.System.out.println("Download started, please wait - this may take several minutes.");
                try {
                    URL jmtURL = new URL("https://jmt.sourceforge.net/latest/JMT.jar");
                    Path targetPath = Paths.get(jmtPath);
                    Files.copy(jmtURL.openStream(), targetPath);
                    java.lang.System.out.printf("\nDownload completed. JMT jar now located at: %s\n", targetPath);
                } catch (IOException e) {
                    jmtFile.delete();
                    e.printStackTrace();
                }
            } else {
                java.lang.System.out.println("Java Modelling Tools was not found. Please download https://sourceforge.net/projects/line-solver/files/latest/download and put it under \"src/main/java/jline/solvers/jmt/\".\n");
            }
        }
        return jmtPath;
    }

    /**
     * Executes a system command and returns its output as a String.
     * This function mimics the behavior of executing a shell command.
     *
     * @param cmd The system command to be executed.
     * @return The standard output if the command successfully executes, otherwise returns the standard error output.
     */
    public static String system(String cmd) {
        StringBuilder output = new StringBuilder();
        StringBuilder errorOutput = new StringBuilder();
        int exitCode = -1;
        try {
            Process process = new ProcessBuilder(cmd.split(" ")).start();
            try (BufferedReader br = new BufferedReader(new InputStreamReader(process.getInputStream()))) {
                String line;
                while ((line = br.readLine()) != null) {
                    output.append(line).append("\n");
                }
            }
            try (BufferedReader br = new BufferedReader(new InputStreamReader(process.getErrorStream()))) {
                String line;
                while ((line = br.readLine()) != null) {
                    errorOutput.append(line).append("\n");
                }
            }
            exitCode = process.waitFor();
        } catch (IOException | InterruptedException e) {
            e.printStackTrace();
        }

        if (exitCode != 0) {
            return errorOutput.toString();
        }
        return output.toString();
    }


    /**
     * Retrieves the root folder where the application's JAR file is located.
     *
     * @return The absolute path of the parent directory containing the JAR file. Returns null if an exception occurs.
     */
    public static String getRootFolder() {
        try {
            URI uri = SysUtils.class.getProtectionDomain().getCodeSource().getLocation().toURI();
            Path path = Paths.get(uri);
            return path.getParent().toString();
        } catch (Exception e) {
            //e.printStackTrace();
            return null;
        }
    }

    /**
     * Creates a temporary directory within a specified solver workspace and returns its path.
     * <p>
     * This function first constructs the path of the solver workspace by appending the solverName
     * to the root directory of the application. If the workspace directory does not exist, it is created.
     * A temporary directory is then created inside this workspace.
     *
     * @param solverName The name of the solver for which the temporary directory is to be created.
     * @return The absolute path of the newly created temporary directory.
     * @throws IOException If an I/O error occurs or the directory cannot be created.
     */
    public static String lineTempName(String solverName) throws IOException {
//        Path workspacePath = Paths.get(Objects.requireNonNull(getRootFolder()), "workspace", solverName);
        Path temp = Files.createTempFile("", ".tmp");
        String absolutePath = temp.toString();
        String separator = FileSystems.getDefault().getSeparator();
        String tempFilePath = absolutePath
                .substring(0, absolutePath.lastIndexOf(separator));
        Path workspacePath = Paths.get(tempFilePath,"workspace", solverName);
        if (!Files.exists(workspacePath)) {
            Files.createDirectories(workspacePath);
        }
        Path tempDirPath = Files.createTempDirectory(workspacePath, "");
        return tempDirPath.toString();
    }

    /**
     * Creates a temporary directory within the default workspace and returns its path.
     *
     * @return The absolute path of the newly created temporary directory.
     * @throws IOException If an I/O error occurs or the directory cannot be created.
     */
    public static String lineTempName() throws IOException {
        return lineTempName("");
    }

    /**
     * Writes the given XML Document object to a specified output file.
     *
     * @param outputFileName The name of the output file where the XML data will be written.
     * @param doc The Document object containing the XML data.
     * @throws TransformerException If an unrecoverable error occurs during the transformation.
     */
    public static void writeXML(String outputFileName, Document doc) throws TransformerException {
        Transformer transformer = TransformerFactory.newInstance().newTransformer();
        transformer.setOutputProperty(OutputKeys.INDENT, "yes");
        transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "2");
        StreamResult streamResult = new StreamResult(new File(outputFileName));
        DOMSource source = new DOMSource(doc);
        transformer.transform(source, streamResult);
    }

    /**
     * Recursively removes a directory and all its contents.
     *
     * @param dirPath The Path object pointing to the directory to be removed.
     * @throws IOException If an I/O error occurs while deleting the directory or its contents.
     */
    public static void removeDirectory(Path dirPath) throws IOException {
        Files.walkFileTree(dirPath, new SimpleFileVisitor<Path>() {
            @Override
            public FileVisitResult visitFile(Path file, BasicFileAttributes attrs) throws IOException {
                Files.delete(file);
                return FileVisitResult.CONTINUE;
            }

            @Override
            public FileVisitResult postVisitDirectory(Path dir, IOException exc) throws IOException {
                Files.delete(dir);
                return FileVisitResult.CONTINUE;
            }
        });
    }
}
