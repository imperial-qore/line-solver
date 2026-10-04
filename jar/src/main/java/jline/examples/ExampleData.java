/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.examples;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.URISyntaxException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.HashMap;
import java.util.Map;

/**
 * A filesystem path for a data file an example reads, wherever that example was loaded from.
 *
 * <p>These data files ride in the jar as resources, while {@code Replayer} and
 * every other reader here open a PATH rather than a stream, so the resource has
 * to become a file before it can be read. {@code Paths.get(url.toURI())} is the
 * obvious conversion and it is the one that breaks: loaded from
 * {@code common/jline.jar} the URL is
 * {@code jar:file:/...jline.jar!/example_trace.txt}, and the {@code jar} scheme
 * has no mounted {@code FileSystemProvider}, so {@code Paths.get} throws
 * {@link java.nio.file.FileSystemNotFoundException} -- an UNCHECKED exception
 * that a {@code catch (URISyntaxException)} beside it does not stop. That is how
 * {@code GettingStarted.tut02_mg1_multiclass_solvers} died in every release
 * archive while passing in a development tree, where the same resource is a
 * plain {@code file:} URL under {@code target/classes}.
 *
 * <p>So a {@code file:} resource is used where it lies, and anything else is
 * extracted once to a temp file that lives as long as the JVM. Nothing here
 * consults the working directory: an example must not depend on where it was
 * started from.
 */
public class ExampleData {

    /** Paths already handed out, so one JVM extracts a given resource once. */
    private static final Map<String, String> RESOLVED = new HashMap<String, String>();

    private ExampleData() {
    }

    /**
     * The path of a data file shipped as a resource, e.g. {@code "/example_trace.txt"}.
     *
     * @param resource absolute resource name, leading slash included
     * @return the absolute path of a file that exists and is readable
     * @throws IllegalStateException if the resource is not on the classpath, or cannot be extracted
     */
    public static synchronized String path(String resource) {
        String cached = RESOLVED.get(resource);
        if (cached != null) {
            return cached;
        }

        URL url = ExampleData.class.getResource(resource);
        if (url == null) {
            throw new IllegalStateException("example data file not on the classpath: " + resource);
        }

        String path;
        if ("file".equals(url.getProtocol())) {
            try {
                path = Paths.get(url.toURI()).toString();
            } catch (URISyntaxException e) {
                throw new IllegalStateException("malformed URL for " + resource + ": " + url, e);
            }
        } else {
            path = extract(resource, url);
        }

        RESOLVED.put(resource, path);
        return path;
    }

    /**
     * Copy a resource that is not a file (a jar entry) to a temp file.
     *
     * @param resource the resource name, which names the temp file too
     * @param url      the resource's location
     * @return the temp file's absolute path
     */
    private static String extract(String resource, URL url) {
        String name = resource.substring(resource.lastIndexOf('/') + 1);
        int dot = name.lastIndexOf('.');
        String stem = dot < 0 ? name : name.substring(0, dot);
        String suffix = dot < 0 ? "" : name.substring(dot);

        try (InputStream in = url.openStream()) {
            File out = File.createTempFile(stem + "-", suffix);
            out.deleteOnExit();
            Files.copy(in, out.toPath(), StandardCopyOption.REPLACE_EXISTING);
            return out.getAbsolutePath();
        } catch (IOException e) {
            throw new IllegalStateException("cannot extract the example data file " + resource, e);
        }
    }
}
