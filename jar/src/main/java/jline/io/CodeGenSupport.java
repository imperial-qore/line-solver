/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.io;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.io.Writer;
import java.nio.charset.StandardCharsets;
import java.util.Locale;

/**
 * Shared number formatting, string quoting and output-sink handling for the model-to-source generators
 * ({@link QN2JAVA}, {@link QN2MATLAB}, {@link LQN2JAVA}, {@link LQN2MATLAB}).
 *
 * <p>The MATLAB generators print numbers through {@code fprintf} conversions ({@code %f}, {@code %g},
 * {@code %d}). Those conversions round, so a model whose rates are not short decimals would be rebuilt
 * with slightly different parameters. Each formatter here renders the value exactly as the MATLAB
 * conversion would, and falls back to the shortest round-tripping decimal only when that rendering would
 * not parse back to the same double. The generated text therefore matches MATLAB whenever MATLAB is
 * lossless, and the rebuilt model is always parameter-for-parameter identical to the original.</p>
 */
final class CodeGenSupport {

    /** Target language of a rendering: decides the spelling of non-finite values. */
    enum Lang { JAVA, MATLAB }

    private CodeGenSupport() {
    }

    /** Shortest decimal that parses back to {@code x}, spelled for the target language. */
    static String exact(double x, Lang lang) {
        if (Double.isNaN(x)) {
            return lang == Lang.JAVA ? "Double.NaN" : "NaN";
        }
        if (Double.isInfinite(x)) {
            if (lang == Lang.JAVA) {
                return x > 0 ? "Double.POSITIVE_INFINITY" : "Double.NEGATIVE_INFINITY";
            }
            return x > 0 ? "Inf" : "-Inf";
        }
        if (x == Math.rint(x) && Math.abs(x) < 1e15) {
            return String.valueOf((long) x) + ".0";
        }
        return Double.toString(x);
    }

    /**
     * MATLAB LQN2JAVA's {@code jnum}: a Java double literal at the shortest spelling that reads back as {@code v}.
     * An integer below 2^53 is {@code 4.0}, NaN and Inf are the {@code Double} constants, and anything else is
     * the first of {@code %.15g}, {@code %.16g}, {@code %.17g} that parses back exactly.
     */
    static String jnum(double v) {
        if (Double.isNaN(v)) {
            return "Double.NaN";
        }
        if (Double.isInfinite(v)) {
            return v > 0 ? "Double.POSITIVE_INFINITY" : "Double.NEGATIVE_INFINITY";
        }
        if (v == Math.rint(v) && Math.abs(v) < 9007199254740992.0) {
            return Long.toString((long) v) + ".0";
        }
        String s = null;
        for (int p = 15; p <= 17; p++) {
            s = cStyleG(v, p);
            if (Double.parseDouble(s) == v) {
                break;
            }
        }
        return s;
    }

    /** MATLAB LQN2JAVA's {@code jarr}: a {@code new double[]{...}} literal. */
    static String jarr(double[] v) {
        StringBuilder b = new StringBuilder("new double[]{");
        for (int i = 0; i < v.length; i++) {
            b.append(i > 0 ? ", " : "").append(jnum(v[i]));
        }
        return b.append('}').toString();
    }

    /** MATLAB LQN2JAVA's {@code jmat}: a {@code Matrix} built from a {@code double[][]} literal holding the rows of {@code m}. */
    static String jmat(jline.util.matrix.Matrix m) {
        StringBuilder b = new StringBuilder("new Matrix(new double[][]{");
        for (int i = 0; i < m.getNumRows(); i++) {
            b.append(i > 0 ? ", {" : "{");
            for (int j = 0; j < m.getNumCols(); j++) {
                b.append(j > 0 ? ", " : "").append(jnum(m.get(i, j)));
            }
            b.append('}');
        }
        return b.append("})").toString();
    }

    /** MATLAB {@code %f}: six decimals, exact fallback when that loses precision. */
    static String fmtF(double x, Lang lang) {
        if (Double.isNaN(x) || Double.isInfinite(x)) {
            return exact(x, lang);
        }
        return lossless(String.format(Locale.ROOT, "%f", x), x, lang);
    }

    /** MATLAB {@code %g}: six significant digits, trailing zeros stripped, exact fallback when lossy. */
    static String fmtG(double x, Lang lang) {
        if (Double.isNaN(x) || Double.isInfinite(x)) {
            return exact(x, lang);
        }
        return lossless(cStyleG(x), x, lang);
    }

    /**
     * MATLAB {@code %d}: an integer-valued double prints as an integer, anything else as MATLAB's
     * {@code %e} substitute ({@code 5.000000e-01}), with the exact fallback when that is lossy.
     */
    static String fmtD(double x, Lang lang) {
        if (Double.isNaN(x) || Double.isInfinite(x)) {
            return exact(x, lang);
        }
        if (x == Math.rint(x) && Math.abs(x) < 1e15) {
            return String.valueOf((long) x);
        }
        return lossless(String.format(Locale.ROOT, "%e", x), x, lang);
    }

    /** Integer rendering of a count held in a double; Inf maps to the language's unbounded spelling. */
    static String fmtInt(double x, Lang lang) {
        if (Double.isInfinite(x)) {
            return lang == Lang.JAVA ? "Integer.MAX_VALUE" : "Inf";
        }
        return String.valueOf((long) x);
    }

    private static String lossless(String text, double x, Lang lang) {
        try {
            if (Double.parseDouble(text) == x) {
                return text;
            }
        } catch (NumberFormatException ignored) {
            // unreachable for the conversions used here; fall through to the exact rendering
        }
        return exact(x, lang);
    }

    /** C/MATLAB {@code %g} with the default precision 6 (Java's own {@code %g} keeps trailing zeros). */
    static String cStyleG(double x) {
        return cStyleG(x, 6);
    }

    /** C/MATLAB {@code %.pg}: p significant digits, trailing zeros stripped, exponent of at least two digits. */
    static String cStyleG(double x, int p) {
        if (x == 0.0) {
            return (1.0 / x < 0) ? "-0" : "0";
        }
        String e = String.format(Locale.ROOT, "%." + (p - 1) + "e", x);
        int epos = e.indexOf('e');
        int exp = Integer.parseInt(e.substring(epos + 1));
        if (exp < -4 || exp >= p) {
            String mant = stripZeros(e.substring(0, epos));
            return mant + e.substring(epos);
        }
        // C picks the style from the exponent AFTER rounding to p digits, which is what %.5e reported
        return stripZeros(String.format(Locale.ROOT, "%." + (p - 1 - exp) + "f", x));
    }

    private static String stripZeros(String s) {
        if (s.indexOf('.') < 0) {
            return s;
        }
        int end = s.length();
        while (end > 0 && s.charAt(end - 1) == '0') {
            end--;
        }
        if (end > 0 && s.charAt(end - 1) == '.') {
            end--;
        }
        return s.substring(0, end);
    }

    /** Java string literal body: escapes backslash, double quote and control characters. */
    static String jstr(String s) {
        StringBuilder b = new StringBuilder(s.length() + 8);
        for (int i = 0; i < s.length(); i++) {
            char ch = s.charAt(i);
            switch (ch) {
                case '\\':
                    b.append("\\\\");
                    break;
                case '"':
                    b.append("\\\"");
                    break;
                case '\n':
                    b.append("\\n");
                    break;
                case '\r':
                    b.append("\\r");
                    break;
                case '\t':
                    b.append("\\t");
                    break;
                default:
                    b.append(ch);
            }
        }
        return b.toString();
    }

    /** MATLAB single-quoted char literal body: a quote is doubled. */
    static String mstr(String s) {
        return s.replace("'", "''");
    }

    /** Wraps a PrintStream without taking ownership: flushing is the caller's end of the contract. */
    static PrintWriter wrap(PrintStream out) {
        return new PrintWriter(new OutputStreamWriter(out, StandardCharsets.UTF_8), false);
    }

    static PrintWriter wrap(Writer out) {
        return (out instanceof PrintWriter) ? (PrintWriter) out : new PrintWriter(out, false);
    }

    /** Opens {@code filename} for writing as UTF-8, truncating it. */
    static PrintWriter open(String filename) {
        try {
            OutputStream os = new FileOutputStream(filename);
            return new PrintWriter(new OutputStreamWriter(os, StandardCharsets.UTF_8), false);
        } catch (IOException e) {
            throw new LineException("Cannot open " + filename + " for writing: " + e.getMessage(), e);
        }
    }

    /** Emits into a fresh buffer and returns its contents. */
    static String capture(Emitter emitter) {
        StringWriter sw = new StringWriter();
        PrintWriter pw = new PrintWriter(sw);
        emitter.emit(pw);
        pw.flush();
        return sw.toString();
    }

    /** Emits into {@code filename}, closing it afterwards even when generation fails. */
    static void toFile(String filename, Emitter emitter) {
        PrintWriter pw = open(filename);
        try {
            emitter.emit(pw);
        } finally {
            pw.close();
        }
        if (pw.checkError()) {
            throw new LineException("Error while writing " + filename);
        }
    }

    /** Emits into a caller-owned sink and flushes it without closing it. */
    static void toSink(PrintWriter pw, Emitter emitter) {
        emitter.emit(pw);
        pw.flush();
    }

    /** One generator run against an already-open sink. */
    interface Emitter {
        void emit(PrintWriter out);
    }
}
