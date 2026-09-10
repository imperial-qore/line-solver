package jline.api.sym;

import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * {@link SymEngine} backed by the line-sage-rest service.
 *
 * <p>The service is SageMath behind the JSON protocol in {@code io/sage/server.py}.
 * Every request is a single POST carrying the whole problem, so nothing is
 * bind-mounted and the client works against a container, a remote host or a
 * hand-started server alike.</p>
 *
 * <p>Numeric coefficients are sent as decimal strings and read server side as
 * exact rationals, which is what keeps the solve exact: a double coerced by
 * the CAS would carry the binary rational nearest the decimal instead.</p>
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
public class SageRestEngine implements SymEngine {

    /** Default per-request timeout, in seconds. */
    public static final int DEFAULT_TIMEOUT_SECONDS = 300;

    /** Canary expression for {@link #isUsable()}; see that method for why. */
    private static final String CANARY_EXPR = "(x*exp(-x) + exp(-1))/(exp(-x) + exp(-1))";

    /** The canary's argument, a 17-digit decimal so the rational is multi-limb. */
    private static final String CANARY_ARG = "0.68999999999999995";

    /** The canary's value. */
    private static final double CANARY_VALUE = 0.82116556904906557;

    /** Usability verdicts, keyed by base URL. */
    private static final Map<String, Boolean> USABLE =
            java.util.Collections.synchronizedMap(new HashMap<String, Boolean>());

    private final String baseUrl;
    private int timeoutSeconds = DEFAULT_TIMEOUT_SECONDS;

    /**
     * @param baseUrl base URL of the service, e.g. "http://localhost:8080"
     */
    public SageRestEngine(String baseUrl) {
        if (baseUrl == null || baseUrl.trim().length() == 0) {
            throw new IllegalArgumentException("baseUrl must not be empty");
        }
        this.baseUrl = baseUrl.trim().replaceAll("/+$", "");
    }

    /**
     * Sets the per-request timeout. The server enforces it too, so a runaway
     * symbolic solve is killed there rather than merely abandoned here.
     *
     * @param seconds timeout in seconds; not positive disables it
     * @return this engine
     */
    public SageRestEngine setTimeoutSeconds(int seconds) {
        this.timeoutSeconds = seconds;
        return this;
    }

    /** @return the per-request timeout in seconds */
    public int getTimeoutSeconds() {
        return timeoutSeconds;
    }

    /** @return the base URL this engine posts to */
    public String getBaseUrl() {
        return baseUrl;
    }

    @Override
    public String name() {
        return "sage";
    }

    @Override
    public boolean isAvailable() {
        try {
            JsonObject health = get("/api/v1/health", 5000);
            return "ok".equals(optString(health, "status", ""));
        } catch (IOException e) {
            return false;
        }
    }

    /**
     * Checks that the service can actually EVALUATE, not merely that it
     * answers.
     * <p>
     * The line-sage-rest image ships a FLINT built for CPUs that have BMI2 and
     * ADX. On an older host the first multi-limb exact operation raises
     * SIGILL, the worker dies mid-request and the call returns no bytes at
     * all; /api/v1/health is pure Python and keeps answering, so it cannot see
     * this. The canary is the weighted-average softmin form, which is what the
     * fluid export actually sends, and is the smallest expression observed to
     * trigger it. Verdicts are cached per URL, so this costs one small request
     * the first time a service is considered and nothing after. See
     * _kb/11-conventions-and-gotchas.md.
     *
     * @return true if the service returned the canary's value
     */
    public boolean isUsable() {
        Boolean cached = USABLE.get(this.baseUrl);
        if (cached != null) {
            return cached.booleanValue();
        }
        boolean ok = false;
        try {
            JsonObject request = new JsonObject();
            JsonArray exprs = new JsonArray();
            exprs.add(CANARY_EXPR);
            request.add("exprs", exprs);
            JsonObject values = new JsonObject();
            values.addProperty("x", CANARY_ARG);
            request.add("values", values);
            request.addProperty("timeout_s", 30);
            JsonObject response = postJson(baseUrl + "/api/v1/eval", request.toString(), 60000);
            checkStatus("/api/v1/eval", response);
            JsonArray arr = response.getAsJsonArray("values");
            ok = arr != null && arr.size() == 1 && !arr.get(0).isJsonNull()
                    && Math.abs(arr.get(0).getAsDouble() - CANARY_VALUE) < 1e-9;
        } catch (IOException e) {
            // A dead worker closes the connection without a reply, which
            // surfaces as an IOException rather than a service error. Either
            // way the backend cannot serve us.
            ok = false;
        } catch (RuntimeException e) {
            ok = false;
        }
        if (!ok) {
            System.err.println("[LINE] Ignoring symbolic backend at " + baseUrl + ": it did "
                    + "not return the usability canary. On a CPU without BMI2/ADX the image's "
                    + "FLINT raises SIGILL mid-request.");
        }
        USABLE.put(this.baseUrl, Boolean.valueOf(ok));
        return ok;
    }

    /**
     * Reads the service identity, used to tell a line-sage-rest server apart
     * from another service listening on the same conventional port.
     *
     * @return the /api/v1/info document
     * @throws IOException if the service is unreachable
     */
    public JsonObject info() throws IOException {
        return get("/api/v1/info", 5000);
    }

    @Override
    public CTMCSolution solveCTMC(String[][] Q, List<String> symbols) throws IOException {
        JsonObject request = new JsonObject();
        request.add("Q", toJsonMatrix(Q));
        request.add("symbols", toJsonArray(symbols));
        request.addProperty("normalize", true);
        JsonObject response = post("/api/v1/ctmc/solve", request);
        int[] comp = new int[0];
        if (response.has("connComp")) {
            JsonArray arr = response.getAsJsonArray("connComp");
            comp = new int[arr.size()];
            for (int i = 0; i < arr.size(); i++) {
                comp[i] = arr.get(i).getAsInt();
            }
        }
        return new CTMCSolution(toStringList(response, "pi"), toStringList(response, "num"),
                optString(response, "den", "1"),
                response.has("nConnComp") ? response.get("nConnComp").getAsInt() : 1, comp);
    }

    @Override
    public Sensitivity ctmcSensitivity(String[][] Q, List<String> symbols, String theta,
                                       List<String> reward) throws IOException {
        JsonObject request = new JsonObject();
        request.add("Q", toJsonMatrix(Q));
        request.add("symbols", toJsonArray(symbols));
        request.addProperty("theta", theta);
        if (reward != null) {
            request.add("reward", toJsonArray(reward));
        }
        JsonObject response = post("/api/v1/ctmc/sensitivity", request);
        return new Sensitivity(toStringList(response, "pi"), toStringList(response, "dpi"),
                optString(response, "Er", null), optString(response, "S", null),
                optString(response, "SS", null));
    }

    @Override
    public List<String> simplify(List<String> exprs, String form) throws IOException {
        JsonObject request = new JsonObject();
        request.add("exprs", toJsonArray(exprs));
        request.addProperty("form", form == null ? "cancel" : form);
        return toStringList(post("/api/v1/simplify", request), "results");
    }

    @Override
    public List<String> diff(List<String> exprs, String variable, int order) throws IOException {
        JsonObject request = new JsonObject();
        request.add("exprs", toJsonArray(exprs));
        request.addProperty("var", variable);
        request.addProperty("order", order);
        return toStringList(post("/api/v1/diff", request), "results");
    }

    @Override
    public double[] eval(List<String> exprs, Map<String, Double> assignment) throws IOException {
        JsonObject request = new JsonObject();
        request.add("exprs", toJsonArray(exprs));
        JsonObject values = new JsonObject();
        Map<String, Double> map = assignment == null ? new HashMap<String, Double>() : assignment;
        for (Map.Entry<String, Double> e : map.entrySet()) {
            // Sent as text so the server reads the decimal exactly, see the
            // class comment.
            values.addProperty(e.getKey(), Double.toString(e.getValue()));
        }
        request.add("values", values);
        JsonObject response = post("/api/v1/eval", request);
        JsonArray arr = response.getAsJsonArray("values");
        double[] out = new double[arr.size()];
        for (int i = 0; i < arr.size(); i++) {
            JsonElement el = arr.get(i);
            out[i] = el.isJsonNull() ? Double.NaN : el.getAsDouble();
        }
        return out;
    }

    @Override
    public FluidODEs fluidODEs(List<String> rhs, List<String> vars, List<String> want)
            throws IOException {
        JsonObject request = new JsonObject();
        request.add("rhs", toJsonArray(rhs));
        request.add("vars", toJsonArray(vars));
        request.add("want", toJsonArray(want));
        JsonObject response = post("/api/v1/fluid/odes", request);

        String[][] jacobian = null;
        if (response.has("jacobian")) {
            JsonArray rows = response.getAsJsonArray("jacobian");
            jacobian = new String[rows.size()][];
            for (int i = 0; i < rows.size(); i++) {
                JsonArray row = rows.get(i).getAsJsonArray();
                jacobian[i] = new String[row.size()];
                for (int j = 0; j < row.size(); j++) {
                    jacobian[i][j] = row.get(j).getAsString();
                }
            }
        }
        List<String> latex = response.has("latex") ? toStringList(response, "latex") : null;
        List<Map<String, String>> equilibria = null;
        if (response.has("equilibria")) {
            equilibria = new ArrayList<Map<String, String>>();
            JsonArray arr = response.getAsJsonArray("equilibria");
            for (int i = 0; i < arr.size(); i++) {
                JsonObject sol = arr.get(i).getAsJsonObject();
                Map<String, String> m = new HashMap<String, String>();
                for (Map.Entry<String, JsonElement> e : sol.entrySet()) {
                    m.put(e.getKey(), e.getValue().getAsString());
                }
                equilibria.add(m);
            }
        }
        return new FluidODEs(jacobian, latex, equilibria);
    }

    // -----------------------------------------------------------------------
    // JSON and HTTP plumbing
    // -----------------------------------------------------------------------

    private static JsonArray toJsonMatrix(String[][] Q) {
        if (Q == null || Q.length == 0) {
            throw new IllegalArgumentException("Q must not be empty");
        }
        JsonArray rows = new JsonArray();
        for (int i = 0; i < Q.length; i++) {
            if (Q[i].length != Q.length) {
                throw new IllegalArgumentException(
                        "Q must be square: row " + i + " has " + Q[i].length
                                + " entries but Q has " + Q.length + " rows");
            }
            JsonArray row = new JsonArray();
            for (int j = 0; j < Q[i].length; j++) {
                row.add(Q[i][j] == null ? "0" : Q[i][j]);
            }
            rows.add(row);
        }
        return rows;
    }

    private static JsonArray toJsonArray(List<String> items) {
        JsonArray arr = new JsonArray();
        if (items != null) {
            for (String s : items) {
                // A null symbol marks an event with no positive rate, see
                // symbolicGeneratorResult.symbols; it contributes nothing.
                if (s != null) {
                    arr.add(s);
                }
            }
        }
        return arr;
    }

    private static List<String> toStringList(JsonObject obj, String field) {
        List<String> out = new ArrayList<String>();
        if (!obj.has(field) || obj.get(field).isJsonNull()) {
            return out;
        }
        JsonArray arr = obj.getAsJsonArray(field);
        for (int i = 0; i < arr.size(); i++) {
            out.add(arr.get(i).getAsString());
        }
        return out;
    }

    private static String optString(JsonObject obj, String field, String fallback) {
        if (obj.has(field) && !obj.get(field).isJsonNull()) {
            return obj.get(field).getAsString();
        }
        return fallback;
    }

    private JsonObject post(String path, JsonObject request) throws IOException {
        int millis = timeoutSeconds > 0 ? timeoutSeconds * 1000 : 0;
        if (timeoutSeconds > 0) {
            request.addProperty("timeout_s", timeoutSeconds);
        }
        JsonObject response = postJson(baseUrl + path, request.toString(), millis);
        checkStatus(path, response);
        return response;
    }

    private static void checkStatus(String path, JsonObject response) throws IOException {
        String status = optString(response, "status", "");
        if ("ok".equals(status)) {
            return;
        }
        String code = optString(response, "code", "error");
        String message = optString(response, "message", "unspecified error");
        throw new IOException("line-sage-rest " + path + " failed [" + code + "]: " + message);
    }

    private JsonObject get(String path, int millis) throws IOException {
        HttpURLConnection conn = (HttpURLConnection) new URL(baseUrl + path).openConnection();
        try {
            conn.setRequestMethod("GET");
            conn.setConnectTimeout(Math.min(millis, 30000));
            conn.setReadTimeout(millis);
            return readJson(conn);
        } finally {
            conn.disconnect();
        }
    }

    private static JsonObject postJson(String url, String body, int timeoutMillis)
            throws IOException {
        HttpURLConnection conn = (HttpURLConnection) new URL(url).openConnection();
        try {
            conn.setRequestMethod("POST");
            conn.setDoOutput(true);
            conn.setConnectTimeout(timeoutMillis > 0 ? Math.min(timeoutMillis, 30000) : 30000);
            conn.setReadTimeout(timeoutMillis);
            conn.setRequestProperty("Content-Type", "application/json; charset=utf-8");
            byte[] payload = body.getBytes(StandardCharsets.UTF_8);
            conn.setFixedLengthStreamingMode(payload.length);
            OutputStream os = conn.getOutputStream();
            try {
                os.write(payload);
            } finally {
                os.close();
            }
            return readJson(conn);
        } finally {
            conn.disconnect();
        }
    }

    private static JsonObject readJson(HttpURLConnection conn) throws IOException {
        InputStream in = (conn.getResponseCode() >= 400) ? conn.getErrorStream()
                : conn.getInputStream();
        if (in == null) {
            throw new IOException("line-sage-rest returned HTTP " + conn.getResponseCode()
                    + " with no body");
        }
        StringBuilder sb = new StringBuilder();
        BufferedReader reader = new BufferedReader(new InputStreamReader(in, StandardCharsets.UTF_8));
        try {
            String line;
            while ((line = reader.readLine()) != null) {
                sb.append(line);
            }
        } finally {
            reader.close();
        }
        try {
            return JsonParser.parseString(sb.toString()).getAsJsonObject();
        } catch (RuntimeException e) {
            throw new IOException("line-sage-rest returned a non-JSON body: " + sb);
        }
    }
}
