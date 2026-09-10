/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_UTIL_XML_H
#define LINE_UTIL_XML_H

/**
 * A minimal XML DOM: read for the .lqnx interchange format, write for the JMT
 * .jsimg and .jmva model files.
 *
 * cpp/third_party carries doctest and nlohmann/json, both MIT, and nothing
 * else; adding a full XML library for one input format would be a heavier
 * dependency than the format warrants. The .lqnx grammar LINE emits and
 * consumes uses only elements, attributes, character data and comments, with
 * no namespaces to resolve, no DTD, no entity declarations and no processing
 * instructions other than the leading declaration. That is what this parses.
 *
 * The DOM mirrors the two org.w3c.dom operations MATLAB's parseXML uses:
 * getAttribute (absent attribute reads as the empty string) and
 * getElementsByTagName (a DESCENDANT search, not a child search -- this
 * distinction matters: parseXML relies on it to reach `task` elements nested
 * inside `processor`, and separately guards against it with an explicit
 * parent-name test when collecting `activity` under `task-activities`).
 *
 * REFUSES rather than guesses: an unterminated tag, a mismatched close tag or
 * an unquoted attribute value is an InputError. A silently mis-parsed model
 * would produce a solvable network with the wrong topology, which is worse
 * than no answer.
 *
 * The write side mirrors the four org.w3c.dom operations the JMT writers use --
 * createElement, setAttribute, appendChild and createTextNode -- and serializes
 * with `serialize`. TEXT-ONLY ELEMENTS ARE EMITTED INLINE, with no indentation
 * inside the tags: JMT reads `<value>` bodies with Double.parseDouble and a
 * newline plus leading spaces around the number is not what MATLAB's xmlwrite
 * produces either. Mixed content (text alongside child elements) does not occur
 * in either grammar and is not represented.
 */

#include <cctype>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "line/util/error.h"

namespace line {
namespace xml {

struct Element {
    std::string name;
    std::vector<std::pair<std::string, std::string>> attrs;
    std::vector<std::unique_ptr<Element>> children;
    const Element* parent = nullptr;
    /** Character data of a text-only element, as `createTextNode` supplies it. */
    std::string text;

    /** Attribute value, or the empty string when absent (org.w3c.dom semantics). */
    std::string attr(const std::string& key) const {
        for (const auto& kv : attrs)
            if (kv.first == key) return kv.second;
        return std::string();
    }

    bool has_attr(const std::string& key) const {
        for (const auto& kv : attrs)
            if (kv.first == key) return true;
        return false;
    }

    /** Descendant-or-self search excluding self, in document order. */
    std::vector<const Element*> by_tag(const std::string& tag) const {
        std::vector<const Element*> out;
        collect(tag, out);
        return out;
    }

    /** Direct children with the given tag, in document order. */
    std::vector<const Element*> child_tags(const std::string& tag) const {
        std::vector<const Element*> out;
        for (const auto& c : children)
            if (c->name == tag) out.push_back(c.get());
        return out;
    }

    /**
     * `setAttribute`: replace the value in place when the key already exists,
     * otherwise append. Replacing in place is what keeps attribute order stable
     * across a writer that sets `className` twice -- the JSIM writer sets the
     * LINE section name first and overwrites it with the JMT one.
     */
    Element& set_attr(const std::string& key, const std::string& value) {
        for (auto& kv : attrs)
            if (kv.first == key) {
                kv.second = value;
                return *this;
            }
        attrs.emplace_back(key, value);
        return *this;
    }

    /** `createElement` + `appendChild` in one step; the child is owned here. */
    Element& add_child(const std::string& tag) {
        std::unique_ptr<Element> c(new Element());
        c->name = tag;
        c->parent = this;
        children.push_back(std::move(c));
        return *children.back();
    }

    /** `createTextNode` + `appendChild` on an element with no child elements. */
    Element& add_text(const std::string& value) {
        text += value;
        return *this;
    }

    /** The common shape `<tag>value</tag>`. */
    Element& add_text_child(const std::string& tag, const std::string& value) {
        Element& c = add_child(tag);
        c.add_text(value);
        return c;
    }

private:
    void collect(const std::string& tag, std::vector<const Element*>& out) const {
        for (const auto& c : children) {
            if (c->name == tag) out.push_back(c.get());
            c->collect(tag, out);
        }
    }
};

namespace detail {

inline bool is_space(char c) { return c == ' ' || c == '\t' || c == '\n' || c == '\r'; }

inline void skip_space(const std::string& s, std::size_t& i) {
    while (i < s.size() && is_space(s[i])) ++i;
}

/** Expand the five predefined XML entities; anything else is left verbatim. */
inline std::string unescape(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) {
        if (s[i] != '&') {
            out.push_back(s[i]);
            continue;
        }
        const std::size_t semi = s.find(';', i);
        if (semi == std::string::npos) {
            out.push_back(s[i]);
            continue;
        }
        const std::string ent = s.substr(i + 1, semi - i - 1);
        if (ent == "amp") out.push_back('&');
        else if (ent == "lt") out.push_back('<');
        else if (ent == "gt") out.push_back('>');
        else if (ent == "quot") out.push_back('"');
        else if (ent == "apos") out.push_back('\'');
        else {
            out.push_back('&');
            continue;
        }
        i = semi;
    }
    return out;
}

inline std::string read_name(const std::string& s, std::size_t& i) {
    const std::size_t start = i;
    while (i < s.size() && !is_space(s[i]) && s[i] != '>' && s[i] != '/' && s[i] != '=') ++i;
    if (i == start) throw InputError("xml: expected a name");
    return s.substr(start, i - start);
}

}  // namespace detail

/** Parse a whole XML document and return its root element. */
inline std::unique_ptr<Element> parse(const std::string& text) {
    std::size_t i = 0;
    std::unique_ptr<Element> root;
    std::vector<Element*> stack;

    while (i < text.size()) {
        if (text[i] != '<') {
            // Character data. WHITESPACE-ONLY RUNS ARE DROPPED: they are the
            // pretty-printer's indentation, and keeping them would give every
            // container element a text body and turn a re-serialization into
            // mixed content. The .lqnx reader looks at no text at all; the JMT
            // result readers need the bodies of the leaf elements, which is what
            // survives the trim.
            const std::size_t start = i;
            while (i < text.size() && text[i] != '<') ++i;
            const std::string raw = text.substr(start, i - start);
            std::size_t b = 0, e = raw.size();
            while (b < e && detail::is_space(raw[b])) ++b;
            while (e > b && detail::is_space(raw[e - 1])) --e;
            if (e > b && !stack.empty()) stack.back()->text += detail::unescape(raw.substr(b, e - b));
            continue;
        }
        if (text.compare(i, 4, "<!--") == 0) {
            const std::size_t end = text.find("-->", i + 4);
            if (end == std::string::npos) throw InputError("xml: unterminated comment");
            i = end + 3;
            continue;
        }
        if (text.compare(i, 9, "<![CDATA[") == 0) {
            const std::size_t end = text.find("]]>", i + 9);
            if (end == std::string::npos) throw InputError("xml: unterminated CDATA section");
            i = end + 3;
            continue;
        }
        if (text.compare(i, 2, "<?") == 0) {
            const std::size_t end = text.find("?>", i + 2);
            if (end == std::string::npos) throw InputError("xml: unterminated processing instruction");
            i = end + 2;
            continue;
        }
        if (text.compare(i, 2, "<!") == 0) {
            const std::size_t end = text.find('>', i + 2);
            if (end == std::string::npos) throw InputError("xml: unterminated declaration");
            i = end + 1;
            continue;
        }
        if (text.compare(i, 2, "</") == 0) {
            i += 2;
            detail::skip_space(text, i);
            const std::string name = detail::read_name(text, i);
            detail::skip_space(text, i);
            if (i >= text.size() || text[i] != '>') throw InputError("xml: malformed close tag");
            ++i;
            if (stack.empty() || stack.back()->name != name)
                throw InputError("xml: close tag </" + name + "> does not match the open element");
            stack.pop_back();
            continue;
        }

        // open or self-closing element
        ++i;
        detail::skip_space(text, i);
        std::unique_ptr<Element> owned(new Element());
        Element* el = owned.get();
        el->name = detail::read_name(text, i);

        while (true) {
            detail::skip_space(text, i);
            if (i >= text.size()) throw InputError("xml: unterminated element <" + el->name + ">");
            if (text[i] == '>' || (text[i] == '/' && i + 1 < text.size() && text[i + 1] == '>')) break;
            const std::string key = detail::read_name(text, i);
            detail::skip_space(text, i);
            if (i >= text.size() || text[i] != '=')
                throw InputError("xml: attribute '" + key + "' has no value");
            ++i;
            detail::skip_space(text, i);
            if (i >= text.size() || (text[i] != '"' && text[i] != '\''))
                throw InputError("xml: attribute '" + key + "' value is not quoted");
            const char quote = text[i++];
            const std::size_t vstart = i;
            while (i < text.size() && text[i] != quote) ++i;
            if (i >= text.size()) throw InputError("xml: unterminated attribute value");
            el->attrs.emplace_back(key, detail::unescape(text.substr(vstart, i - vstart)));
            ++i;
        }

        const bool self_closing = text[i] == '/';
        i += self_closing ? 2 : 1;

        el->parent = stack.empty() ? nullptr : stack.back();
        Element* raw = el;
        if (stack.empty()) {
            if (root) throw InputError("xml: more than one root element");
            root = std::move(owned);
        } else {
            stack.back()->children.push_back(std::move(owned));
        }
        if (!self_closing) stack.push_back(raw);
    }

    if (!stack.empty()) throw InputError("xml: unterminated element <" + stack.back()->name + ">");
    if (!root) throw InputError("xml: the document has no root element");
    return root;
}

/** Read and parse a file. */
inline std::unique_ptr<Element> parse_file(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw InputError("xml: cannot open " + path);
    std::ostringstream ss;
    ss << f.rdbuf();
    return parse(ss.str());
}

/** A fresh detached element, the `createElement` of the write side. */
inline std::unique_ptr<Element> element(const std::string& tag) {
    std::unique_ptr<Element> e(new Element());
    e->name = tag;
    return e;
}

namespace detail {

/**
 * Escape for character data and for attribute values alike.
 *
 * `>` is escaped although only `]]>` requires it, and `"`/`'` although only the
 * matching quote requires it: one escaper for both positions cannot be applied
 * in the wrong place, and every XML reader expands all five entities.
 */
inline std::string escape(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) {
        switch (s[i]) {
            case '&': out += "&amp;"; break;
            case '<': out += "&lt;"; break;
            case '>': out += "&gt;"; break;
            case '"': out += "&quot;"; break;
            case '\'': out += "&apos;"; break;
            default: out.push_back(s[i]);
        }
    }
    return out;
}

inline void serialize_element(const Element& e, int depth, std::string& out) {
    const std::string pad(static_cast<std::size_t>(depth) * 2, ' ');
    out += pad;
    out += "<";
    out += e.name;
    for (std::size_t i = 0; i < e.attrs.size(); ++i) {
        out += " ";
        out += e.attrs[i].first;
        out += "=\"";
        out += escape(e.attrs[i].second);
        out += "\"";
    }
    if (e.children.empty() && e.text.empty()) {
        out += "/>\n";
        return;
    }
    out += ">";
    if (e.children.empty()) {
        out += escape(e.text);
        out += "</";
        out += e.name;
        out += ">\n";
        return;
    }
    out += "\n";
    for (std::size_t i = 0; i < e.children.size(); ++i)
        serialize_element(*e.children[i], depth + 1, out);
    out += pad;
    out += "</";
    out += e.name;
    out += ">\n";
}

}  // namespace detail

/**
 * Serialize a document: the XML declaration MATLAB's xmlwrite emits, then the
 * root subtree.
 */
inline std::string serialize(const Element& root) {
    std::string out = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    detail::serialize_element(root, 0, out);
    return out;
}

/** Serialize to a file, creating or truncating it. */
inline void write_file(const std::string& path, const Element& root) {
    std::ofstream f(path, std::ios::binary);
    if (!f) throw InputError("xml: cannot write " + path);
    const std::string text = serialize(root);
    f.write(text.data(), static_cast<std::streamsize>(text.size()));
    if (!f) throw InputError("xml: write failed on " + path);
}

}  // namespace xml
}  // namespace line

#endif  // LINE_UTIL_XML_H
