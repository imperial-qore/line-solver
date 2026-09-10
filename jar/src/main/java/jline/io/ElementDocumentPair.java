/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.io;

import org.w3c.dom.Document;
import org.w3c.dom.Element;

public class ElementDocumentPair {
    public Element simElem;
    public Document simDoc;

    public ElementDocumentPair(Element simElem, Document simDoc) {
        this.simElem = simElem;
        this.simDoc = simDoc;
    }

    public Object getSimElem() {
        return simElem;
    }

    public Object getSimDoc() {
        return simDoc;
    }
}
