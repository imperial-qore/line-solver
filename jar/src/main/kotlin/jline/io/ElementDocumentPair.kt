/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.io

import org.w3c.dom.Document
import org.w3c.dom.Element

class ElementDocumentPair(@JvmField var simElem: Element?, @JvmField var simDoc: Document?) {
    fun getSimElem(): Any? {
        return simElem
    }

    fun getSimDoc(): Any? {
        return simDoc
    }
}
