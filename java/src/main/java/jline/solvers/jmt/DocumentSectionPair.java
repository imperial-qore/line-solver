package jline.solvers.jmt;

import org.w3c.dom.Document;
import org.w3c.dom.Element;

public class DocumentSectionPair {
    public Document simDoc;
    public Element section;

    public DocumentSectionPair(Document simDoc, Element section){
        this.simDoc = simDoc;
        this.section = section;
    }

    public Document getSimDoc() {
        return simDoc;
    }

    public Element getSection() {
        return section;
    }
}
