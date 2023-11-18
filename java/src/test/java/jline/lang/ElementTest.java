package jline.lang;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class ElementTest {
    Element element = new Element("element1");
    @Test
    void getName() {
        assertEquals(element.getName(),"element1");
    }

    @Test
    void setName() {
        element.setName("element2");
        assertEquals(element.getName(),"element2");
    }
}