package jline.lang.layerednetworks;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class ActivityPrecedenceTest {

    @Test
    void getPrecedenceId() {
        assertEquals(ActivityPrecedence.getPrecedenceId("pre"),1);
        assertEquals(ActivityPrecedence.getPrecedenceId("pre-AND"),2);
        assertEquals(ActivityPrecedence.getPrecedenceId("pre-OR"),3);
        assertEquals(ActivityPrecedence.getPrecedenceId("post"),11);
        assertEquals(ActivityPrecedence.getPrecedenceId("post-AND"),12);
        assertEquals(ActivityPrecedence.getPrecedenceId("post-OR"),13);
        assertEquals(ActivityPrecedence.getPrecedenceId("post-LOOP"),14);
        assertEquals(ActivityPrecedence.getPrecedenceId("post-CACHE"),15);

    }
}