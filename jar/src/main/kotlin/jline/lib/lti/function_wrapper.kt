/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lib.lti

import org.apfloat.Apcomplex
import java.util.function.UnaryOperator

/**
 * Function wrapper class, allowing us to add and multiply using the complex high-precision Apcomplex by overriding UnaryOperator
 * Code was adapted from [...](https://github.com/apache/commons-math/blob/1119bdbdefd0f6decb2211fd00cd3480b7e9dc13/commons-math-legacy/src/main/java/org/apache/commons/math4/legacy/analysis/FunctionUtils.java#L104)
 */
object function_wrapper {
    @SafeVarargs
    fun add(vararg `fun`: UnaryOperator<Apcomplex>): UnaryOperator<Apcomplex> {
        return UnaryOperator { value: Apcomplex ->
            var result = `fun`[0].apply(value)
            for (i in 1..<`fun`.size) {
                result = result.add(`fun`[i].apply(value))
            }
            result
        }
    }

    @SafeVarargs
    fun multiply(vararg `fun`: UnaryOperator<Apcomplex>): UnaryOperator<Apcomplex> {
        return UnaryOperator { value: Apcomplex ->
            var result = `fun`[0].apply(value)
            for (i in 1..<`fun`.size) {
                result = result.multiply(`fun`[i].apply(value))
            }
            result
        }
    }
}
