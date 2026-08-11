/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lib.lti;

import java.util.function.UnaryOperator;

import org.apfloat.Apcomplex;

/**
 * Function wrapper class, allowing us to add and multiply using the complex high-precision Apcomplex by overriding UnaryOperator.
 * Code was adapted from
 * https://github.com/apache/commons-math/blob/1119bdbdefd0f6decb2211fd00cd3480b7e9dc13/commons-math-legacy/src/main/java/org/apache/commons/math4/legacy/analysis/FunctionUtils.java#L104
 */
public final class function_wrapper {
    private function_wrapper() {}

    @SafeVarargs
    public static UnaryOperator<Apcomplex> add(final UnaryOperator<Apcomplex>... fun) {
        return new UnaryOperator<Apcomplex>() {
            @Override
            public Apcomplex apply(Apcomplex value) {
                Apcomplex result = fun[0].apply(value);
                for (int i = 1; i < fun.length; i++) {
                    result = result.add(fun[i].apply(value));
                }
                return result;
            }
        };
    }

    @SafeVarargs
    public static UnaryOperator<Apcomplex> multiply(final UnaryOperator<Apcomplex>... fun) {
        return new UnaryOperator<Apcomplex>() {
            @Override
            public Apcomplex apply(Apcomplex value) {
                Apcomplex result = fun[0].apply(value);
                for (int i = 1; i < fun.length; i++) {
                    result = result.multiply(fun[i].apply(value));
                }
                return result;
            }
        };
    }
}
