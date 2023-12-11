package jline.lib.ltinversion;

import org.apfloat.Apcomplex;
import java.util.function.UnaryOperator;

/**
 * Function wrapper class, allowing us to add and multiply using the complex high-precision Apcomplex by overriding UnaryOperator
 * Code was adapted from <a href="https://github.com/apache/commons-math/blob/1119bdbdefd0f6decb2211fd00cd3480b7e9dc13/commons-math-legacy/src/main/java/org/apache/commons/math4/legacy/analysis/FunctionUtils.java#L104">...</a>
 */
public class function_wrapper {
    @SafeVarargs
    public static UnaryOperator<Apcomplex> add(UnaryOperator<Apcomplex>... fun) {
        return value -> {
            Apcomplex result = fun[0].apply(value);
            for (int i = 1; i < fun.length; i++) {
                result = result.add(fun[i].apply(value));
            }
            return result;
        };
    }

    @SafeVarargs
    public static UnaryOperator<Apcomplex> multiply(UnaryOperator<Apcomplex>... fun) {
        return value -> {
            Apcomplex result = fun[0].apply(value);
            for (int i = 1; i < fun.length; i++) {
                result = result.multiply(fun[i].apply(value));
            }
            return result;
        };
    }
}
