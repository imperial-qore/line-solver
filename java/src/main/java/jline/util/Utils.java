package jline.util;

import org.apache.commons.math3.special.Gamma;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

/**
 * Miscellaneous utilities
 */
public class Utils {

    // TODO: polymorphic version that also returns nanMean
    // Return mean absolute percentage error of approx with respect to exact
    public static double mape(Matrix approx, Matrix exact) {

      int numRows = approx.getNumRows();
      double totalAbsolutePercentageError = 0;
      int numExactGreaterThanZero = 0;
      for (int row = 0; row < numRows; row++) {
        if (exact.get(row, 0) > 0) {
          totalAbsolutePercentageError += Math.abs(1 - (approx.get(row, 0) / exact.get(row, 0)));
          numExactGreaterThanZero++;
        }
      }
      return totalAbsolutePercentageError / numExactGreaterThanZero;
    }

    public static <T extends Object> List<T> unique(List<T> list) {
        return new ArrayList<T>(new HashSet<>(list));
    }

    public static int findString(Map<Integer, String> map, String target) {
        int res = 0;
        for(Map.Entry e : map.entrySet()){
            if (e.getValue().equals(target)) res = (int) e.getKey();
        }
        return res;
    }
}
