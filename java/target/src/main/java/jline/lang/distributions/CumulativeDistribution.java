package jline.lang.distributions;

import jline.util.Pair;

import java.io.Serializable;
import java.util.*;

public class CumulativeDistribution<T> implements Serializable, Cloneable {
    protected Collection<Pair<Double, T>> pdf; // T is some random event
    protected Random random;
    protected int numElements;

    public CumulativeDistribution(Random random) {
        pdf = new ArrayList<Pair<Double, T>>();
        this.random = random;
        this.numElements = 0;
    }

    // Copy constructor
    public CumulativeDistribution(CumulativeDistribution<T> other) {
        this.pdf = new ArrayList<Pair<Double, T>>();
        for (Pair<Double, T> pair : other.pdf) {
            this.pdf.add(new Pair<Double, T>(pair));
        }
        this.random = other.random;  // This could also be cloned from other if required
        this.numElements = other.numElements;
    }

    public void addElement(T elem, double prob) {
        Pair<Double, T> elementPair = new Pair<Double, T>(prob, elem);
        this.pdf.add(elementPair);
        this.numElements++;
    }

    public T sample(Random random) {
        double serialProb = random.nextDouble();
        double cumProb = 0;
        Iterator<Pair<Double, T>> pdfIter = this.pdf.iterator();
        while (pdfIter.hasNext()) {
            Pair<Double, T> tPair = pdfIter.next();
            cumProb += (tPair.getLeft());
            if (cumProb >= serialProb) {
                return tPair.getRight();
            }
        }
        return null;
    }

    public ArrayList<T> getPossibleEvents(){
        ArrayList<T> arrayList = new ArrayList<>();
        Iterator<Pair<Double, T>> pdfIter = this.pdf.iterator();
        while (pdfIter.hasNext()) {
            Pair<Double, T> tPair = pdfIter.next();
            arrayList.add(tPair.getRight());
        }
        return arrayList;
    }

    public ArrayList<Pair<Double, T>> getPossibleEventProbability(){
        ArrayList<Pair<Double, T>> arrayList = new ArrayList<>();
        Iterator<Pair<Double, T>> pdfIter = this.pdf.iterator();
        while (pdfIter.hasNext()) {
            Pair<Double, T> tPair = pdfIter.next();
            arrayList.add(tPair);
        }
        return arrayList;
    }


    public void normalize(double factor) {
        Iterator<Pair<Double, T>> pdfIter = this.pdf.iterator();
        while (pdfIter.hasNext()) {
            Pair<Double, T> tPair = pdfIter.next();
            double p0 = tPair.getLeft();
            p0 /= factor;
            tPair.setLeft(p0);
        }
    }

    @Override
    public CumulativeDistribution<T> clone() {
        return new CumulativeDistribution<>(this);
    }
}
