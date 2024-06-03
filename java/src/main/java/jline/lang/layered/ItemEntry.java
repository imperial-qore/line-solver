package jline.lang.layered;


import jline.lang.distributions.Distribution;
import org.apache.commons.lang3.SerializationUtils;

/**
 * A caching service that gives access to items
 */
public class ItemEntry extends Entry {
    protected int cardinality;
    protected Distribution popularity;

    public ItemEntry(String name, LayeredNetwork model, int cardinality, Distribution distribution) {
        super(model, name);
        this.cardinality = cardinality;
        if(distribution.isDiscrete()){
            this.popularity = SerializationUtils.clone(distribution);
        }
    }

    public int getCardinality(){
        return cardinality;
    }

    public Distribution getPopularity(){
        return popularity;
    }

    @Override
    public ItemEntry on(Task newParent){
        newParent.addEntry(this);
        this.parent = newParent;
        return this;
    }
}
