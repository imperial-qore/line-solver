package jline.lang;

import jline.lang.nodes.Cache;

/**
 * A set of cacheable items
 */
public class ItemSet extends NetworkElement{
    private final int nitems;
    private int index;
    private final Cache reference;
    private final boolean replicable;

    public ItemSet(Network model, String name, int nitems, Cache reference) {
        super(name);
        this.nitems = nitems;
        this.index = 0;
        this.replicable = false;
        if(reference == null){
            throw new RuntimeException("ItemSet must be pinned to a Cache");
        }
        this.reference = reference;
        model.addItemSet(this);
    }

    public void setIndex(int index) {
        this.index = index;
    }

    public boolean hasReplicableItems(){
        return this.replicable;
    }

    public int getNumberOfItems(){
        return this.nitems;
    }

    public int getIndex() {
        return index;
    }
}
