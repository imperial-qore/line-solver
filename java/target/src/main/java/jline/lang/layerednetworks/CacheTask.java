package jline.lang.layerednetworks;

/**
 * A task that offers caching services
 */
public class CacheTask extends Task {
    protected int items;

    protected int itemLevelCap;

    protected int replacementPolicy;
    public CacheTask(LayeredNetwork model, String name){
        super(model,name);

    }
}
