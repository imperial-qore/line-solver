// Copyright (c) 2012-2022, Imperial College London
// All rights reserved.

package jline.lang;

import java.util.List;

/**
 * A model defined by a collection of sub-models
 */
public class Ensemble extends Model {

  protected List<Network> ensemble;

  public Ensemble(List<Network> models) {
    super("Ensemble");
    this.ensemble = models;
  }

  public Ensemble(String name) {
    super(name);
  }

  public void setEnsemble(List<Network> ensemble) {
    this.ensemble = ensemble;
  }

  public List<Network> getEnsemble() {
    return this.ensemble;
  }

  public Network getModel(int modelIdx) {
    return this.ensemble.get(modelIdx);
  }

  // NOTE: the following LINE methods have not been migrated to JLINE
  // a) copyElement - overrides an unimplemented method in "Copyable", and is unused
}
