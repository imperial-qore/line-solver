/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `python/examples/advanced/switchoverTimes/`: the time a server needs to
 * switch from serving one class to another.
 *
 * REFUSED, AND THE REFUSAL IS THE ANSWER. The reference declares a switchover
 * per ORDERED PAIR of classes at an FCFS station
 * (`queue.setSwitchover(class1, class2, dist)`); this port's
 * `Network::set_switchover` is per CLASS and is accepted only at a POLLING
 * station, because that is the only place `Station::switchover` is read -- the
 * per-buffer walk of `solver_mva_polling_analyzer`. Dropping the pair structure
 * onto the per-class vector would build a DIFFERENT model under this name, and
 * declaring nothing would make the switchover free, so the example declares its
 * model unrepresentable instead.
 *
 * The reference solves it with JMT alone, which this port does not carry
 * either, so nothing is lost that the port could otherwise have computed.
 */

#include "examples_common.h"

namespace line {
namespace examples {

/** `switchover_basic.py`: pairwise switchover times at an FCFS queue. */
void switchover_basic() {
    // TODO(cpp): avg_table_jmt = JMT(model, seed=23000, keep=True).get_avg_table(); print(avg_table_jmt)
    na("JMT",
       "the MODEL, not the engine: `jmt_avg` drives the same JMT.jar the reference drives, but "
       "the model below cannot be built in this port at all, so there is nothing to hand it");
    note("N/A: Queue.setSwitchover(fromClass, toClass, dist) has no counterpart in this port. "
         "Network::set_switchover declares one distribution per class and only at a POLLING "
         "station, so the pairwise switchover of this model cannot be expressed.");
}

LINE_EXAMPLE("advanced/switchoverTimes", switchover_basic);

}  // namespace examples
}  // namespace line
