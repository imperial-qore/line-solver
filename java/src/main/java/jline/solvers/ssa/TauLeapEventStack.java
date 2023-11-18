package jline.solvers.ssa;

import jline.lang.OutputStrategy;
import jline.lang.distributions.Binomial;
import jline.lang.distributions.Poisson;
import jline.lang.nodes.Node;
import jline.lang.nodes.Source;
import jline.solvers.ssa.events.Event;
import jline.solvers.ssa.events.NodeEvent;
import jline.solvers.ssa.metrics.IllegalTauLeapPercentageMetric;
import jline.solvers.ssa.state.SSAStateMatrix;
import jline.solvers.ssa.strategies.TauLeapingOrderStrategy;
import jline.solvers.ssa.strategies.TauLeapingStateStrategy;
import jline.solvers.ssa.strategies.TauLeapingType;
import jline.solvers.ssa.strategies.TauLeapingVarType;
import jline.util.Pair;

import java.util.*;

public class TauLeapEventStack extends EventStack {
    protected TauLeapingType tauLeapingType;

    protected int tauFailures;

    protected IllegalTauLeapPercentageMetric illegalTauLeapPercentage;

    public TauLeapEventStack() {
        this.tauLeapingType = null;

        this.tauFailures = 0;
        this.illegalTauLeapPercentage = new IllegalTauLeapPercentageMetric();
    }

    public double getIllegalTauLeapPercentage() {
        return this.illegalTauLeapPercentage.getMetric();
    }

    public void configureTauLeap(TauLeapingType tauLeapingType) {
        this.tauLeapingType = tauLeapingType;
    }

    @SuppressWarnings("unchecked")
    private void orderEventList(Random random) {
        /*
                For DirectedGraph and DirectedCycle ordering methods, re-order the event list according to a topological
                    sort.
         */
        List<Event> outEvents = new ArrayList<Event>();
        Queue<Pair<Node, Event>> candidateEvents = new LinkedList<Pair<Node, Event>>();
        List<Pair<Node, Event>> nodeEvents = new LinkedList<Pair<Node, Event>>();

        // Add any events that don't correspond to a node, first
        for (Event event : this.eventList) {
            if (event instanceof NodeEvent) {
                Node node = ((NodeEvent)event).getNode();
                nodeEvents.add(new Pair(node, event));
            } else {
                outEvents.add(event);
            }
        }

        Collections.shuffle(outEvents, random);

        Iterator<Pair<Node,Event>> eventIterator = nodeEvents.iterator();
        // Iterate through each event, add any that correspond to a Source or a reference station to canidateEvents
        while (eventIterator.hasNext()) {
            Pair<Node,Event> iterPair = eventIterator.next();
            Node eventNode = iterPair.getLeft();

            if ((eventNode instanceof Source) || (eventNode.isRefstat())) {
                candidateEvents.add(iterPair);
                eventIterator.remove();
            }
        }

        // Dequeue any events in candidateEvent, add any new events that are further in the topology, and append it
        //   to outEvents
        while(!candidateEvents.isEmpty()) {
            Pair<Node,Event> iterPair = candidateEvents.remove();
            List<OutputStrategy> outputStrategies = iterPair.getLeft().getOutputStrategies();
            List<Node> outNodes = new ArrayList<Node>();

            for (OutputStrategy outputStrategy : outputStrategies) {
                Node dest = outputStrategy.getDestination();
                if (!outNodes.contains(dest)) {
                    outNodes.add(dest);
                }
            }

            Event event = iterPair.getRight();
            outEvents.add(event);
            eventIterator = nodeEvents.iterator();

            while(eventIterator.hasNext()) {
                Pair<Node, Event> iterPair2 = eventIterator.next();
                Node eventNode = iterPair2.getLeft();
                event = iterPair2.getRight();

                if (outNodes.contains(eventNode)) {
                    candidateEvents.add(iterPair2);
                    eventIterator.remove();
                }
            }
        }

        // Add any remaining events
        for (Pair<Node, Event> iterPair : nodeEvents) {
            outEvents.add(iterPair.getRight());
        }

        this.eventList = outEvents;
    }

    public double tauLeapUpdateMultistep(SSAStateMatrix networkState, Timeline ssarunner, double t, Random random) {
        /*
                This is part of an optimization, but I've seperated TauLeapUpdate into two functions, one for TwoTimes
                and one for all other ordering strategies. 99% of this is a simple copy+paste from an earlier iteration.

                This is because TwoTimes needs to keep track of unapplied iterations for all events.
         */
        List<Pair<Event, Integer>> eventPairs = new ArrayList<Pair<Event, Integer>>();

        double tau = this.tauLeapingType.getTau();
        boolean foundEvent = false;

        TauLeapingVarType varType = this.tauLeapingType.getVarType();
        TauLeapingOrderStrategy orderStrategy = this.tauLeapingType.getOrderStrategy();
        TauLeapingStateStrategy stateStrategy = this.tauLeapingType.getStateStrategy();

        // Handle any ordering issues
        if (orderStrategy == TauLeapingOrderStrategy.RandomEventFixed) {
            if (!this.fixedEventOrder) {
                Collections.shuffle(this.eventList, random);
                this.fixedEventOrder = true;
            }
        } else if (orderStrategy == TauLeapingOrderStrategy.DirectedGraph) {
            if (!this.fixedEventOrder) {
                this.orderEventList(random);
                this.fixedEventOrder = true;
            }
        } else if (orderStrategy == TauLeapingOrderStrategy.RandomEvent) {
            Collections.shuffle(this.eventList, random);
        } else if (orderStrategy == TauLeapingOrderStrategy.DirectedCycle) {
            if (!this.fixedEventOrder) {
                this.orderEventList(random);
                this.fixedEventOrder = true;
            }

            // cycle through directed graph
            Event first = this.eventList.get(0);
            this.eventList.remove(0);
            this.eventList.add(first);
        }

        // Build a list of all events, as well as firing counts
        for (Event event : this.eventList) {
            double eventRate = event.getRate(networkState);
            if (eventRate == Double.POSITIVE_INFINITY) {
                this.handleImmediate(networkState, ssarunner, random);
                return t;
            } else if (Double.isNaN(eventRate)) {
                continue;
            }

            foundEvent = true;

            int eCount = 0;

            final int MAX_BINOMIAL_REPS = 1000;

            TauLeapingVarType effectiveVarType = varType;

            if (effectiveVarType == TauLeapingVarType.Binomial) {
                int maxReps = event.getMaxRepetitions(networkState);
                if (maxReps > MAX_BINOMIAL_REPS) {
                    effectiveVarType = TauLeapingVarType.Poisson;
                }
            }

            if (effectiveVarType == TauLeapingVarType.Poisson) {
                Poisson pDist = new Poisson(eventRate*tau);
                eCount = pDist.getRealization(random);
            } else if (effectiveVarType == TauLeapingVarType.Binomial) {
                int maxJump = event.getMaxRepetitions(networkState);
                double prob = (eventRate*tau)/((double)maxJump);

                prob = Math.min(prob, 1.0);
                if (prob != 0) {
                    Binomial bDist = new Binomial(prob, maxJump);
                    eCount = bDist.getRealization(random);
                }
            }

            if (eCount == 0) {
                continue;
            }

            eventPairs.add(new Pair<Event, Integer>(event, eCount));
        }

        if (!foundEvent) {
            System.out.println("No event found!");
            return t;
        }

        t += tau;
        this.curT = t;
        ssarunner.setNextTime(this.curT);


        // We know that this will only run twice, but it might be
        //      a bit better to make this a variable for future development
        int nTimes = 1;
        if (stateStrategy == TauLeapingStateStrategy.TwoTimes) {
            nTimes = 2;
        }

        for (int timeIter = 0; timeIter < nTimes; timeIter++) {
            Iterator<Pair<Event, Integer>> eventPairIterator = eventPairs.listIterator();
            List<Pair<Event, Integer>> nextList = new ArrayList<Pair<Event, Integer>>();
            while (eventPairIterator.hasNext()) {
                Pair <Event,Integer> eventPair = eventPairIterator.next();
                Event event = eventPair.getLeft();
                int eCount = eventPair.getRight();


                int rem = ssarunner.afterEventN(eCount, event, networkState);

                if ((stateStrategy == TauLeapingStateStrategy.TwoTimes) && (timeIter==0)) {
                    // build another list here - it's a bit faster to build two lists rather than keeping track of
                    //   the end-point
                    nextList.add(new Pair<Event, Integer>(event, rem));
                }

                if ((eCount-rem) == 0) {
                    continue;
                }

                if (rem == 0) {
                    this.illegalTauLeapPercentage.addSuccessful();
                } else {
                    this.illegalTauLeapPercentage.addIllegal();
                }
            }

            // reset the list..
            eventPairs = nextList;
        }

        ssarunner.recordCache();

        return t;
    }

    public double tauLeapUpdateOnestep(SSAStateMatrix networkState, Timeline ssarunner, double t, Random random) {
        double tau = this.tauLeapingType.getTau();
        boolean foundEvent = false;
        final int MAX_BINOMIAL_REPS = 1000;

        TauLeapingVarType varType = this.tauLeapingType.getVarType();
        TauLeapingOrderStrategy orderStrategy = this.tauLeapingType.getOrderStrategy();
        TauLeapingStateStrategy stateStrategy = this.tauLeapingType.getStateStrategy();

        // apply ordering changes
        if (orderStrategy == TauLeapingOrderStrategy.RandomEventFixed) {
            if (!this.fixedEventOrder) {
                Collections.shuffle(this.eventList, random);
                this.fixedEventOrder = true;
            }
        } else if ((orderStrategy == TauLeapingOrderStrategy.DirectedGraph) ||
                (orderStrategy == TauLeapingOrderStrategy.DirectedCycle)) {
            if (!this.fixedEventOrder) {
                this.orderEventList(random);
                this.fixedEventOrder = true;
            }
        } else if (orderStrategy == TauLeapingOrderStrategy.RandomEvent) {
            Collections.shuffle(this.eventList, random);
        }

        if (orderStrategy == TauLeapingOrderStrategy.DirectedCycle) {
            // cycle through directed graph;
            Event first = this.eventList.get(0);
            this.eventList.remove(0);
            this.eventList.add(first);
        }

        /*if (stateStrategy == TauLeapingStateStrategy.CycleCutoff) {
            stateMatrix.allowIllegalStates();
        } else*/ if ((stateStrategy == TauLeapingStateStrategy.TimeWarp) || (stateStrategy == TauLeapingStateStrategy.TauTimeWarp)) {
            networkState.cacheState();
        }

        t += tau;
        this.curT = t;
        ssarunner.setNextTime(this.curT);

        // loop through each event, calculate number of reps, and fire accordingly
        for (Event event : this.eventList) {
            double eventRate = event.getRate(networkState);
            if (eventRate == Double.POSITIVE_INFINITY) {
                this.handleImmediate(networkState, ssarunner, random);
                return t;
            } else if (Double.isNaN(eventRate)) {
                continue;
            }

            int eCount = 0;

            if (varType == TauLeapingVarType.Poisson) {
                Poisson pDist = new Poisson(eventRate*tau);
                eCount = pDist.getRealization(random);
            } else if (varType == TauLeapingVarType.Binomial) {
                int maxJump = event.getMaxRepetitions(networkState);
                if (maxJump > MAX_BINOMIAL_REPS) {
                    Poisson pDist = new Poisson(eventRate*tau);
                    eCount = pDist.getRealization(random);
                } else {
                    double prob = (eventRate * tau) / ((double) maxJump);

                    prob = Math.min(prob, 1.0);
                    if (prob != 0) {
                        Binomial bDist = new Binomial(prob, maxJump);
                        eCount = bDist.getRealization(random);
                    }
                }
            }

            if (eCount == 0) {
                continue;
            }

            // Find the number of unapplied reps. Note the max - this is necessary because on some occasions
            //    you'll receive a negative number when CycleCutoff is applied
            int rem = Math.max(ssarunner.afterEventN(eCount, event, networkState),0);

            if (stateStrategy == TauLeapingStateStrategy.TimeWarp) {
                if (rem != 0) {
                    if (this.tauFailures == 2) {
                        this.illegalTauLeapPercentage.addIllegal();
                        this.tauFailures = 0;
                    } else {
                        this.tauFailures++;
                        networkState.revertToCache();
                        ssarunner.clearCache();
                        return t - tau;
                    }
                } else {
                    this.illegalTauLeapPercentage.addSuccessful();
                    this.tauFailures = 0;
                }
            } else if (stateStrategy == TauLeapingStateStrategy.TauTimeWarp) {
                if (rem != 0) {
                    if (this.tauFailures == 2) {
                        this.illegalTauLeapPercentage.addIllegal();
                        this.tauFailures = 0;
                        tauLeapingType.resetTau();
                    } else {
                        networkState.revertToCache();
                        ssarunner.clearCache();
                        tauLeapingType.setTau(Math.max(tau / 2, 0.0001));
                        this.tauFailures++;
                        return t - tau;
                    }
                } else {
                    this.illegalTauLeapPercentage.addSuccessful();
                    this.tauFailures = 0;
                    tauLeapingType.resetTau();
                }
            } else if ((eCount-rem) == 0) {
                continue;
            } else if (rem != 0) {
                this.illegalTauLeapPercentage.addIllegal();
            } else {
                this.illegalTauLeapPercentage.addSuccessful();
            }

        }

        /*if (stateStrategy == TauLeapingStateStrategy.CycleCutoff) {
            stateMatrix.forbidIllegalStates();
        }*/

        ssarunner.recordCache();

        return t;
    }


    public double tauLeapUpdate(SSAStateMatrix networkState, Timeline Timeline, double t, Random random) {
        // Wrapper function for a cleaner API
        if (tauLeapingType.getStateStrategy() == TauLeapingStateStrategy.TwoTimes) {
            return tauLeapUpdateMultistep(networkState, Timeline, t, random);
        } else {
            return tauLeapUpdateOnestep(networkState, Timeline, t, random);
        }
    }

}
