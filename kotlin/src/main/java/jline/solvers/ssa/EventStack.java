package jline.solvers.ssa;

import jline.lang.nodes.Node;
import jline.lang.nodes.Source;
import jline.solvers.ctmc.EventData;
import jline.solvers.ssa.events.Event;
import jline.solvers.ssa.events.NodeEvent;
import jline.solvers.ssa.state.SSAStateMatrix;
import jline.lang.distributions.CumulativeDistribution;
import jline.util.Pair;
import jline.lang.OutputStrategy;

import java.util.*;

public class EventStack {
    /*
        EventStack -
            Manages the following:
            1. Which event(s) fire(s) next?
            2. Firing the next event(s)
            3. Tracking progress of time
            4. Configuration of tau leaping
     */
    protected List<Event> eventList;

    protected boolean fixedEventOrder;

    protected double curT;

    public EventStack() {
        this.eventList = new ArrayList<Event>();
        //this.eventList = new LinkedList<Event>();
        this.fixedEventOrder = false;
    }

    public void addEvent(Event event) {
        this.eventList.add(event);
    }

     void handleImmediate(SSAStateMatrix networkState, Timeline ssarunner, Random random) {
        /*
            Once an immediate event has been found, try to find any others. Build a cdf, and fire.
         */
        CumulativeDistribution<Event> immediateCumulativeDistribution = new CumulativeDistribution<Event>(random);
        double totalImmediate = 0;

        for (Event event : this.eventList) {
            if (event.getRate(networkState) == Double.POSITIVE_INFINITY) {
                totalImmediate += 1.0;
                immediateCumulativeDistribution.addElement(event, 1.0);
            }
        }

        immediateCumulativeDistribution.normalize(totalImmediate);

        Event e = immediateCumulativeDistribution.sample(random);
        ssarunner.afterEvent(e, networkState);
        if (ssarunner != null) {
            //timeline.record(t, e, stateMatrix);
        }
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

    public double updateState(SSAStateMatrix networkState, Timeline ssarunner, double t, Random random) {
        /*
            This uses the generic Gillespie algorithm to determine and fire the next event
         */
        CumulativeDistribution<Event> eventCumulativeDistribution = new CumulativeDistribution<Event>(random);
        double totalRate = 0;

        boolean foundEvent = false;

        for (Event event : this.eventList) {
            double eventRate = event.getRate(networkState);
            if (eventRate == Double.POSITIVE_INFINITY) {
                this.handleImmediate(networkState, ssarunner, random);
                return t;
            } else if (Double.isNaN(eventRate)) {
                continue;
            }

            foundEvent = true;

            totalRate += eventRate;
            eventCumulativeDistribution.addElement(event, eventRate);
        }

        eventCumulativeDistribution.normalize(totalRate);

        if (!foundEvent) {
            System.out.println("No event found!");
            return t;
        }

        double timeDelta = Math.log(1-random.nextDouble())/(-totalRate);
        t += timeDelta;
        this.curT = t;
        ssarunner.setTime(this.curT);

        Event chosenEvent = eventCumulativeDistribution.sample(random);
        ssarunner.afterEvent(chosenEvent,networkState);

        return t;
    }

    public void updateStateSpace(ArrayList<SSAStateMatrix> stateSpace, Queue<SSAStateMatrix> queue, Set<SSAStateMatrix> stateSet) {
        while (!queue.isEmpty()){

            SSAStateMatrix networkState1 = new SSAStateMatrix(queue.remove());
            ArrayList<Event> eventArrayList = new ArrayList<>();
            for (Event event : this.eventList) {
                double eventRate = event.getRate(networkState1);
                //TODO handleImmediate
//                if (eventRate == Double.POSITIVE_INFINITY) {
//                    this.handleImmediate(stateMatrix1, timeline, t, random);
//                    return t;
//                }
                if (Double.isNaN(eventRate)) {
                    continue;
                }
                eventArrayList.add(event);
            }
            for (Event event : eventArrayList) {
                event.getNextState(networkState1, stateSpace, queue, stateSet);
            }
        }
    }

    public void updateEventSpace(ArrayList<EventData> eventSpace, Queue<SSAStateMatrix> queue, Set<EventData> eventSet) {
        while (!queue.isEmpty()){

            SSAStateMatrix networkState1 = new SSAStateMatrix(queue.remove());

            ArrayList<Event> eventArrayList = new ArrayList<>();

            for (Event event : this.eventList) {
                double eventRate = event.getRate(networkState1);

                if (Double.isNaN(eventRate) || eventRate==0) {
                    continue;
                }


                eventArrayList.add(event);
            }

            for (Event event : eventArrayList) {
               event.getNextEventState(networkState1, eventSpace,event,queue, networkState1, eventSet);
            }
        }
    }

}
