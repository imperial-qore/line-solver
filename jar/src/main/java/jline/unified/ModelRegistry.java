/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.unified;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.examples.java.advanced.CDFRespTModel;
import jline.examples.java.advanced.CyclicPollingModel;
import jline.examples.java.basic.ForkJoinModel;
import jline.examples.java.advanced.InitStateModel;
import jline.examples.java.advanced.LoadDependentModel;
import jline.examples.java.basic.MixedModel;
import jline.examples.java.basic.PrioModel;
import jline.examples.java.advanced.StateDepRoutingModel;
import jline.examples.java.advanced.StateProbabilitiesModel;
import jline.examples.java.basic.StochPetriNetModel;
import jline.examples.java.basic.ClosedModel;
import jline.examples.java.basic.OpenModel;
import jline.lang.Network;

import java.util.function.Supplier;

/**
 * Registry that maps model names to builder functions for unified cross-language testing.
 */
public final class ModelRegistry {
    public static final ModelRegistry INSTANCE = new ModelRegistry();

    private final Map<String, Supplier<Network>> modelMap = new HashMap<String, Supplier<Network>>();

    private ModelRegistry() {
        registerAllModels();
    }

    private void registerAllModels() {
        // Open networks
        register("oqn_basic", new Supplier<Network>() { public Network get() { return OpenModel.oqn_basic(); } });
        register("oqn_fourqueues", new Supplier<Network>() { public Network get() { return OpenModel.oqn_fourqueues(); } });
        register("oqn_oneline", new Supplier<Network>() { public Network get() { return OpenModel.oqn_oneline(); } });
        register("oqn_cs_routing", new Supplier<Network>() { public Network get() { return OpenModel.oqn_cs_routing(); } });
        register("oqn_vsinks", new Supplier<Network>() { public Network get() { return OpenModel.oqn_vsinks(); } });

        // Closed networks
        register("cqn_repairmen", new Supplier<Network>() { public Network get() { return ClosedModel.cqn_repairmen(); } });
        register("cqn_bcmp_theorem", new Supplier<Network>() { public Network get() { return ClosedModel.cqn_bcmp_theorem_ps(); } });
        register("cqn_multiserver", new Supplier<Network>() { public Network get() { return ClosedModel.cqn_multiserver(); } });
        register("cqn_scheduling_dps", new Supplier<Network>() { public Network get() { return ClosedModel.cqn_scheduling_dps(); } });
        register("cqn_twoclass_erl", new Supplier<Network>() { public Network get() { return ClosedModel.cqn_twoclass_erl(); } });
        register("cqn_twoqueues_multi", new Supplier<Network>() { public Network get() { return ClosedModel.cqn_twoqueues_multi(); } });
        register("cqn_oneline", new Supplier<Network>() { public Network get() { return ClosedModel.cqn_oneline(); } });
        register("cqn_repairmen_multi", new Supplier<Network>() { public Network get() { return ClosedModel.cqn_repairmen_multi(); } });
        register("cqn_mmpp2_service", new Supplier<Network>() { public Network get() { return ClosedModel.cqn_mmpp2_service(); } });
        register("cqn_twoclass_hyperl", new Supplier<Network>() { public Network get() { return ClosedModel.cqn_twoclass_hyperl(); } });
        register("cqn_threeclass_hyperl", new Supplier<Network>() { public Network get() { return ClosedModel.cqn_threeclass_hyperl(); } });

        // Mixed networks
        register("mqn_basic", new Supplier<Network>() { public Network get() { return MixedModel.mqn_basic(); } });
        register("mqn_singleserver_fcfs", new Supplier<Network>() { public Network get() { return MixedModel.mqn_singleserver_fcfs(); } });
        register("mqn_singleserver_ps", new Supplier<Network>() { public Network get() { return MixedModel.mqn_singleserver_ps(); } });
        register("mqn_multiserver_fcfs", new Supplier<Network>() { public Network get() { return MixedModel.mqn_multiserver_fcfs(); } });
        register("mqn_multiserver_ps", new Supplier<Network>() { public Network get() { return MixedModel.mqn_multiserver_ps(); } });

        // Priority models
        register("prio_identical", new Supplier<Network>() { public Network get() { return PrioModel.prio_identical(); } });
        register("prio_hol_open", new Supplier<Network>() { public Network get() { return PrioModel.prio_hol_open(); } });
        register("prio_hol_closed", new Supplier<Network>() { public Network get() { return PrioModel.prio_hol_closed(); } });
        register("prio_psprio", new Supplier<Network>() { public Network get() { return PrioModel.prio_psprio(); } });

        // Fork-Join models
        register("fj_asymm", new Supplier<Network>() { public Network get() { return ForkJoinModel.fj_asymm(); } });
        register("fj_basic_closed", new Supplier<Network>() { public Network get() { return ForkJoinModel.fj_basic_closed(); } });
        register("fj_basic_nesting", new Supplier<Network>() { public Network get() { return ForkJoinModel.fj_basic_nesting(); } });
        register("fj_basic_open", new Supplier<Network>() { public Network get() { return ForkJoinModel.fj_basic_open(); } });
        register("fj_complex_serial", new Supplier<Network>() { public Network get() { return ForkJoinModel.fj_complex_serial(); } });
        register("fj_cs_multi_visits", new Supplier<Network>() { public Network get() { return ForkJoinModel.fj_cs_multi_visits(); } });
        register("fj_cs_postfork", new Supplier<Network>() { public Network get() { return ForkJoinModel.fj_cs_postfork(); } });
        register("fj_cs_prefork", new Supplier<Network>() { public Network get() { return ForkJoinModel.fj_cs_prefork(); } });
        register("fj_deep_nesting", new Supplier<Network>() { public Network get() { return ForkJoinModel.fj_deep_nesting(); } });
        register("fj_delays", new Supplier<Network>() { public Network get() { return ForkJoinModel.fj_delays(); } });
        register("fj_nojoin", new Supplier<Network>() { public Network get() { return ForkJoinModel.fj_nojoin(); } });
        register("fj_route_overlap", new Supplier<Network>() { public Network get() { return ForkJoinModel.fj_route_overlap(); } });
        register("fj_serialfjs_closed", new Supplier<Network>() { public Network get() { return ForkJoinModel.fj_serialfjs_closed(); } });
        register("fj_serialfjs_open", new Supplier<Network>() { public Network get() { return ForkJoinModel.fj_serialfjs_open(); } });
        register("fj_threebranches", new Supplier<Network>() { public Network get() { return ForkJoinModel.fj_threebranches(); } });
        register("fj_twoclasses_forked", new Supplier<Network>() { public Network get() { return ForkJoinModel.fj_twoclasses_forked(); } });

        // Polling models
        register("polling_exhaustive_det", new Supplier<Network>() { public Network get() { return CyclicPollingModel.polling_exhaustive_det(); } });
        register("polling_exhaustive_exp", new Supplier<Network>() { public Network get() { return CyclicPollingModel.polling_exhaustive_exp(); } });
        register("polling_gated", new Supplier<Network>() { public Network get() { return CyclicPollingModel.polling_gated(); } });
        register("polling_klimited", new Supplier<Network>() { public Network get() { return CyclicPollingModel.polling_klimited(); } });
        register("polling_decrementing", new Supplier<Network>() { public Network get() { return CyclicPollingModel.polling_decrementing(); } });

        // Stochastic Petri Net models
        register("spn_basic_closed", new Supplier<Network>() { public Network get() { return StochPetriNetModel.spn_basic_closed(); } });
        register("spn_basic_open", new Supplier<Network>() { public Network get() { return StochPetriNetModel.spn_basic_open(); } });
        register("spn_closed_fourplaces", new Supplier<Network>() { public Network get() { return StochPetriNetModel.spn_closed_fourplaces(); } });
        register("spn_closed_twoplaces", new Supplier<Network>() { public Network get() { return StochPetriNetModel.spn_closed_twoplaces(); } });
        register("spn_fourmodes", new Supplier<Network>() { public Network get() { return StochPetriNetModel.spn_fourmodes(); } });
        register("spn_inhibiting", new Supplier<Network>() { public Network get() { return StochPetriNetModel.spn_inhibiting(); } });
        register("spn_open_sevenplaces", new Supplier<Network>() { public Network get() { return StochPetriNetModel.spn_open_sevenplaces(); } });
        register("spn_twomodes", new Supplier<Network>() { public Network get() { return StochPetriNetModel.spn_twomodes(); } });

        // State-dependent routing
        register("sdroute_open", new Supplier<Network>() { public Network get() { return StateDepRoutingModel.sdroute_open(); } });
        register("sdroute_closed", new Supplier<Network>() { public Network get() { return StateDepRoutingModel.sdroute_closed(); } });
        register("sdroute_twoclasses_closed", new Supplier<Network>() { public Network get() { return StateDepRoutingModel.sdroute_twoclasses_closed(); } });

        // Initial state models
        register("init_state_ps", new Supplier<Network>() { public Network get() { return InitStateModel.init_state_ps(); } });
        register("init_state_fcfs_exp", new Supplier<Network>() { public Network get() { return InitStateModel.init_state_fcfs_exp(); } });
        register("init_state_fcfs_nonexp", new Supplier<Network>() { public Network get() { return InitStateModel.init_state_fcfs_nonexp(); } });

        // Load-dependent models
        register("ld_class_dependence", new Supplier<Network>() { public Network get() { return LoadDependentModel.ld_class_dependence(); } });
        register("ld_multiserver_fcfs", new Supplier<Network>() { public Network get() { return LoadDependentModel.ld_multiserver_fcfs(); } });
        register("ld_multiserver_ps", new Supplier<Network>() { public Network get() { return LoadDependentModel.ld_multiserver_ps(); } });
        register("ld_multiserver_ps_twoclasses", new Supplier<Network>() { public Network get() { return LoadDependentModel.ld_multiserver_ps_twoclasses(); } });

        // State probabilities models
        register("statepr_aggr", new Supplier<Network>() { public Network get() { return StateProbabilitiesModel.statepr_aggr(); } });
        register("statepr_aggr_large", new Supplier<Network>() { public Network get() { return StateProbabilitiesModel.statepr_aggr_large(); } });
        register("statepr_sys_aggr", new Supplier<Network>() { public Network get() { return StateProbabilitiesModel.statepr_sys_aggr(); } });
        register("statepr_sys_aggr_large", new Supplier<Network>() { public Network get() { return StateProbabilitiesModel.statepr_sys_aggr_large(); } });
        register("statepr_allprobs_ps", new Supplier<Network>() { public Network get() { return StateProbabilitiesModel.statepr_allprobs_ps(); } });
        register("statepr_allprobs_fcfs", new Supplier<Network>() { public Network get() { return StateProbabilitiesModel.statepr_allprobs_fcfs(); } });

        // CDF Response Time models
        register("cdf_respt_closed", new Supplier<Network>() { public Network get() { return CDFRespTModel.cdf_respt_closed(); } });
        register("cdf_respt_closed_threeclasses", new Supplier<Network>() { public Network get() { return CDFRespTModel.cdf_respt_closed_threeclasses(); } });
        register("cdf_respt_open_twoclasses", new Supplier<Network>() { public Network get() { return CDFRespTModel.cdf_respt_open_twoclasses(); } });
        register("cdf_respt_distrib", new Supplier<Network>() { public Network get() { return CDFRespTModel.cdf_respt_distrib(); } });
        register("cdf_respt_populations", new Supplier<Network>() { public Network get() { return CDFRespTModel.cdf_respt_populations(); } });
    }

    private void register(String name, Supplier<Network> builder) {
        modelMap.put(name, builder);
    }

    /**
     * Get the singleton instance (for Java interop)
     */
    public static ModelRegistry getInstance() {
        return INSTANCE;
    }

    /**
     * Get a model by name.
     */
    public static Network getModel(String name) {
        Supplier<Network> builder = INSTANCE.modelMap.get(name);
        if (builder == null) {
            StringBuilder sb = new StringBuilder();
            boolean first = true;
            for (String key : INSTANCE.modelMap.keySet()) {
                if (!first) sb.append(", ");
                sb.append(key);
                first = false;
            }
            throw new IllegalArgumentException(
                    "Model '" + name + "' is not registered. Available models: " + sb.toString());
        }
        return builder.get();
    }

    /**
     * Get list of all registered model names
     */
    public static List<String> getAvailableModels() {
        return new ArrayList<String>(INSTANCE.modelMap.keySet());
    }

    /**
     * Check if a model is registered
     */
    public static boolean hasModel(String name) {
        return INSTANCE.modelMap.containsKey(name);
    }
}
