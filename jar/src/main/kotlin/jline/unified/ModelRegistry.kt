/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.unified

import jline.lang.Network
import jline.examples.java.basic.*
import jline.examples.java.advanced.*

/**
 * Registry that maps model names to builder functions for unified cross-language testing.
 *
 * This singleton class provides a registry that maps model names (as used in the unified
 * JSON test definitions) to Java/Kotlin functions that create the corresponding Network objects.
 * It delegates to the existing model implementations in jline.examples.java.basic.* and
 * jline.examples.java.advanced.* packages.
 *
 * Usage:
 * ```kotlin
 * val registry = ModelRegistry.getInstance()
 * val model = registry.getModel("oqn_basic")
 * ```
 */
object ModelRegistry {
    private val modelMap: MutableMap<String, () -> Network> = mutableMapOf()

    init {
        registerAllModels()
    }

    private fun registerAllModels() {
        // =======================================================================
        // Open networks (from OpenModel.java)
        // =======================================================================
        register("oqn_basic") { OpenModel.oqn_basic() }
        register("oqn_fourqueues") { OpenModel.oqn_fourqueues() }
        register("oqn_oneline") { OpenModel.oqn_oneline() }
        register("oqn_cs_routing") { OpenModel.oqn_cs_routing() }
        register("oqn_vsinks") { OpenModel.oqn_vsinks() }

        // =======================================================================
        // Closed networks (from ClosedModel.java)
        // =======================================================================
        register("cqn_repairmen") { ClosedModel.cqn_repairmen() }
        register("cqn_bcmp_theorem") { ClosedModel.cqn_bcmp_theorem_ps() }
        register("cqn_multiserver") { ClosedModel.cqn_multiserver() }
        register("cqn_scheduling_dps") { ClosedModel.cqn_scheduling_dps() }
        register("cqn_twoclass_erl") { ClosedModel.cqn_twoclass_erl() }
        register("cqn_twoqueues") { ClosedModel.cqn_twoqueues() }
        register("cqn_twoqueues_multi") { ClosedModel.cqn_twoqueues_multi() }
        register("cqn_oneline") { ClosedModel.cqn_oneline() }
        register("cqn_repairmen_multi") { ClosedModel.cqn_repairmen_multi() }
        register("cqn_mmpp2_service") { ClosedModel.cqn_mmpp2_service() }
        register("cqn_twoclass_hyperl") { ClosedModel.cqn_twoclass_hyperl() }
        register("cqn_threeclass_hyperl") { ClosedModel.cqn_threeclass_hyperl() }

        // =======================================================================
        // Mixed networks (from MixedModel.java)
        // =======================================================================
        register("mqn_basic") { MixedModel.mqn_basic() }
        register("mqn_singleserver_fcfs") { MixedModel.mqn_singleserver_fcfs() }
        register("mqn_singleserver_ps") { MixedModel.mqn_singleserver_ps() }
        register("mqn_multiserver_fcfs") { MixedModel.mqn_multiserver_fcfs() }
        register("mqn_multiserver_ps") { MixedModel.mqn_multiserver_ps() }

        // =======================================================================
        // Priority models (from PrioModel.java)
        // =======================================================================
        register("prio_identical") { PrioModel.prio_identical() }
        register("prio_hol_open") { PrioModel.prio_hol_open() }
        register("prio_hol_closed") { PrioModel.prio_hol_closed() }
        register("prio_psprio") { PrioModel.prio_psprio() }

        // =======================================================================
        // Fork-Join models (from ForkJoinModel.java)
        // =======================================================================
        register("fj_asymm") { ForkJoinModel.fj_asymm() }
        register("fj_basic_closed") { ForkJoinModel.fj_basic_closed() }
        register("fj_basic_nesting") { ForkJoinModel.fj_basic_nesting() }
        register("fj_basic_open") { ForkJoinModel.fj_basic_open() }
        register("fj_complex_serial") { ForkJoinModel.fj_complex_serial() }
        register("fj_cs_multi_visits") { ForkJoinModel.fj_cs_multi_visits() }
        register("fj_cs_postfork") { ForkJoinModel.fj_cs_postfork() }
        register("fj_cs_prefork") { ForkJoinModel.fj_cs_prefork() }
        register("fj_deep_nesting") { ForkJoinModel.fj_deep_nesting() }
        register("fj_delays") { ForkJoinModel.fj_delays() }
        register("fj_nojoin") { ForkJoinModel.fj_nojoin() }
        register("fj_route_overlap") { ForkJoinModel.fj_route_overlap() }
        register("fj_serialfjs_closed") { ForkJoinModel.fj_serialfjs_closed() }
        register("fj_serialfjs_open") { ForkJoinModel.fj_serialfjs_open() }
        register("fj_threebranches") { ForkJoinModel.fj_threebranches() }
        register("fj_twoclasses_forked") { ForkJoinModel.fj_twoclasses_forked() }

        // =======================================================================
        // Polling models (from CyclicPollingModel.java)
        // =======================================================================
        register("polling_exhaustive_det") { CyclicPollingModel.polling_exhaustive_det() }
        register("polling_exhaustive_exp") { CyclicPollingModel.polling_exhaustive_exp() }
        register("polling_gated") { CyclicPollingModel.polling_gated() }
        register("polling_klimited") { CyclicPollingModel.polling_klimited() }

        // =======================================================================
        // Stochastic Petri Net models (from StochPetriNetModel.java)
        // =======================================================================
        register("spn_basic_closed") { StochPetriNetModel.spn_basic_closed() }
        register("spn_basic_open") { StochPetriNetModel.spn_basic_open() }
        register("spn_closed_fourplaces") { StochPetriNetModel.spn_closed_fourplaces() }
        register("spn_closed_twoplaces") { StochPetriNetModel.spn_closed_twoplaces() }
        register("spn_fourmodes") { StochPetriNetModel.spn_fourmodes() }
        register("spn_inhibiting") { StochPetriNetModel.spn_inhibiting() }
        register("spn_open_sevenplaces") { StochPetriNetModel.spn_open_sevenplaces() }
        register("spn_twomodes") { StochPetriNetModel.spn_twomodes() }

        // =======================================================================
        // State-dependent routing (from StateDepRoutingModel.java)
        // =======================================================================
        register("sdroute_open") { StateDepRoutingModel.sdroute_open() }
        register("sdroute_closed") { StateDepRoutingModel.sdroute_closed() }
        register("sdroute_twoclasses_closed") { StateDepRoutingModel.sdroute_twoclasses_closed() }

        // =======================================================================
        // Initial state models (from InitStateModel.java)
        // =======================================================================
        register("init_state_ps") { InitStateModel.init_state_ps() }
        register("init_state_fcfs_exp") { InitStateModel.init_state_fcfs_exp() }
        register("init_state_fcfs_nonexp") { InitStateModel.init_state_fcfs_nonexp() }

        // =======================================================================
        // Load-dependent models (from LoadDependentModel.java)
        // =======================================================================
        register("ld_class_dependence") { LoadDependentModel.ld_class_dependence() }
        register("ld_multiserver_fcfs") { LoadDependentModel.ld_multiserver_fcfs() }
        register("ld_multiserver_ps") { LoadDependentModel.ld_multiserver_ps() }
        register("ld_multiserver_ps_twoclasses") { LoadDependentModel.ld_multiserver_ps_twoclasses() }

        // =======================================================================
        // State probabilities models (from StateProbabilitiesModel.java)
        // =======================================================================
        register("statepr_aggr") { StateProbabilitiesModel.statepr_aggr() }
        register("statepr_aggr_large") { StateProbabilitiesModel.statepr_aggr_large() }
        register("statepr_sys_aggr") { StateProbabilitiesModel.statepr_sys_aggr() }
        register("statepr_sys_aggr_large") { StateProbabilitiesModel.statepr_sys_aggr_large() }
        register("statepr_allprobs_ps") { StateProbabilitiesModel.statepr_allprobs_ps() }
        register("statepr_allprobs_fcfs") { StateProbabilitiesModel.statepr_allprobs_fcfs() }

        // =======================================================================
        // CDF Response Time models (from CDFRespTModel.java)
        // =======================================================================
        register("cdf_respt_closed") { CDFRespTModel.cdf_respt_closed() }
        register("cdf_respt_closed_threeclasses") { CDFRespTModel.cdf_respt_closed_threeclasses() }
        register("cdf_respt_open_twoclasses") { CDFRespTModel.cdf_respt_open_twoclasses() }
        register("cdf_respt_distrib") { CDFRespTModel.cdf_respt_distrib() }
        register("cdf_respt_populations") { CDFRespTModel.cdf_respt_populations() }
    }

    private fun register(name: String, builder: () -> Network) {
        modelMap[name] = builder
    }

    /**
     * Get the singleton instance (for Java interop)
     */
    @JvmStatic
    fun getInstance(): ModelRegistry = this

    /**
     * Get a model by name.
     *
     * @param name Model name as string (e.g., "oqn_basic")
     * @return Network object configured for the model
     * @throws IllegalArgumentException if model name is not registered
     */
    @JvmStatic
    fun getModel(name: String): Network {
        val builder = modelMap[name]
            ?: throw IllegalArgumentException(
                "Model '$name' is not registered. Available models: ${modelMap.keys.joinToString(", ")}"
            )
        return builder()
    }

    /**
     * Get list of all registered model names
     */
    @JvmStatic
    fun getAvailableModels(): List<String> = modelMap.keys.toList()

    /**
     * Check if a model is registered
     */
    @JvmStatic
    fun hasModel(name: String): Boolean = modelMap.containsKey(name)
}
