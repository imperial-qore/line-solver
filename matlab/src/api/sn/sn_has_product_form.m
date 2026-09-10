%{ @file sn_has_product_form.m
 %  @brief Checks if the network has product-form solution
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network has product-form solution
 %
 % @details
 % Returns true if the network satisfies conditions for product-form equilibrium distribution.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_product_form(sn)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>bool<td>True if the network has product-form solution
 % </table>
%}
function bool = sn_has_product_form(sn)

bool = all(sn.sched==SchedStrategy.INF | sn.sched==SchedStrategy.PS | sn.sched==SchedStrategy.FCFS | sn.sched==SchedStrategy.LCFSPR | sn.sched==SchedStrategy.LCFS | sn.sched==SchedStrategy.EXT);
bool = bool & ~sn_has_multi_class_heter_fcfs(sn);
bool = bool & ~sn_has_priorities(sn);
bool = bool & ~sn_has_fork_join(sn);
bool = bool & ~sn_has_sd_routing(sn);
% BCMP asks for infinite buffers. Nothing here read sn.cap/sn.classcap/
% sn.droprule, so a BAS-blocked station (cqn_bas_blocking) or any binding
% finite buffer passed the gate and the console reported hasProductForm=1
% on a network whose truncation couples the station occupancies.
bool = bool & ~sn_has_blocking(sn);
% BCMP type 1 asks the FCFS service to be exponential. SN_HAS_MULTI_CLASS_HETER_FCFS
% compares the class MEANS only, so a class-homogeneous Erlang, hyper-exponential
% or deterministic FCFS station passed this gate and was dispatched to exact MVA,
% which reads the means alone and returns the exponential answer with no warning.
iset = find(sn.sched == SchedStrategy.FCFS);
for idx=1:length(iset)
    i = iset(idx);
    icset = isfinite(sn.scv(i,:)) & sn.scv(i,:)>0;
    bool = bool & all(sn.scv(i,icset) > 1-GlobalConstants.FineTol) & all(sn.scv(i,icset) < 1+GlobalConstants.FineTol);
end
end