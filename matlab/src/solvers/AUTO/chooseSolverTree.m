function choice = chooseSolverTree(sn)
% CHOICE = CHOOSESOLVERTREE(SN)
%
% Learned solver selection: returns the name of the solver the fitted
% decision tree recommends, or '' when it abstains and the caller must
% fall back on the heuristic.
%
% GENERATED FILE - DO NOT EDIT BY HAND.
% Regenerate with line-test.git/algsel/export_tree.py after retraining;
% the feature map it indexes is sn_algsel_features.m, which must stay
% field for field identical to algsel/features.py.
%
% Trained on 500 instances, families cqn,fjc,fjo,mqn,oqn, learner dt_clf, stamp 2026-07-25.
%   gap closed vs single-best solver, leave-one-family-out: -13.936
%   gap closed vs single-best solver, stratified 5-fold by family: -5.465
% Abstention thresholds: leaf support >= 20, leaf purity >= 0.60,
% and every split feature inside its training range.

f = sn_algsel_features(sn);
choice = '';

    if f.f_log_state_space < 0.0 || f.f_log_state_space > 259.2710486794479
        choice = '';  % f_log_state_space outside training support
    elseif f.f_log_state_space <= 2.936129927635193
        if f.f_rho_max_open < 0.0 || f.f_rho_max_open > 10.0
            choice = '';  % f_rho_max_open outside training support
        elseif f.f_rho_max_open <= 1.710871160030365
            if f.f_n_closed < 0.0 || f.f_n_closed > 3.0
                choice = '';  % f_n_closed outside training support
            elseif f.f_n_closed <= 0.5
                if f.f_bottleneck_ratio < 1.0004359670723832 || f.f_bottleneck_ratio > 32.59065337673579
                    choice = '';  % f_bottleneck_ratio outside training support
                elseif f.f_bottleneck_ratio <= 1.2381284832954407
                    if f.f_demand_cv < 0.0132959652844532 || f.f_demand_cv > 0.441739047614697
                        choice = '';  % f_demand_cv outside training support
                    elseif f.f_demand_cv <= 0.3455698937177658
                        % leaf abstains: support=10 purity=0.98
                        choice = '';
                    else
                        % leaf abstains: support=13 purity=0.55
                        choice = '';
                    end
                else
                    if f.f_bottleneck_ratio < 1.0004359670723832 || f.f_bottleneck_ratio > 32.59065337673579
                        choice = '';  % f_bottleneck_ratio outside training support
                    elseif f.f_bottleneck_ratio <= 1.4162747263908386
                        % leaf abstains: support=10 purity=0.84
                        choice = '';
                    else
                        % leaf abstains: support=86 purity=0.54
                        choice = '';
                    end
                end
            else
                if f.f_log_state_space < 0.0 || f.f_log_state_space > 259.2710486794479
                    choice = '';  % f_log_state_space outside training support
                elseif f.f_log_state_space <= 1.7299033999443054
                    if f.f_log_state_space < 0.0 || f.f_log_state_space > 259.2710486794479
                        choice = '';  % f_log_state_space outside training support
                    elseif f.f_log_state_space <= 1.0168417096138
                        % support=62 purity=0.99
                        choice = 'ctmc';
                    else
                        % support=26 purity=0.79
                        choice = 'ctmc';
                    end
                else
                    % leaf abstains: support=17 purity=0.67
                    choice = '';
                end
            end
        else
            if f.f_max_servers < 1.0 || f.f_max_servers > 5.0
                choice = '';  % f_max_servers outside training support
            elseif f.f_max_servers <= 2.5
                if f.f_demand_mean < 0.4782704888654595 || f.f_demand_mean > 17.979001499617286
                    choice = '';  % f_demand_mean outside training support
                elseif f.f_demand_mean <= 1.549913763999939
                    % leaf abstains: support=12 purity=0.97
                    choice = '';
                else
                    % leaf abstains: support=11 purity=0.74
                    choice = '';
                end
            else
                if f.f_arrival_rate < 0.0 || f.f_arrival_rate > 3.5478173726428768
                    choice = '';  % f_arrival_rate outside training support
                elseif f.f_arrival_rate <= 1.2968124151229858
                    if f.f_frac_hypo < 0.0 || f.f_frac_hypo > 1.0
                        choice = '';  % f_frac_hypo outside training support
                    elseif f.f_frac_hypo <= 0.24444445222616196
                        % leaf abstains: support=11 purity=0.79
                        choice = '';
                    else
                        % leaf abstains: support=10 purity=0.84
                        choice = '';
                    end
                else
                    if f.f_frac_sched_siro < 0.0 || f.f_frac_sched_siro > 0.8333333333333334
                        choice = '';  % f_frac_sched_siro outside training support
                    elseif f.f_frac_sched_siro <= 0.0714285746216774
                        % leaf abstains: support=19 purity=0.97
                        choice = '';
                    else
                        % support=20 purity=0.88
                        choice = 'mva_exact';
                    end
                end
            end
        end
    else
        if f.f_frac_sched_ps < 0.0 || f.f_frac_sched_ps > 0.75
            choice = '';  % f_frac_sched_ps outside training support
        elseif f.f_frac_sched_ps <= 0.0714285746216774
            if f.f_frac_sched_lcfspr < 0.0 || f.f_frac_sched_lcfspr > 1.0
                choice = '';  % f_frac_sched_lcfspr outside training support
            elseif f.f_frac_sched_lcfspr <= 0.0833333358168602
                if f.f_jobs_per_chain < 0.0 || f.f_jobs_per_chain > 40.0
                    choice = '';  % f_jobs_per_chain outside training support
                elseif f.f_jobs_per_chain <= 2.416666626930237
                    if f.f_mean_phases < 1.0 || f.f_mean_phases > 5.0
                        choice = '';  % f_mean_phases outside training support
                    elseif f.f_mean_phases <= 2.109890103340149
                        % leaf abstains: support=13 purity=0.82
                        choice = '';
                    else
                        % leaf abstains: support=11 purity=0.97
                        choice = '';
                    end
                else
                    if f.f_mean_phases < 1.0 || f.f_mean_phases > 5.0
                        choice = '';  % f_mean_phases outside training support
                    elseif f.f_mean_phases <= 2.450000047683716
                        % support=30 purity=1.00
                        choice = 'ssa';
                    else
                        % leaf abstains: support=10 purity=0.59
                        choice = '';
                    end
                end
            else
                if f.f_frac_sched_siro < 0.0 || f.f_frac_sched_siro > 0.8333333333333334
                    choice = '';  % f_frac_sched_siro outside training support
                elseif f.f_frac_sched_siro <= 0.125
                    % leaf abstains: support=10 purity=0.85
                    choice = '';
                else
                    % leaf abstains: support=13 purity=0.79
                    choice = '';
                end
            end
        else
            if f.f_rho_max_open < 0.0 || f.f_rho_max_open > 10.0
                choice = '';  % f_rho_max_open outside training support
            elseif f.f_rho_max_open <= 9.983254432678223
                if f.f_demand_max < 0.8745904620935151 || f.f_demand_max > 42.85554925970361
                    choice = '';  % f_demand_max outside training support
                elseif f.f_demand_max <= 2.4841485023498535
                    % leaf abstains: support=10 purity=0.77
                    choice = '';
                else
                    if f.f_frac_sched_inf < 0.0 || f.f_frac_sched_inf > 0.8333333333333334
                        choice = '';  % f_frac_sched_inf outside training support
                    elseif f.f_frac_sched_inf <= 0.0714285746216774
                        % support=28 purity=0.64
                        choice = 'fld';
                    else
                        % leaf abstains: support=53 purity=0.54
                        choice = '';
                    end
                end
            else
                % leaf abstains: support=15 purity=0.76
                choice = '';
            end
        end
    end
end
