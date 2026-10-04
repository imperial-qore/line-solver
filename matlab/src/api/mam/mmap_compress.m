%{ @file mmap_compress.m
 %  @brief Compresses an MMAP into a smaller representation
 %
 %  @author LINE Development Team
%}

%{
 % @brief Compresses a Marked MAP into a smaller representation
 %
 % @details
 % This function compresses an MMAP (Marked Markovian Arrival Process) into
 % a smaller representation using various compression methods including
 % mixture fitting, MAMAP2, and M3PP approaches.
 %
 % @par Syntax:
 % @code
 % MMAP = mmap_compress(MMAP)
 % MMAP = mmap_compress(MMAP, config)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>MMAP<td>Original Marked Markovian Arrival Process
 % <tr><td>config<td>(Optional) Configuration struct with compression method
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>MMAP<td>Compressed MMAP
 % </table>
%}
function MMAP = mmap_compress(MMAP, config)
if nargin<2 % config not supplied: use the default method
    config = struct;
end
if ~isfield(config,'method')
    config.method = 'default';
end

K = length(MMAP)-2;
switch config.method
    case {'default','mixture','mixture.order1'}
        % Order-1 mixture (M3A). K components, one per class, recombined by
        % mmap_mixture with the class probabilities p_c as mixing weights.
        %
        % Component c must carry the law of the inter-arrival time CONDITIONED
        % ON THE ARRIVAL THAT ENDS IT BEING OF CLASS c: mmap_mixture marks the
        % arrival LEAVING component c with class c (mmap_mixture.m:24-31), so
        % the class of an arrival and the interval preceding it are both
        % governed by the component active during that interval. That
        % conditional law is the class-c BACKWARD moment set B(c,1:3), i.e.
        % E[T^k | class of the ending arrival = c] (mmap_backward_moment.m).
        % It is NOT the forward moment (E[T^k | class of the STARTING arrival
        % = c]) and it is NOT the class-c marginal MAP from mmap_maps (whose
        % mean is 1/lambda_c, the time between successive class-c arrivals --
        % mixing those with weights lambda_c/Lambda inflates the mean to
        % K/Lambda). This mirrors 'mixture.order2', which conditions on the
        % (last,next) class PAIR and fits the CROSS moments
        % (mmap_mixture_fit.m:1-5), order 1 simply drops the "last" index.
        %
        % PRESERVED exactly: aggregate moments 1..3, via the M3A mixture law
        % M_k = sum_c B(k,c)*p_c stated in mmap_backward_moment.m:5-11 (M1 is
        % always exact because aph2_adjust never alters M1; M2/M3 are exact
        % when APH(2)-feasible); the class probabilities p_c and hence the
        % per-class rates lambda_c = p_c/M1; marking consistency
        % D1 = sum_c D1^(c); and MAP feasibility (mmap_mixture normalizes).
        %
        % LOST by construction: every autocorrelation. mmap_mixture re-enters
        % each component at its map_pie on every arrival (mmap_mixture.m:20-22),
        % so the intervals are i.i.d. and the result is a RENEWAL process:
        % acf -> 0, IDC -> the SCV-determined renewal value, and the class
        % sequence becomes i.i.d. (sigma(i,j) -> p_j). Retaining sigma is
        % exactly what 'mixture.order2' buys with its K^2 components.
        p = mmap_pc(MMAP);
        B = mmap_backward_moment(MMAP, [1 2 3], 1);
        AMAPs = cell(1,K);
        for k=1:K
            if p(k) <= GlobalConstants.Zero
                % Class c never arrives, so B(k,:) is an 0/0 normalization.
                % The component carries zero mixture weight: any proper MAP
                % leaves the result unchanged.
                AMAPs{k} = map_exponential(1);
            else
                AMAPs{k} = aph2_fit(B(k,1), B(k,2), B(k,3));
            end
        end
        MMAP = mmap_mixture(p, AMAPs);
    case 'mixture.order2'
        MMAP = mmap_mixture_fit_mmap(MMAP);
    case 'mamap2'
        MMAP = mmap_normalize(MMAP);
        MMAP = mamap2m_fit_mmap(MMAP);
    case 'mamap2.fb'
        MMAP = mmap_normalize(MMAP);
        MMAP = mamap2m_fit_gamma_fb_mmap(MMAP);
    case 'm3pp.approx_cov'
        MMAP = m3pp2m_fitc_theoretical(MMAP, 'approx_cov', 1, 1e6); %derivest
    case 'm3pp.approx_ag'
        MMAP = m3pp2m_fitc_theoretical(MMAP, 'approx_ag', 1, 1e6); %derivest
    case 'm3pp.exact_delta'
        MMAP = m3pp2m_fitc_theoretical(MMAP, 'exact_delta', 1, 1e6); %derivest
    case 'm3pp.approx_delta'
        MMAP = m3pp2m_fitc_theoretical(MMAP, 'approx_delta', 1, 1e6); %derivest
end
MMAP = mmap_normalize(MMAP);
end
