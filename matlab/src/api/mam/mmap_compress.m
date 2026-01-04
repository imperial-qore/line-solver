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
if nargin<1%~exist('config','var')
    config = struct;
end
if ~isfield(config,'method')
    config.method = 'default';
end

K = length(MMAP)-2;
switch config.method
    case {'default','mixture','mixture.order1'}
        lambda = mmap_lambda(MMAP);
        AMAPs = mmap_maps(MMAP);
        for k=1:K
            AMAPs{k} =  mmpp2_fit1(map_mean(AMAPs{k}),map_scv(AMAPs{k}),map_skew(AMAPs{k}),map_idc(AMAPs{k}));
        end
        MMAP = mmap_mixture(lambda/sum(lambda), AMAPs);
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
