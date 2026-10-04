function tf = jmtJmvaIsClosedOnly(method)
% TF = JMTJMVAISCLOSEDONLY(METHOD)
%
% True for the JMVA algorithms that solve a CLOSED product-form network only:
% RECAL, CoMoM, Chow, Bard-Schweitzer (both 'jmva.bs' and 'jmva.amva'), AQL,
% Linearizer and De Souza-Muntz Linearizer. Measured against JMT 1.2.x: each
% answers an open or a mixed model with
%   jmt.common.exception.UnsupportedModelException: The selected solver cannot
%   handle open classes, please choose another.
% and a load-dependent one with the same exception naming load-dependent
% stations, while the exact MVA engine behind 'jmva' and 'jmva.mva' serves both.
%
% The three consequences of that split are declared in one place each:
% getMethodFeatureSet drops OpenClass and LoadDependence for these names,
% jmtMethodRefusal refuses a multi-server station for them (a server count has
% no feature name), and writeJMVA raises the same sentence when asked by name.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

tf = any(strcmpi(method, {'jmva.amva','jmva.recal','jmva.comom','jmva.chow', ...
    'jmva.bs','jmva.aql','jmva.lin','jmva.dmlin'}));
end
