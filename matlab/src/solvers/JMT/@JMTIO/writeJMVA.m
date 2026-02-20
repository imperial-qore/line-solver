function outputFileName = writeJMVA(sn, outputFileName, options)
% FNAME = WRITEJMVA(SN, FNAME, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%if ~self.model.hasProductFormSolution
%    line_error(mfilename,'JMVA requires the model to have a product-form solution.');
%end

%if self.model.hasClassSwitching
%    line_error(mfilename,'JMVA does not support class switching.');
%end

mvaDoc = com.mathworks.xml.XMLUtils.createDocument('model');
mvaElem = mvaDoc.getDocumentElement;

xmlnsXsi = 'http://www.w3.org/2001/XMLSchema-instance';
mvaElem.setAttribute('xmlns:xsi', xmlnsXsi);
mvaElem.setAttribute('xsi:noNamespaceSchemaLocation', 'JMTmodel.xsd');

algTypeElement = mvaDoc.createElement('algType');
switch options.method
    case {'jmva.recal'}
        if max(sn.nservers(isfinite(sn.nservers))) > 1
            line_error(mfilename,sprintf('%s does not support multi-server stations.',options.method));
        end
        algTypeElement.setAttribute('name','RECAL');
    case {'jmva.chow'}
        if max(sn.nservers(isfinite(sn.nservers))) > 1
            line_error(mfilename,sprintf('%s does not support multi-server stations.',options.method));
        end
        algTypeElement.setAttribute('name','Chow');
    case {'jmva.bs','jmva.amva'}
        if max(sn.nservers(isfinite(sn.nservers))) > 1
            line_error(mfilename,sprintf('%s does not support multi-server stations.',options.method));
        end
        algTypeElement.setAttribute('name','Bard-Schweitzer');
    case {'jmva.aql'}
        if max(sn.nservers(isfinite(sn.nservers))) > 1
            line_error(mfilename,sprintf('%s does not support multi-server stations.',options.method));
        end
        algTypeElement.setAttribute('name','AQL');
    case {'jmva.lin'}
        if max(sn.nservers(isfinite(sn.nservers))) > 1
            line_error(mfilename,sprintf('%s does not support multi-server stations.',options.method));
        end
        algTypeElement.setAttribute('name','Linearizer');
    case {'jmva.dmlin'}
        if max(sn.nservers(isfinite(sn.nservers))) > 1
            line_error(mfilename,sprintf('%s does not support multi-server stations.',options.method));
        end
        algTypeElement.setAttribute('name','De Souza-Muntz Linearizer');
    %case {'jmva.ls'}
    %    algTypeElement.setAttribute('name','Logistic Sampling');
    otherwise
        algTypeElement.setAttribute('name','MVA');
end
algTypeElement.setAttribute('tolerance','1.0E-7');
algTypeElement.setAttribute('maxSamples',num2str(options.samples));

%%%%%%%%%%
M = sn.nstations;    %number of stations
NK = sn.njobs';  % initial population per class
C = sn.nchains;
SCV = sn.scv;

% determine service times
ST = 1./sn.rates;
ST(isnan(sn.rates))=0;
SCV(isnan(SCV))=1;

[Lchain,STchain,Vchain,alpha,Nchain,SCVchain] = sn_get_demands_chain(sn);

%%%%%%%%%%
parametersElem = mvaDoc.createElement('parameters');
classesElem = mvaDoc.createElement('classes');
classesElem.setAttribute('number',num2str(sn.nchains));
stationsElem = mvaDoc.createElement('stations');
stationsElem.setAttribute('number',num2str(sn.nstations - sum(sn.nodetype == NodeType.Source)));
refStationsElem = mvaDoc.createElement('ReferenceStation');
refStationsElem.setAttribute('number',num2str(sn.nchains));
algParamsElem = mvaDoc.createElement('algParams');

sourceid = sn.nodetype == NodeType.Source;
for c=1:sn.nchains
    if isfinite(sum(sn.njobs(sn.chains(c,:))))
        classElem = mvaDoc.createElement('closedclass');
        classElem.setAttribute('population',num2str(Nchain(c)));
        classElem.setAttribute('name',sprintf('Chain%02d',c));
    else
        classElem = mvaDoc.createElement('openclass');
        classElem.setAttribute('rate',num2str(sum(sn.rates(sourceid,sn.chains(c,:)))));
        classElem.setAttribute('name',sprintf('Chain%02d',c));
    end
    classesElem.appendChild(classElem);
end

isLoadDep = false(1,sn.nstations);
for i=1:sn.nstations
    switch sn.nodetype(sn.stationToNode(i))
        case NodeType.Delay
            statElem = mvaDoc.createElement('delaystation');
            statElem.setAttribute('name',sn.nodenames{sn.stationToNode(i)});
        case NodeType.Queue
            if sn.nservers(i) == 1
                isLoadDep(i) = false;
                statElem = mvaDoc.createElement('listation');
                statElem.setAttribute('name',sn.nodenames{sn.stationToNode(i)});
                statElem.setAttribute('servers',num2str(1));
            else
                isLoadDep(i) = true;
                statElem = mvaDoc.createElement('ldstation');
                statElem.setAttribute('name',sn.nodenames{sn.stationToNode(i)});
                statElem.setAttribute('servers',num2str(1));
            end
        otherwise
            continue
    end
    srvTimesElem = mvaDoc.createElement('servicetimes');
    for c=1:sn.nchains
        if isLoadDep(i)
            statSrvTimeElem = mvaDoc.createElement('servicetimes');
            statSrvTimeElem.setAttribute('customerclass',sprintf('Chain%02d',c));
            ldSrvString = num2str(STchain(i,c));
            if any(isinf(NK))
                line_error(mfilename,'JMVA does not support open classes in load-dependent models;');
            end

            for n=2:sum(NK)
                ldSrvString = sprintf('%s;%s',ldSrvString,num2str(STchain(i,c)/min( n, sn.nservers(i) )));
            end
            statSrvTimeElem.appendChild(mvaDoc.createTextNode(ldSrvString));
            srvTimesElem.appendChild(statSrvTimeElem);
        else
            statSrvTimeElem = mvaDoc.createElement('servicetime');
            statSrvTimeElem.setAttribute('customerclass',sprintf('Chain%02d',c));
            statSrvTimeElem.appendChild(mvaDoc.createTextNode(num2str(STchain(i,c))));
            srvTimesElem.appendChild(statSrvTimeElem);
        end
    end
    statElem.appendChild(srvTimesElem);
    visitsElem = mvaDoc.createElement('visits');
    for c=1:sn.nchains
        statVisitElem = mvaDoc.createElement('visit');
        statVisitElem.setAttribute('customerclass',sprintf('Chain%02d',c));
        if STchain(i,c) > 0
            val = Lchain(i,c) ./ STchain(i,c);
        else
            val = 0;
        end
        statVisitElem.appendChild(mvaDoc.createTextNode(num2str(val)));
        visitsElem.appendChild(statVisitElem);
    end
    statElem.appendChild(visitsElem);

    stationsElem.appendChild(statElem);
end

refstatchain = zeros(C,1);
for c=1:sn.nchains
    inchain = sn.inchain{c};
    refstatchain(c) = sn.refstat(inchain(1));
end

for c=1:sn.nchains
    classRefElem = mvaDoc.createElement('Class');
    classRefElem.setAttribute('name',sprintf('Chain%02d',c));
    % For open chains, the refstat is Source which is excluded from JMVA output.
    % Use the first non-Source station as the reference station instead.
    refIdx = refstatchain(c);
    refNodeType = sn.nodetype(sn.stationToNode(refIdx));
    if refNodeType == NodeType.Source
        for ii=1:sn.nstations
            nt = sn.nodetype(sn.stationToNode(ii));
            if nt ~= NodeType.Source && nt ~= NodeType.Sink
                refIdx = ii;
                break;
            end
        end
    end
    classRefElem.setAttribute('refStation',sn.nodenames{sn.stationToNode(refIdx)});
    refStationsElem.appendChild(classRefElem);
end

compareAlgsElem = mvaDoc.createElement('compareAlgs');
compareAlgsElem.setAttribute('value','false');
algParamsElem.appendChild(algTypeElement);
algParamsElem.appendChild(compareAlgsElem);

parametersElem.appendChild(classesElem);
parametersElem.appendChild(stationsElem);
parametersElem.appendChild(refStationsElem);
mvaElem.appendChild(parametersElem);
mvaElem.appendChild(algParamsElem);

xmlwrite(outputFileName, mvaDoc);
end
