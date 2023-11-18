function simodel = getStateIndependentModel(self)
sn = self.getStruct;
if any(sn.nodetype == NodeType.Cache)
    [~,Pnodes] = self.getRoutingMatrix;

    simodel = Network('simodel');
    for ind=1:self.getNumberOfNodes
        switch class(self.nodes{ind})
            case 'Cache'
                Pcs = zeros(self.getNumberOfClasses);
                hitClass = full(self.nodes{ind}.getHitClass);
                missClass = full(self.nodes{ind}.getMissClass);

                actualHitProb = self.nodes{ind}.getHitRatio;
                actualMissProb = self.nodes{ind}.getMissRatio;

                for r=1:length(hitClass)
                    if hitClass(r)>0
                        Pcs(r,hitClass(r)) = actualHitProb(r);
                    end
                end
                for r=1:length(missClass)
                    if missClass(r)>0
                        Pcs(r,missClass(r)) = actualMissProb(r);
                    end
                end
                for r=1:size(Pcs,1)
                    if sum(Pcs(r,:)) == 0
                        Pcs(r,r) = 1;
                    end
                end
                staticcache = ClassSwitch(simodel, 'StaticCache', Pcs);
            otherwise
                simodel.addNode(self.nodes{ind});
        end
    end

    for r=1:self.getNumberOfClasses
        simodel.addJobClass(self.classes{r});
    end

    simodel.linkFromNodeRoutingMatrix(Pnodes);

    % sanitize disabled classes
    for ind=1:self.getNumberOfNodes
        switch class(self.nodes{ind})
            case 'Cache'
                for r=1:length(staticcache.output.outputStrategy)
                    if isempty(staticcache.output.outputStrategy{r})
                        staticcache.output.outputStrategy{r} = {[],'Disabled'};
                    end
                end
                for r=(length(staticcache.output.outputStrategy)+1) : simodel.getNumberOfClasses
                    staticcache.output.outputStrategy{r} = {[],'Disabled'};
                end
        end
    end
end
end