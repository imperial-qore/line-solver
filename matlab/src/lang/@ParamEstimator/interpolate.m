function interpolate(self)
% interpolate all data across all the available timestamps
tunion = [];

for i=1:size(self.samples,1)
    for r=1:size(self.samples,2)
        for d=1:length(self.samples{i,r})
            if ~isempty(self.samples{i,r}{d})
                tunion = union(tunion, self.samples{i,r}{d}.t);
            end
        end
    end
end

for i=1:length(self.samplesAggr)
    for d=1:length(self.samplesAggr{i})
        if ~isempty(self.samplesAggr{i}{d})
            tunion = union(tunion, self.samplesAggr{i}{d}.t);
        end
    end
end

for i=1:size(self.samples,1)
    for r=1:size(self.samples,2)
        for d=1:length(self.samples{i,r})
            if ~isempty(self.samples{i,r}{d})
                % spline interpolation appears to be less
                % sensitive to the number of samples
                self.samples{i,r}{d}.data = interp1(self.samples{i,r}{d}.t, self.samples{i,r}{d}.data, tunion, 'spline');
                self.samples{i,r}{d}.t = tunion;
            end
        end
    end
end

for i=1:length(self.samplesAggr)
    for d=1:length(self.samplesAggr{i})
        if ~isempty(self.samplesAggr{i}{d})
            % spline interpolation appears to be less
            % sensitive to the number of samples
            self.samplesAggr{i}{d}.data = interp1(self.samplesAggr{i}{d}.t, self.samplesAggr{i}{d}.data, tunion, 'spline');
            self.samplesAggr{i}{d}.t = tunion;
        end
    end
end

end
