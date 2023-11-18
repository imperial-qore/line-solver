function [TF] = mmap_issym(MMAP)
% Returns 1 if the MMAP is symbolic, 0 otherwise.

TF = 0;

for i = 1:length(MMAP)
    if strcmp(class(MMAP{i}),'sym')
    	TF = 1;
    	return
    end
end