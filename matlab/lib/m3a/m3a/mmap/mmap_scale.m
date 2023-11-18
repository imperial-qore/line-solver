function [SCALED] = mmap_scale(MMAP, M, maxIter)
% Changes the mean inter-arrival time of an MMAP.
% INPUT
% - MMAP: MMAP to be scaled
% - M: new mean
if nargin<3
    maxIter=30;
end
C = length(MMAP)-2;
if length(M)==1

    MOLD = map_mean(MMAP);

    ratio = MOLD/M;

    SCALED = cell(1,2+C);
    SCALED{1} = MMAP{1} * ratio;
    SCALED{2} = MMAP{2} * ratio;

    for c = 1:C
        SCALED{2+c} = MMAP{2+c} * ratio;
    end
else
    SCALED{1} = MMAP{1};
    SCALED{2} = 0*SCALED{1};
    l = mmap_count_lambda(MMAP);
    for c = 1:C
        if l(c) > 0
            SCALED{2+c} = MMAP{2+c} * (1/M(c))/l(c);
            SCALED{2} = SCALED{2} + SCALED{2+c};
        else
            SCALED{2+c} = MMAP{2+c} * 0;
            SCALED{2} = SCALED{2} + SCALED{2+c};            
        end
    end
    MMAP = mmap_normalize(SCALED);
    return
    % the previous assignment is heuristic because it also affects the
    % other classes, we now refine it
    try
        if maxiter > 0
            options = optimset();
            options.Display = 'off';
            options.tolFun = 1e-2;
            options.MaxIter = maxIter;
            %x=fmincon(@(X) objfun(X,M,MMAP),ones(1,C),[],[],[],[],1e-6+zeros(1,C),[],[],options);
            x = fminsearchbnd(@(X) objfun(X,M,MMAP),ones(1,C),1e-6+zeros(1,C),[],options);
            SCALED{1} = MMAP{1};
            SCALED{2} = 0*SCALED{1};
            for c = 1:C
                SCALED{2+c} = MMAP{2+c} * x(c);
                SCALED{2} = SCALED{2} + SCALED{2+c};
            end
            SCALED = mmap_normalize(SCALED);
        end
    catch
        error('The input MMAP is invalid.');
    end
end
end

function f = objfun(x,M,MMAP)
f = 0;
C = length(MMAP)-2;
SCALED{1} = MMAP{1};
SCALED{2} = 0*SCALED{1};
for c = 1:C
    SCALED{2+c} = MMAP{2+c} * x(c);
    SCALED{2} = SCALED{2} + SCALED{2+c};
end
SCALED = mmap_normalize(SCALED);
l = mmap_count_lambda(SCALED);
for c=1:C
    f= f + norm((1/M(c))-l(c));
end
end

