function b = isTransitionGreater(IM1, IM2, P1, P2, W1, W2, D1, D2)
% ISTRANSITIONGREATER determines whether T1 has a higher priority than T2
% b is true if T1 > T2, false otherwise
switch(IM1)
    % T1 is immediate     
    case(1)
        switch(IM2)
            % T2 is immediate             
            case(1)
                if P1 > P2
                    b = true;
                elseif P1 < P2
                    b = false;
                else
                    b = W1 > W2;
                end
            % T2 is timed
            case(0)
                b = true;
        end
    % T1 is timed         
    case(0)
        switch(IM2)
            % T2 is immediate             
            case(1)
                b = false;
            % T2 is timed
            case(0)
                
                
        end
end

end