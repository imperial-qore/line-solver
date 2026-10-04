classdef ItemEntry < Entry
    
    properties
        
        cardinality;
        popularity;
        
    end
    
    
    methods
        %Constructor
        function self = ItemEntry(model, name, cardinality, distribution)
            % SELF = LAYEREDNETWORKELEMENT(NAME)
            
            self@Entry(model, name);
            self.cardinality = cardinality;  % number of items
            if distribution.isDiscrete          
                self.popularity = distribution.copy;
            else
                line_error(mfilename,'A discrete popularity distribution is required.');                   
            end            
            
           
        end
        
        function self = on(self, parent)
            % SELF = ON(SELF, PARENT)
            
            parent.addEntry(self);
            self.parent = parent;
            
        end
       
       
        
    end
    
end