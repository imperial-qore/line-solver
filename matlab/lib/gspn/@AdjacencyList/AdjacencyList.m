classdef AdjacencyList < handle
    
    properties (Access = private)
        edges;
    end
    
    methods (Access = public)
        function this = AdjacencyList(edge_count)
            this.edges = cell(1,edge_count);
        end
        
        function insert(this, source, destination, value)
            if source > length(this.edges) || isempty(this.edges{source})
                this.edges{source} = dlnode([destination,value]);
            else
                insertAfter(dlnode([destination,value]), this.edges{source});
            end
        end
        
        function edges = getEdges(this)
           edges = this.edges; 
        end
    end
    
end