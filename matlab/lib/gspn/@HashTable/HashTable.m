classdef HashTable < handle
    
    properties
        key_value_pair_count;
    end
    
    properties (Access = private)
        row_count;
        rows;
    end
    
    methods (Access = public)
        function this = HashTable(row_count)
            this.row_count = row_count;
            this.key_value_pair_count = 0;
            this.rows = cell(1,row_count);
        end
        
        function [exists,sequence_number] = search(this, marking)
            primary_key = this.primary_key(marking);
            if ~isempty(this.rows{primary_key})
                secondary_key = this.secondary_key(marking);
                list = this.rows{primary_key};
                while ~isempty(list)
                    if list.data(1) == secondary_key
                        exists = true;
                        sequence_number = list.data(2);
                        return
                    end
                    list = list.next;
                end
            end
            exists = false;
            sequence_number = -1;
        end
        
        function sequence_number = insert(this, marking)
            primary_key = this.primary_key(marking);
            secondary_key = this.secondary_key(marking);
            list = this.rows{primary_key};
            sequence_number = this.key_value_pair_count + 1;
            if ~isempty(list)
                insertAfter(dlnode([secondary_key, sequence_number]), list);
            else
                this.rows{primary_key} = dlnode([secondary_key,sequence_number]);
            end
            sprintf('Inserted at position %dx', primary_key);
            this.key_value_pair_count = sequence_number;
        end
        
        
        function key = primary_key(this, marking)
           key = mod(hf1(marking,17),this.row_count) + 1; 
        end
        
        function key = secondary_key(this, marking)
%            key = hf2(marking,17,32);
            key = hf4(marking);
        end
        
        function vector = vev(this)
           vector = this.rows; 
        end
    end
    
end