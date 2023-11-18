classdef TypedQueue < handle
   
   % Queue - strongly-typed Queue collection
   %
   % Properties:
   %
   %   Type (string)
   %
   % Methods:
   %
   %   Queue(type)
   %   display
   %   size
   %   isempty
   %   clear
   %   contains(obj)
   %   offer(obj)
   %   remove(obj)
   %   peek   - returns [] if queue is empty
   %   poll   - returns [] if queue is empty
   %   values - returns contents in a cell array
   %
   % Notes:
   %
   % Compatible classes must overload eq() for object-to-object comparisons.
   %
   % Example:
   %
   % q = Queue('Widget')   
   %
   % q.offer(RedWidget(1))
   % q.offer(RedWidget(3))
   % q.offer(RedWidget(2))
   % q.offer(BlueWidget(2))
   % q.offer(BlueWidget(2))
   %
   % q.size
   % q.remove(RedWidget(2));
   % q.size
   % q.remove(BlueWidget(2));
   % q.size
   %
   % q.peek
   % q.poll
   % q.size
   %
   % Author: dimitri.shvorob@gmail.com, 3/15/09 
   
   properties (GetAccess = protected, SetAccess = protected, Hidden = true)
       Elements
   end
   
   properties (SetAccess = protected)
       Type
   end
         
   methods 
       
       function[obj] = TypedQueue(type)
           if ~ischar(type)
              throw(MException('Queue:constructorInvalidType','??? ''type'' must be a valid class name.'))
           end   
           obj.Elements = {};
           obj.Type = type;
       end
       
       function disp(obj)
           disp([class(obj) '<' obj.Type '> (head on top)'])
           if ~obj.isempty
              for i = 1:obj.size
                  disp(obj.Elements{i})
              end   
           else
              disp([])
           end
       end
              
       function[out] = size(obj)
           out = length(obj.Elements);
       end
       
       function[out] = values(obj)
           out = obj.Elements;
       end
              
       function[out] = isempty(obj)
           out = obj.size == 0;
       end
       
       function[obj] = clear(obj)
           obj.Elements = {};
       end
       
       function[out] = contains(obj,e)
           out = false;
           for i = 1:obj.size
               if e == obj.Elements{i}
                  out = true;
                  break
               end
           end           
       end
       
       function[obj] = offer(obj,e) 
           if ~isa(e,obj.Type)
              throw(MException('Queue:offerInvalidType','??? Invalid type.'))
           end
           if isempty(obj.Elements)
              obj.Elements = {e};
           else
              obj.Elements{end+1} = e;
           end
       end   
       
       function[obj] = remove(obj,e)
           if length(e) > 1
              throw(MException('Queue:removeMultiple','??? Cannot remove multiple elements at once.'))
           end 
           if ~isa(e,obj.Type)
              throw(MException('Queue:removeInvalidType','??? Invalid type.'))
           end
           if ~isempty(obj.Elements)
              k = [];
              for i = 1:obj.size
                  if e == obj.Elements{i}
                     k = [k i];  %#ok
                  end
              end
              if ~isempty(k)
                  obj.Elements(k) = [];
              end
           end
       end
       
       function[out] = peek(obj)
           if ~obj.isempty
               out = obj.Elements{1};
           else
               out = [];
           end    
       end
       
       function[out] = poll(obj)
           if ~obj.isempty
               out = obj.Elements{1};
               obj.Elements(1) = [];
           else
               out = [];
           end    
       end
           
   end   
    
end