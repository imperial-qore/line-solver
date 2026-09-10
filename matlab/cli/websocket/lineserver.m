classdef lineserver < WebSocketServer
    
    properties
    end
    
    methods
        function obj = lineserver(varargin)
            %Constructor
            obj@WebSocketServer(varargin{:});
        end
    end
    
    methods (Access = protected)
        
        function onOpen(obj,conn,message)
            fprintf('%s\n',message)
        end
        
        function onTextMessage(obj,conn,message)           
            fname = tempname;
            r = regexp(message, '[\n]');
            args = strtrim(message(1:r)); % first row supplies args
            message = message((r+1):end);
            fid = fopen(fname,'w+');
            fprintf(fid, '%s', message); 
            fclose(fid);
            fprintf('Client model saved in: %s\n',fname);
            splitargs = split(args,',');   
            splitargs{1} = '--file'; % override target IP
            splitargs{2} = fname; % override target PORT            
            ret = linemcr(splitargs{:});
            conn.send(ret);
            conn.close();
        end
        
        function onBinaryMessage(obj,conn,bytearray)            
            fname = tempname;
            r = regexp(message, '[\n]');
            args = strtrim(message(1:r)); % first row supplies args
            message = message((r+1):end);
            fid = fopen(fname,'w+');
            fprintf(fid, '%s', message); 
            fclose(fid);
            fprintf('Client archive saved in: %s\n',fname);
            splitargs = split(args,',');   
            splitargs{1} = '--file'; % override target IP
            splitargs{2} = fname; % override target PORT            
            ret = linemcr(splitargs{:});
            conn.send(ret);
            conn.close();            
        end
        
        function onError(obj,conn,message)
            fprintf('%s\n',message)
        end
        
        function onClose(obj,conn,message)
            fprintf('%s\n',message)
            fprintf('Enter Q to stop the server at any time.\n')   
        end
    end
end

