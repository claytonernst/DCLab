classdef DCObject
    %This is an abstract class for DClab objects
    
    properties
        userData;
    end
    
    methods (Abstract)
        [list,sz] = displayProps(obj)
    end
    
    methods
        function display(obj)
            [list,sz] = displayProps(obj);
            if isequal(get(0,'FormatSpacing'),'compact')
                disp([inputname(1) ' =']);
                fprintf('    %s %s object\n',sz,class(obj))
                for i1=1:length(list)
                    disp(['       ' list{i1}]);
                end
            else
                disp(' ');
                disp([inputname(1) ' =']);
                disp(' ');
                fprintf('    %s %s object\n',sz,class(obj))
                for i1=1:length(list)
                    disp(['       ' list{i1}]);
                end
                disp(' ');
            end
        end
        function varargout = horzcat(varargin)
            error('Cannot horizontally concatenate %s objects',class(varargin{2}))
        end
        function varargout = vertcat(varargin)
            error('Cannot vertically concatenate %s objects',class(varargin{2}))
        end
        function out = cat(varargin)
            if isdouble(varargin{1}) && isequal(varargin{1},1)
                out = vertcat(varargin{2:end});
            elseif isdouble(varargin{1}) && isequal(varargin{1},2)
                out = horzcat(varargin{2:end});
            else
                error('Cannot concatenate %s objects in this way',class(varargin{2}))
            end
        end
        function obj = set(obj,varargin)
            ni = nargin;
            if mod(ni,2)==0
                error('Incorrect input.  Must be property/value pairs.')
            end
            for i1=1:(ni-1)/2
                obj.(varargin{2*i1-1}) = varargin{2*i1};
            end
        end
    end
    
end