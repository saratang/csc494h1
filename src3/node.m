classdef node
    properties
        val
        next
    end
    methods
        function obj = node(varargin)
            if nargmin == 1
                obj.val = varargin{1};
            else
                obj.val = NaN;
            end
            obj.next = NaN;
        end
    end
end