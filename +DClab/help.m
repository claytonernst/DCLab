function t = help(varargin)

if nargin>1
    error('Help only supports one topic')
elseif nargin ==0
    builtin('help')
end

if ~ischar(varargin{1})
    error('Argument to help must be a string.')
end

topic = varargin{1};

if length(topic)>6 && all(lower(topic(1:6))=='dclab.')
    tmp = help(topic);
else

    tmp = help(['DClab.' topic]);

    if isempty(tmp)
        tmp = help(topic);
    end
    if isempty(tmp) && nargout==0
        help(topic);
        return
    end

end

if nargout==1
    t = tmp;
else
    disp(tmp);
end