function [bool message] = isLegalDomain(modelDomain,nameInMessage)
%ISLEGALDOMAIN checks that DOMAIN has appropriate .name and .range property/values
%
%   BOOL = ISLEGALDOMAIN(DOMAIN) returns a true if there are no problems
%   with DOMAIN. We check that DOMAIN has fields .name and .range. That
%   .range contains 1-by-2 vectors that describe nonempty intervals, and
%   that the names are unique.
%
%   [BOOL MESSAGE] = ISLEGALDOMAIN(DOMAIN) additionally returns a string
%   describing any found problems with DOMAIN. If BOOL==true, MESSAGE will
%   be empty.
%
%   [...] = ISLEGALDOMAIN(DOMAIN,NAMEINMESSAGE) allows you to supply a name
%   that will be used to refer to the domain in any generated messages. For
%   example, this may be 'NEWDOMAIN', or 'model''s domain'. If not
%   supplied, 'DOMAIN' will be used.  

if nargin == 1
  nameInMessage = 'DOMAIN';
end

message = '';
bool = true;

if ~isstruct(modelDomain) || ~isequal(numel(modelDomain),length(modelDomain)) || (numel(modelDomain) > 0 && size(modelDomain,2)~=1)
  bool = false;
  message = [nameInMessage ' must an n-by-1 structure array'];
elseif ~isfield(modelDomain,'name') || ~isfield(modelDomain,'range')
  bool = false;
  message = ['Insufficient fields: ' nameInMessage ' must have fields ''name'' and ''range'''];
else
  %check the contents of the structure array
  names = {modelDomain.name}';
  if ~iscellstr(names) || ~isequal(size(names),[numel(modelDomain),1])
    bool = false;
    message = ['The .name field of ' nameInMessage ' must contain single row chars'];
  end
  if length(names) ~= length(unique(names))
    bool = false;
    message = [nameInMessage ' must contain parameters with distinct names'];
  end
  
  ranges = {modelDomain.range}';
  if ~all(cellfun('size',ranges,2)==2)
    bool = false;
    message = ['The .range field of ' nameInMessage ' must contain 1x2 vectors'];
  else
    ranges = vertcat(modelDomain.range);
    if ~isnumeric(ranges) || ~isequal(size(ranges),[numel(modelDomain),2])
      bool = false;
      message = ['The .range field of ' nameInMessage ' must contain n 1x2 numeric vectors'];
    end
    if any(ranges(:,2) < ranges(:,1))
      bool = false;
      message = [nameInMessage ' must describe a nonempty set'];
    end
  end
end
