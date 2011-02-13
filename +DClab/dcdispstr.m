function dcdispstr(str,hand,iserr)
% str is a single row char array or an error structure
% hand is empty or a valid graphics handle
% iserr is a bool, indicating whether this is an error or standard display

if ~iserr
  if isempty(hand)
    if isstruct(str)
      disp(str.message)
    else
      disp(str)
    end
  else
    oldstr = get(hand,'String');
    if isstruct(str)
      set(hand,'String',[oldstr;{str.message}]);
    else
      set(hand,'String',[oldstr;{str}]);
    end
    set(hand,'listboxtop',max(1,length(oldstr)-3))
    drawnow
  end
else
  if isempty(hand)
    if isstruct(str)
      DClab.disperror(str)
      rethrow(str)
    elseif isa(str,'MException')
        rethrow(str)
    else
      error(str)
    end
  else
    if isstruct(str)
      errordlg(str.message,'Internal Error','modal')
      DClab.disperror(str) %throw error, gui calls to functions that use this should have try catch to avoid full crash
      rethrow(str)
    else
      errordlg(str,'Internal Error','modal')
      DClab.disperror(str) %throw error, gui calls to functions that use this should have try catch to avoid full crash
      rethrow(str)
    end
  end
end
