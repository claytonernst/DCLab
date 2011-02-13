function out = stopGui
% StopGui  creates a stop button for functions
%
%  Call StopGui somewhere in your code.  It will return the
%  function handle to figure, for your use.  If there is no other figures
%  open at the time, its figure handle will be 1.
%
%  At the point where you'd like to see if you function should stop, it
%  should make the following call:
%               getappdata(figHandle,'stopIt')
%  If this returns true, then the button has been pressed, otherwise, it
%  hasn't.  You may do with information as you please.
%
%  You cannot close the stopGui figure until you've presed the Stop!
%  button.  If you want to close the figure, before the user clicked stop
%  (i.e. you don't want them to stop anymore), use:
%               delete(figHandle);
%

ss = get(0,'screensize');

fh = figure('position', [0.5*ss(3:4) 105 40], ...
            'menubar', 'none', ...
            'name', 'StopGUI', ...
            'NumberTitle', 'off', ...
            'Resize', 'off', ...
            'CloseRequestFcn', 'error(''Not until you''''ve stopped...'')');
        
bh = uicontrol('string', 'Stop!', ...
               'position', [27 5 51 30], ...
               'CallBack', {@LOCAL_stopIt,fh});

setappdata(fh,'stopIt',false);

if nargout==1
    out = fh;
end

function LOCAL_stopIt(obj,eventH,fh)
setappdata(fh,'stopIt',true);
set(fh,'CloseRequestFcn',['delete(' num2str(fh) ');']);
