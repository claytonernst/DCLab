function disperror(err)
%DISPERROR Display latest error message and where the error occured.
%   DISPERROR displays the message of the last error and also where this
%   error occured. This is helpful when the error was captured by a try-
%   catch block and then rethrown. Matlab will show the error message, but
%   unfortunately not where the error occured. DISPERROR displays the error
%   message and the error stack, complete with hyperlinks.
%
%   Another use for DISPERROR is with callbacks containing try-catch
%   block(s). In this case, however, the call to disperror must be in the
%   catch part of a try-catch block, such as:
%
%     set(gcf, 'ButtonDownFcn', 'try;fun;catch disperror;end');
%   or
%     set(gcf, 'ButtonDownFcn', 'try;fun;catch disperror;rethrow(lasterror);end');
%
%   The latter usage will cause the error message to be repeated, but will also
%   give information on which callback that was executed.
%
%   Example:
%
%      function testdisperror
%         try
%            a=[1 2 3];
%            b=a(4711);  % Index error
%         catch
%            % possibly do some housekeeping here
%            rethrow(lasterror);
%            end
%
%      >> testdisperror
%      ??? Attempted to access a(4711); index out of bounds because numel(a)=3.
%
%      >> disperror
%      ??? Attempted to access a(4711); index out of bounds because numel(a)=3.
%      Error in ==> testdisperror at 4
%                   ------------------  % <- Hyperlink
%
%   =======================
%   Version: 1.0 2006-02-27
%   Author: Jerker W?gberg, More Research, SWEDEN
%   email: char(hex2dec(reshape('6A65726B65722E77616762657267406D6F72652E7365',2,[])')') 

%rpf added err as an input and removed err = lasterror
%     err = lasterror;

     % See if it is a syntax error
      match=regexp(err.message ...
                   , {'(?<=File: *)([^ ].*)(?= *Line:)','(?<=Line: *)\d+'} ...
                   , 'match');
      sp=1;
      mf=[];
      if any(cellfun('isempty',match))
         % Not syntax error, mimic Matlabs first error line
         disp(['??? ' err.message]);
         if length(err.stack)>0
           [mf,fn,line]=editem(err.stack(1));
           sp=2;
         end
      else
        % Syntax error, first line of message contains info for hyperlink
        mf=char(match{1});
        line=str2num(char(match{2}));        % Show the actual error
        disp(['??? ' char(regexp(err.message,'(?<=\n)(.*)','match'))]);
        [fn,fn]=fileparts(mf);
      end
      if ~isempty(mf)
        % Show the error stack
        disp(sprintf('Error in ==> <a href="matlab:opentoline(''%s'',%d)">%s at %d</a>',mf,line,fn,line));
        for i=sp:length(err.stack)
          [mf,fn,line]=editem(err.stack(i));
          disp(sprintf('  In <a href="matlab:opentoline(''%s'',%d)">%s at %d</a>',mf,line,fn,line));
        end
      end 

function [mf,fn,line]=editem(x)
      % Build the function name
      [p,mnam]=fileparts(x.file);
      fn=x.name;
      if ~strcmpi(regexp(fn,'^[^/]+', 'match'),mnam)
        fn = [mnam '>' fn];
      end
      mf=x.file;
      line=x.line;
