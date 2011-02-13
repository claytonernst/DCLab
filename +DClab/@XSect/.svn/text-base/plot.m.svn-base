function handles = plot(obj,options)


if obj.ndims==2
    if nargin==1
        options.plotAllBool = 1;
        options.distinctFeasSetBool = 1;
    end
    handles = LOCALplot2(obj,options);
elseif obj.ndims==3
    if nargin==1
        options.plotAllBool = 1;
        options.distinctFeasSetBool = 1;
        options.smoothDataBool = 1;
        options.plotColorType = 1;
    end
    handles = LOCALplot3(obj,options);
end

function handles = LOCALplot2(obj,options)
% plot2dxsect1: method to plot feasible set cross-sections
%
% function [handles] = plot2dxsect1(twodxsectObj,options)
%
% INPUTS: 
%   obj: twoxsect object that supplies data to plot
%
%   options[optional]: Structure with the following fields
%     plotAllBool[default is 1]: Logical taking the following values
%       0: Plot just the combined feasible set from the experiments
%         in twoxsectObj.trg2use.
%       1: Plot the combined feasible set and individual feasible 
%         sets for each experiment listed in  twodxsectObj.trg2use
%
%     distinctFeasSetBool[default is 1]: Used only if
%       plotAllBool = 1.  If true(1), the individual feasible set 
%       plots will show the combined feasible set in a distinct color.
%
%=======finish previous========

% OUTPUTS
%  handles: Structure with the following fields:
%  FINISH ME            

% 1/9/04:  rpf, initial coding

ni = nargin;                           
no = nargout;

error( nargchk(1,2,ni) );
error( nargchk(0,1,no) );

% Initialize inputs
switch ni
 case 1
  options = [];
 otherwise
  %do nothing
end

% Initialize outputs
handles = [];

% Input error checking
if ~isempty(options) && ~isstruct(options)
  error('Inputs: invalid options')  
end

% Initialize options
try
  [plotAllBool,distinctFeasSetBool] = localInitOpts2(options);
catch
  error('Invalid options for plot2dxsect \n%s', lasterr)  
end

Nunits = length(obj.units);
resultSize = [obj.Npts obj.Npts];

% Reshape feasible set from vector to Npts-by-Npts
% matrix for plotting.
feasSet = cell(Nunits,1);
ofeasSet = cell(Nunits,1);
if distinctFeasSetBool && plotAllBool
  for i1 = 1:Nunits
    % zero out points that are in the combined outer feasible set.
    % These will be plotted in a different color
    feasSet{i1} = reshape(obj.indvFeasSet{i1}&~obj.combOFeasSet,resultSize);
    ofeasSet{i1} = reshape(obj.indvOFeasSet{i1}&~obj.combOFeasSet,resultSize);
  end
elseif plotAllBool
  for i1 = 1:Nunits
    feasSet{i1} = reshape(obj.indvFeasSet{i1},resultSize);
    ofeasSet{i1} = reshape(obj.indvOFeasSet{i1},resultSize);
  end
else
  %do nothing
end
% Reshape combined feasible set for plotting
combFeasSet = reshape(obj.combFeasSet,resultSize);
combOFeasSet = reshape(obj.combOFeasSet,resultSize);

% Generate grid locations
x1 = linspace(obj.paramBnds(1,1),obj.paramBnds(1,2),obj.Npts+1);
x1 = x1(2:end)-0.5*diff(x1(1:2));
x2 = linspace(obj.paramBnds(2,1),obj.paramBnds(2,2),obj.Npts+1);
x2 = x2(2:end)-0.5*diff(x2(1:2));

%=======Plot feasible set data==============================

%create a colormap for plots, elements of fullFeasSet vs. color:
%     1: white
%     2: blue
%     3: dark blue
%     4: red
%     5: dark red
cmat = [1 1 1;0 0 1;0 0 0.5;1 0 0;0.5 0 0];

%find 2nd subplot dimension
dim2 = ceil(Nunits/2);

if plotAllBool
  figure
  colormap(cmat)
  for i1 = 1:Nunits
    
    if distinctFeasSetBool
      %want to plot white outside the outer feas set
      %dark blue on the outer feas set / feas / OuterAll
      %blue on feas/ OuterAll
      %dark red on OuterAll/All
      %red on All
      tmp = 3*ofeasSet{i1}-feasSet{i1} + 5*combOFeasSet-combFeasSet;
    else
      tmp = 3*ofeasSet{i1}-feasSet{i1};
    end
  
    subplot(2,dim2,i1);
    image(x1,x2,tmp);
    axis xy; % get out of ij matrix plotting mode
    xlabel(['Param ' obj.params{1}]);
    ylabel(['Param ' obj.params{2}]);
    title(['DatasetUnit ' obj.units{i1}]);
    grid on
    axis square

  end
end

%plot combined feasible set
figure;
colormap(cmat)

% place 4 on combFeasSet and 5 on the extra from the outer feasible set

tmp = combOFeasSet*5-combFeasSet; 
image(x1,x2,tmp)

axis xy
xlabel(['Param ' obj.params{1}]);
ylabel(['Param ' obj.params{2}]);
title('Feasible Set of Combined Experiments');
grid on
axis square


%last line of plot2dfeasset
%==================== local functions =====================

function [plotAllBool,distinctFeasSetBool] = localInitOpts2(opt)

if ~isfield(opt,'plotAllBool');
  plotAllBool = true;
elseif ~isempty( intersect(opt.plotAllBool,[0 1]) )
  plotAllBool = opt.plotAllBool;
else
  warning(['Value in opt.plotAllBool was not 0 or 1, using ' ...
           'defaults.']);
  plotAllBool = true;
end

if ~isfield(opt,'distinctFeasSetBool');
  distinctFeasSetBool = true;
elseif ~isempty( intersect(opt.distinctFeasSetBool,[0 1]) )
  distinctFeasSetBool =  opt.distinctFeasSetBool;
else
  warning(['Value in opt.distinctFeasSetBool was not 0 or 1, using ' ...
           'defaults.']);
  distinctFeasSetBool = true;
end
   
%if ~isfield(opt,'plotColorType');
%  plotColorType = 1;
%elseif ~isempty( intersect(opt.plotColorType,[1 2 3 4 5]) )
%  plotColorType =  opt.plotColorType;
%else
%  warning(['Value in opt.plotColorType was not 1-5, using ' ...
%           'defaults.']);
%  plotColorType = 1;
%end

% last line of localInitializeOptions
%=====================================================

%to do: finish help, all relevant figure handles as output
% make sure we don't get keys and experiment numbers crossed up anywhere

% think about projections instead of just freezing all other
% variables to zero, also work on positioning plots in display.

function handles = LOCALplot3(obj,options)
% plot3dxsect1: method to plot feasible set cross-sections
%
% function [handles] = plot3dxsect1(threedxsectObj,options)
%
% INPUTS: 
%   obj: threedxsect object that supplies data to plot
%
%   options[optional]: Structure with the following fields
%     plotAllBool[default is 1]: Logical taking the following values
%       0: Plot just the combined feasible set from the experiments
%         in threedxsectObj.trg2use.
%       1: Plot the combined feasible set and individual feasible 
%         sets for each experiment listed in  threedxsectObj.trg2use
%
%     distinctFeasSetBool[default is 1]: Used only if
%       plotAllBool = 1.  If true(1), the individual feasible set 
%       plots will show the combined feasible set in a distinct color.
%
%     smoothDataBool[default is 1]: Indicates if we should attempt
%       to smooth the surface of the resulting plots. 
%       False(0) for no smoothing, (1) for smoothing.
%
%     plotColorType[default is 1]: A scalar that determines how the 
%       resulting plots are colored. Valid values are 1-5.
%       1:
%       2:
%=======finish previous========

% OUTPUTS
%  handles: Structure with the following fields:
%  FINISH ME

% 1/3/04:  rpf, initial coding

ni = nargin;                           
no = nargout;

error( nargchk(1,2,ni) );
error( nargchk(0,1,no) );

% Initialize inputs
switch ni
 case 1
  options = [];
 otherwise
  %do nothing
end

% Initialize outputs
handles = [];

% Input error checking
if ~isempty(options) && ~isstruct(options)
  error('Inputs: invalid options')  
end

% Initialize options
try
  [plotAllBool,distinctFeasSetBool,smoothDataBool,plotColorType] = localInitOpts3(options);
catch
  error('Invalid options for plot2dxsect \n%s', lasterr)  
end

Nunits = length(obj.units);
Npts = obj.Npts;
resultSize = Npts*ones(1,3);

% Reshape feasible set from vector to Npts-by-Npts-by-Npts
% matrix for plotting.
feasSet = cell(Nunits,1);
ofeasSet = cell(Nunits,1);
if distinctFeasSetBool && plotAllBool
  for i1 = 1:Nunits
    % zero out points that are in the combined outer feasible set.
    % These will be plotted in a different color
    feasSet{i1} = reshape(obj.indvFeasSet{i1}&~obj.combOFeasSet,resultSize);
    ofeasSet{i1} = reshape(obj.indvOFeasSet{i1}&~obj.combOFeasSet,resultSize);
  end
elseif plotAllBool
  for i1 = 1:Nunits
    feasSet{i1} = reshape(obj.indvFeasSet{i1},resultSize);
    ofeasSet{i1} = reshape(obj.indvOFeasSet{i1},resultSize);
  end
else
  %do nothing
end

% Reshape combined feasible set for plotting
combFeasSet = reshape(obj.combFeasSet,resultSize);
combOFeasSet = reshape(obj.combOFeasSet,resultSize);
  
if smoothDataBool
  combFeasSet = smooth3(combFeasSet,'box',3);
  combOFeasSet = smooth3(combOFeasSet,'box',3);
  if plotAllBool
    for i1 = 1:Nunits
      feasSet{i1} = smooth3(feasSet{i1},'box',3);
      ofeasSet{i1} = smooth3(ofeasSet{i1},'box',3);
    end
  end
end

% Generate grid locations
x1 = linspace(obj.paramBnds(1,1),obj.paramBnds(1,2),Npts);
%x1 = x1(2:end)-0.5*diff(x1(1:2));
x2 = linspace(obj.paramBnds(2,1),obj.paramBnds(2,2),Npts);
%x2 = x2(2:end)-0.5*diff(x2(1:2));
x3 = linspace(obj.paramBnds(3,1),obj.paramBnds(3,2),Npts);
%x3 = x3(2:end)-0.5*diff(x3(1:2));

%===========reshape delete me!!!!!!!!=========
%disp('Reshaping plots...')
%par2use = par2use([1 3 2]);
%for i1 = 1:Ntrgs
%  feasSet{i1} = permute(feasSet{i1},[3 2 1]);
%end
%%%%126 is right, 127 and 97 reversed
%%%%combFeasSet = permute(combFeasSet,[2 3 1]);
%combFeasSet = permute(combFeasSet,[3 2 1]);


%=======Plot feasible set data==============================
if ~isempty(intersect(plotColorType,[2 3 4 5]))
  [X Y Z] = meshgrid(x1,x2,x3);
end

if plotAllBool
  for i1 = 1:Nunits
    figure
    % plot feasible sets
    % Note on plotting. Each matrix of feasible sets contains
    % logical values (1 or 0). Isosurface will create a surface
    % of points whenever the matrix is equal to 0.5. This corresponds
    % to the surfaces on the interior of the unit cube because here the
    % logical values transition from 0 to 1.  Isocaps fills in the
    % the remaining surfaces where the feasible set is bounded by the
    % unit cube by enclosing any logical values greater than 0.5.
	
    %Setting camlight right or camlight left helps appearance in some
    %cases.  
    %To switch view angles use e.g. set(gca,'XDir','reverse').
    %This preserves axis labels better than manually rotating plot

    % plot surface
    p = patch( isosurface(x1,x2,x3,feasSet{i1},0.5) );
		
    %generate remaining faces
    pp = patch( isocaps(x1,x2,x3,feasSet{i1},0.5) );
		
    if distinctFeasSetBool
      %plot feasible set for combined experiments
      ppp = patch( isosurface(x1,x2,x3,combFeasSet,0.5) );
      pppp = patch( isocaps(x1,x2,x3,combFeasSet,0.5) );
    end
		
    switch plotColorType
     case 1
      set(p,'FaceColor','green','EdgeColor','black');
      set(pp,'FaceColor','green','EdgeColor','black');
      if distinctFeasSetBool
        set(ppp,'FaceColor','blue','EdgeColor','black');
        set(pppp,'FaceColor','blue','EdgeColor','black');
      end
        
     case 2
      [r g b] = meshgrid(Npts:-1:1,Npts:-1:1,1:Npts);
      isocolors(X,Y,Z,r/Npts,g/Npts,b/Npts,p);
      isocolors(X,Y,Z,r/Npts,g/Npts,b/Npts,pp);
      set(p,'FaceColor','interp','EdgeColor','none')
      set(pp,'FaceColor','interp','EdgeColor','none')
      if distinctFeasSetBool
        isocolors(X,Y,Z,r/Npts,g/Npts,b/Npts,ppp);
        isocolors(X,Y,Z,r/Npts,g/Npts,b/Npts,pppp);
        set(ppp,'FaceColor','white','EdgeColor','interp');
        set(pppp,'FaceColor','white','EdgeColor','interp');
      end
        
     case 3
      [r g b] = meshgrid(Npts:-1:1,1:Npts,1:Npts);
      isocolors(X,Y,Z,g/Npts,r/Npts,b/Npts,p);
      isocolors(X,Y,Z,g/Npts,r/Npts,b/Npts,pp);
      set(p,'FaceColor','interp','EdgeColor','none')
      set(pp,'FaceColor','interp','EdgeColor','none')
      if distinctFeasSetBool
        isocolors(X,Y,Z,g/Npts,r/Npts,b/Npts,ppp);
        isocolors(X,Y,Z,g/Npts,r/Npts,b/Npts,pppp);
        set(ppp,'FaceColor','white','EdgeColor','interp');
        set(pppp,'FaceColor','white','EdgeColor','interp');
      end
            
     case 4
      [r g b] = meshgrid(Npts:-1:1,Npts:-1:1,1:Npts);
      isocolors(X,Y,Z,g/Npts,r/Npts,b/Npts,p);
      isocolors(X,Y,Z,g/Npts,r/Npts,b/Npts,pp);
      set(p,'FaceColor','interp','EdgeColor','none')
      set(pp,'FaceColor','interp','EdgeColor','none')
      if distinctFeasSetBool
        isocolors(X,Y,Z,g/Npts,r/Npts,b/Npts,ppp);
        isocolors(X,Y,Z,g/Npts,r/Npts,b/Npts,pppp);
        set(ppp,'FaceColor','yellow','EdgeColor','none');
        set(pppp,'FaceColor','yellow','EdgeColor','none');
      end
                    
     case 5
        
      %create greyscale colormap
      %plot will transition from light to dark
      temp = (1:Npts)';
      temp2 = zeros(Npts,Npts,Npts);
      for i2 = 1:Npts
        for i3 = 1:Npts
          temp2(:,i2,i3) = temp+i2-1+i3-1;
        end
      end
      r = temp2/(Npts*3-2);
      r = r*0.75 + 0.1; %red, blue, and green equal for greyscale
      b = r;
      g = r;
              
      isocolors(X,Y,Z,g,r,b,p);
      isocolors(X,Y,Z,g,r,b,pp);
      % Matlab default for 'AmbientStrength' is 0.3
      % Valid range is [0, 1].  Add a bit more light.
      % Currently Ambient Strength is different for the feasible set vs.
      % the combined feasible set to highlight different sets.
      set(p,'FaceColor','interp','EdgeColor','none','AmbientStrength',0.5)
      set(pp,'FaceColor','interp','EdgeColor','none','AmbientStrength',0.5)
      if distinctFeasSetBool
        isocolors(X,Y,Z,g,r,b,ppp);
        isocolors(X,Y,Z,g,r,b,pppp);
        %set(ppp,'FaceColor',[0.9 0.9 0.9],'EdgeColor','none','AmbientStrength',0.5);
        %set(pppp,'FaceColor',[0.9 0.9 0.9],'EdgeColor','none','AmbientStrength',0.5);
        set(ppp,'FaceColor','interp','EdgeColor','none','AmbientStrength',0.3);
        set(pppp,'FaceColor','interp','EdgeColor','none', ...
                 'AmbientStrength',0.3);
      end
            
    end %switch
  
    %daspect([1 1 1])
    xlabel(['Param ' obj.params{1}]);
    ylabel(['Param ' obj.params{2}]);
    zlabel(['Param ' obj.params{3}]);
    title(['DatasetUnit ' obj.units{i1}]); 
    
    axis(reshape(obj.paramBnds',1,6));
    view(3); %axis tight
    camlight %left
    lighting phong%flat
    grid on
  end %for i1 = 1:Ntrg
end %if plotAllBool

% Plot combined feasible set
hall = figure;

hinner = subplot(1,2,2);
%plot surfaces
p = patch( isosurface(x1,x2,x3,combFeasSet,0.5,'noshare') );

%generate remaining faces
pp = patch( isocaps(x1,x2,x3,combFeasSet,0.5) );

houter = subplot(1,2,1);
%plot surfaces
p2 = patch( isosurface(x1,x2,x3,combOFeasSet,0.5,'noshare') );

%generate remaining faces
pp2 = patch( isocaps(x1,x2,x3,combOFeasSet,0.5) );

switch plotColorType
 case 1
  set([p pp p2 pp2],'FaceColor','blue','EdgeColor','black');
 case {2,4}
  [r g b] = meshgrid(Npts:-1:1,Npts:-1:1,1:Npts);
  isocolors(X,Y,Z,r/Npts,g/Npts,b/Npts,p);
  isocolors(X,Y,Z,r/Npts,g/Npts,b/Npts,pp);
  isocolors(X,Y,Z,r/Npts,g/Npts,b/Npts,p2);
  isocolors(X,Y,Z,r/Npts,g/Npts,b/Npts,pp2);
  
  set([p pp p2 pp2],'FaceColor','interp','EdgeColor','none');
%  set(pp,'FaceColor','interp','EdgeColor','none');
 case {3,5}
  [r g b] = meshgrid(Npts:-1:1,1:Npts,1:Npts);
  isocolors(X,Y,Z,g/Npts,r/Npts,b/Npts,p);
  isocolors(X,Y,Z,g/Npts,r/Npts,b/Npts,pp);
  isocolors(X,Y,Z,g/Npts,r/Npts,b/Npts,p2);
  isocolors(X,Y,Z,g/Npts,r/Npts,b/Npts,pp2);
  set([p pp p2 pp2],'FaceColor','interp','EdgeColor','none');
end

%daspect([1 1 1])
axes(hinner)
xlabel(['Param ' obj.params{1}]);
ylabel(['Param ' obj.params{2}]);
zlabel(['Param ' obj.params{3}]);
title([obj.innerOrTrue 'Feasible Set of Combined Experiments']);
axis(reshape(obj.paramBnds',1,6));
view(3); %axis tight
camlight left
lighting phong
grid on

axes(houter)
xlabel(['Param ' obj.params{1}]);
ylabel(['Param ' obj.params{2}]);
zlabel(['Param ' obj.params{3}]);
title('outerFeasible Set of Combined Experiments');
axis(reshape(obj.paramBnds',1,6));
view(3); %axis tight
camlight left
lighting phong
grid on

p = get(hall,'Position');
p(1) = 300;
p(3) = 900;
set(hall,'Position',p)

%lighting flat





%last line of plot3dfeasset.m====================================


% for i1 = 1:Ntrgs
%   %fullFeasSet{i1} = feasSet{i1}*-35 + combFeasSet*20 +40;
%   fullFeasSet{i1} = feasSet{i1}+combFeasSet*2 +1;
% end

% figure
% colormap([1 1 1;0 0 1;1 0 0])
% clear feasSet combFeasSet
% asdf = fullFeasSet{1};
% asdf = permute(asdf,[3 1 2]);
% image(x1,x2,asdf(:,:,41))
% axis xy



%clear all

%==================== local functions =====================

function [plotAllBool,distinctFeasSetBool,smoothDataBool,plotColorType] ...
    = localInitOpts3(opt)

if ~isfield(opt,'plotAllBool');
  plotAllBool = true;
elseif ~isempty( intersect(opt.plotAllBool,[0 1]) )
  plotAllBool = opt.plotAllBool;
else
  warning(['Value in opt.plotAllBool was not 0 or 1, using ' ...
           'defaults.']);
  plotAllBool = true;
end

if ~isfield(opt,'distinctFeasSetBool');
  distinctFeasSetBool = true;
elseif ~isempty( intersect(opt.distinctFeasSetBool,[0 1]) )
  distinctFeasSetBool =  opt.distinctFeasSetBool;
else
  warning(['Value in opt.distinctFeasSetBool was not 0 or 1, using ' ...
           'defaults.']);
  distinctFeasSetBool = true;
end
   
if ~isfield(opt,'smoothDataBool');
  smoothDataBool = true;
elseif ~isempty( intersect(opt.smoothDataBool,[0 1]) )
  smoothDataBool =  opt.smoothDataBool;
else
  warning(['Value in opt.smoothDataBool was not 0 or 1, using ' ...
           'defaults.']);
  smoothDataBool = true;
end

if ~isfield(opt,'plotColorType');
  plotColorType = 1;
elseif ~isempty( intersect(opt.plotColorType,[1 2 3 4 5]) )
  plotColorType =  opt.plotColorType;
else
  warning(['Value in opt.plotColorType was not 1-5, using ' ...
           'defaults.']);
  plotColorType = 1;
end

% last line of localInitializeOptions
%=====================================================

%to do: finish help, all relevant figure handles as output
% make sure we don't get keys and experiment numbers crossed up anywhere

%most of the ploting options don't work because X, Y, and Z don't exist

