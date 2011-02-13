% A quick script to test the installation of the Data Collaboration
% software. This in not intended to thoroughly test the software
% itself. This script simply employs the basic functionality to
% ensure that needed optimization routines are available.

disp('Testing software installation...')

% Check the matlab version. Has to be at least 7.1 (SP3) in order for Data Collaborationto run
[y m] = datevec(version('-date'));
if y < 2008 && m < 2
    disp('Installation failed!')
    disp('This Data Collaboration package is not compatible with Matlab releases prior to 7.6.0 (R2008a)')
    return
end

% Test if sedumi is installed
try
  whichcpx(struct([])); %call a random sedumi function that will be mexed if it is installed
catch
  disp(' ')
  disp('Installation failed!')
  disp('Launch install_sedumi from $/DClabV2/SeDuMi_1_1') 
  disp('  This will compile several mex files that are employed by the SeDuMi optimization package.')
  disp('  MATLAB''s mex function must be properly configured to work with your local C compiler')
  disp('  in order to successfully install SeDuMi.')
  disp('  If the call to install_sedumi fails, you have more work to do...')
  disp('  Typing mex -setup or help mex at the command window should get you started.')
  disp(' ')
  disp('Due to a MATLAB limitation, Unix users with very up-to-date installations')
  disp('  may need to install an older version of gcc to complete the operation.')
  disp(' ')
  disp('After a successful call to the install_sedumi script, run dctest.m again.')
  %  disp(['You can ignore any messages from SeDuMi about modifying ' ...
%        'your path'])
  return
end

% Test if the new version of sedumi is required because the MATLAB version is >= 2006B
%my_fprintf is a file included in SeDuMi 1.1R3 to avoid fprintf problem
if y >= 2006 && m >= 8 && ~exist('my_fprintf','file') 
    disp('Installation failed!')
    disp('The version of SeDuMi on your path is not 1.1R3 or later!')
    disp(['  Typically this occurs if DCsetup detected an existing ' ...
          'installation of'])
    disp(['  SeDuMi on your path, and thus did not add the version ' ...
          'included with this'])
    disp('  toolbox to your MATLAB path.')
    disp('Please modify your MATLAB startup procedure or use RMPATH to remedy this, and run dctest.m again.')
    return
end

% Create domain of quadratic model.
dom(1).name = 'p1';
dom(1).range = [-10 -5];
dom(2,1).name = 'p2';
dom(2).range = [3 5];

% Create assertion objects
MA = ResponseModel([1 3 -7; 3 5 -1;-7 -1 0.5],dom);
EA = ResponseObservation(150,30);
PA = [FreeParameter('p1',-6.5,[-1.5 1.5]); FreeParameter('p2',3.5,[-0.5 0.5])];

% Create dataset and analysis options
Unit = ModelAndObservationPair(EA,MA,'test1');
Dset = DCDataset(Unit,PA);
opt = DCOptions('display','off','maxBranchBoundIter',2);

% Test consistency call
ConsistencyTest(Dset,opt);

% Test prediction call
Mnot = ResponseModel(eye(3),dom);
ResponsePrediction(Mnot,Dset,opt);

disp('Installation successful!')
