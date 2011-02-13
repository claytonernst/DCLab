function PD = subdivideDomain(PD,varargin)
%SUBDIVIDEDOMAIN subdivides the domain of the PolyDataset
%
%   This method simply calls the identically named method of the
%   PiecewiseSurrogateModelTree class. Perhaps PolyDataset should just
%   multiply inherit...
%
%   PD = SUBDIVIDEDOMAIN(PD,LEAF,PARAMNAME,LOC) subdivides the subdivision
%   LEAF of the subdivided domain of the PolyDataset PD in the dimension
%   PARAMNAME at the location LOC. This causes new DCSurfaces for the
%   response models of each MOPair in PD to be constructed over the new
%   subdivisions, effectively refining the surrogate models for these
%   response models. The display and fitting options used to construct the
%   new DCSurfaces will be inherited from those of LEAF.
%
%   PD = SUBDIVIDEDOMAIN(PD,LEAF,PARAMNAME,LOC,INHERIT) allows you to
%   specify if the DCSurfaces defined over the resulting subdivisions should
%   be recomputed or if they should be inherited from their parents. See
%   input description below.
%
%   PD = SUBDIVIDEDOMAIN(PD,LEAF,PARAMNAME,LOC,INHERIT,OPT) allows you to
%   supply a DCOptions object OPT that will dictate the display and fitting
%   options used to construct the new DCSurfaces.
%
%   Inputs
%      PD: A PolyDataset object.
%      LEAF: An integer that specifies which subdivision will be
%        subdivided. It must be a member of
%        PD.PiecewiseSurrogateModelTree.treeLeafIndices 
%      PARAMNAME: This is a char array specifying the dimension that the
%         partitioning plane passes though. PARAMNAME must be the name of a
%         parameter in PD.
%      LOC: This scalar specifies the location of the partitioning plane
%         that passes though the dimension indicated by PARAMNAME.
%      INHERIT: A column vector consisting of ones and zeros. Its length
%         must be equal to the number of DCSurfaces that are defined on
%         LEAF. I.e., equal to get(PD.Tree,'nSurfaces',LEAF). All
%         DCSurfaces that exist on LEAF and correspond to components of
%         INHERIT that equal 1 will be inherited by the new subdivisions.
%         For components of INHERIT that equal 0, this function determines
%         which MOPairs they correspond to and then creates DCSurfaces over
%         the new subdivisions for the response models of these pairs. This
%         process might involve the construction of multiple DCSurfaces for
%         any given response model, and their I/O transformations may differ
%         from those of the previous DCSurfaces that existed on LEAF for a
%         given MOPair and that were not inherited. The sum total of these
%         events mean that a given MOPair may have several DCSurfaces over
%         a new subdivision, some of which were simply inherited, and some
%         of which were newly created. In the rare case where I/O
%         transformations for which a DCSurface was inherited are also used
%         by a newly created DCSurface, the inherited surface will be
%         thrown away and only the new surface will remain. If INHERIT is a
%         scalar 0, it will be interpreted as zeros(nSurf,1).
%      OPT: A DCOptions object. 
%
%   Ouputs
%      PD: The resulting PolyDataset object.
%
% See also PolyDataset, PiecewiseSurrogateModelTree
%
%   Function reference page in Help browser
%      <a href="matlab:web(fullfile(fileparts(which('DCsetup')), ...
%      'docs', 'html', 'functions', 'PDsplit.html'), ...
%      '-helpbrowser')">PolyDataset/split</a>


m = PD.nPairs;
RMCell = cell(m,1);
for i1 = 1:m
  RMCell{i1} = PD.ModelAndObservationPair(i1).ResponseModel;
end

PD.PiecewiseSurrogateModelTree = subdivideDomain(PD.PiecewiseSurrogateModelTree,RMCell,varargin{:});
