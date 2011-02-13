classdef DCContainer < DClab.DCObject
    % DCContainer object
    %
    %   obj = DCContainer;  Creates an empty DCContainer object.  Properties
    %   can then be set using set(obj,propertyname,propertyvalue);
    %
    %   obj = DCContainer(propname1,propvalue1,propname2,propvalue2,...);
    %   Creates a DCContainer object with the specified property name/value
    %   pairs.
    %
    %   Properties:
    %       parameterKey: A numParams-by-1 cell array of strings denoting the
    %           parameter keys
    %       parameterID: A numParams-by-1 cell array of strings denoting the
    %           parameter IDs
    %       parameterRange: A numParams-by-1 cell array.  Each entry contains a
    %           1-by-2 double array that are the lower and upper bounds on the
    %           parameter value
    %       paramInitialValue: A numParams-by-1 cell array.  Each entry contains
    %           a 1-by-1 double array that is the initial parameter value.
    %       paramTrans: A numParams-by-1 cell array of strings denoting
    %           the transformation to return actual variables from parameters.
    %       centerSpan: A numParams-by-1 cell array. Each entry contains a
    %           1-by-2 double array that are the center and span in the
    %           transformation expression.
    %       optVarUnits: A numParams-by-1 cell array of strings denoting
    %           the units of the actual variables.
    %       parameterLinks: A numParams-by-1 cell array of strings denoting
    %           the relevant primeIDs for parameters.
    %       modelKey: A numModels-by-1 cell array of strings denoting the model
    %           keys
    %       modelID: A numModels-by-1 cell array of strings denoting the model
    %           IDs
    %       modelCoeffs: A numModels-by-1 cell array.  If the i-th model
    %           depends upon m parameters, then the i-th cell contains an
    %           (m+1)-by-(m+1) double array representing the quadratic
    %           coefficients.  If the matrix is Q, then the model is
    %           y=[1;x]'*Q*[1;x] for a m-by-1 input x.
    %       modelParamIDs: A numModels-by-1 cell array.  If the i-th model
    %           depends upon m parameters, then the i-th cell contains an
    %           m-by-1 cell array of strings, denoting the parameter IDs of the
    %           parameters the cooresponding model depends on.  These do not
    %           need to be in the same order as listed in the parameterID
    %           field, but they do need to be in the same order as the
    %           models defined the modelCoeff field expects them.
    %       targetKey: A numModels-by-1 cell array (number of models should
    %           be the same as the number of targets) of strings denoting the
    %           target keys.
    %       targetID: A numModels-by-1 cell array of strings denoting the
    %           target IDs.
    %       targetLabel: A numModels-by-1 cell array of strings describing the
    %           target.  Will be used as an axis label.  E.g. 'log of target
    %           range'.
    %       uncRange: A numModels-by-1 cell array.  Each entry contains a
    %           1-by-2 double denoting the range of values predicted by the
    %           experiment (e.g. if d is the measured value and there is a
    %           plus/minus uncertainty u, then the range would be [d-u, d+u])
    %       targetValue: A numParams-by-1 cell array.  Each entry contains
    %           a 1-by-1 double array denoting the measured value d.
    %       targetTrans: A numParams-by-1 cell array of strings denoting
    %           the transformation to return actual values from targets.
    %       targetUnits: A numParams-by-1 cell array of strings denoting
    %           the units of the actual targets.
    %       targetLinks: A numParams-by-1 cell array of strings denoting
    %           the relevant primeIDs for targets.
    %       trialModelID: A numModels-by-1 cell array of strings denoting the
    %           trial model ID
    %       predictionModelKey: A numPredModels-by-1 cell array of
    %           strings denoting the prediction model keys
    %       predictionModelID: A numPredModels-by-1 cell array of
    %           strings denoting the prediction model IDs
    %       predictionModelCoeffs: A numPredModels-by-1 cell array.
    %           If the i-th model depends upon m parameters, then
    %           the i-th cell contains an (m+1)-by-(m+1) double
    %           array representing the quadratic coefficients.  If
    %           the matrix is Q, then the model is y=[1;x]'*Q*[1;x]
    %           for a m-by-1 input x. 
    %       predictionModelParamIDs: A numPredModels-by-1 cell
    %           array.  If the i-th model depends upon m
    %           parameters, then the i-th cell contains an m-by-1
    %           cell array of strings, denoting the parameter IDs
    %           of the parameters the cooresponding model depends
    %           on.  These do not need to be in the same order as
    %           listed in the parameterID field, but they do need
    %           to be in the same order as the models defined the
    %           modelCoeff field expects them. 
    %       callback2show: A 1-by-1 function handle. The associated function takes a
    %           PrIMe ID as input and opens the associated viewer.
    %           e.g. callback2show('a00000001') %targets/dataAttributes
    %                callback2show('m00000003') %trial model
    %                callback2show('v00000001') %parameters or optVariables
    %       callback2show2: A 1-by-1 function handle. The associated function
    %           takes two PrIMe IDs as input and opens the associated viewer.
    %           e.g. callback2show('r00012186','rk00000005') %rate constants
    %                callback2show('m00000003','sm00000001') %surrogate models
    %
    %       The order of each field should be the same: the parameterKey,
    %       parameterID, and parameterRange fields should all be in the same
    %       order; the modelKey, modelID, modelCoeffs, and modelParamIDs should
    %       all be in the same order.
    %

    properties

        parameterKey = {};
        parameterID = {};
        parameterRange = {};
        paramInitialValue = {};
        paramTrans = {};
        centerSpan = {};
        optVarUnits = {};
        parameterLinks = {};
        modelKey = {};
        modelID = {};
        modelCoeffs = {};
        modelParamIDs = {};
        targetKey = {};
        targetID = {};
        targetLabel = {};
        uncRange = {};
        targetValue = {};
        targetTrans = {};
        targetUnits = {};
        targetLinks = {};
        trialModelID = {};
        predictionModelKey = {};
        predictionModelID = {};
        predictionModelCoeffs = {};
        predictionModelParamIDs = {};
        callback2show;
        callback2show2;
    end

    methods
        function obj = DCContainer(varargin)
            ni = nargin;
            if ni>0
                if mod(ni,2)~=0
                    error('Incorrect number of inputs')
                else
                    for i1=1:(ni/2)
                        obj.(varargin{2*i1-1}) = varargin{2*i1};
                    end
                end
            end
        end
        
        function [list,sz] = displayProps(obj)
            list = [];
            sz = '';
        end
    end %public methods
    
end %classdef