classdef Guisession < DClab.DCObject
    % function guiobj = Guisession(handles)
    %
    % inputs: handles in the guihandles object from the main gui (e.g., gui6)
    %
    % Specific to Ryan's GUI....maybe obsolete...could use DCContainer?

    properties
        options = DClab.DCOptions;
        DatasetUnit_name = {};
        DatasetUnit_value = {};
        DatasetUnitInfo = struct('def_d',{},'def_uncVect',{},'isdef',{},'curwt',{});
        Model_name = {};
        Model_value = {};
        ModP_name = {'default'};
        ModP_value = {DClab.FreeParameter};
        ConsistencyTest_name = {};
        ConsistencyTest_value = {};
        ResponsePrediction_name = {};
        ResponsePrediction_value = {};
        ParameterOptimization_name = {};
        ParameterOptimization_value = {};
        UnitArray = DClab.ModelAndObservationPair;
        FreeParameter = DClab.FreeParameter;
        RMcell = {};
    end

    methods

        function guiobj = Guisession(handles)

            if nargin >0
                guiobj.options = handles.options;
                guiobj.DatasetUnit_name = get(handles.DsetList(1),'String');
                guiobj.DatasetUnit_value = handles.DatasetUnit;
                guiobj.DatasetUnitInfo = handles.DatasetUnitInfo;
                guiobj.Model_name = get(handles.ModelList(1),'String');
                guiobj.Model_value = handles.Model;
                guiobj.ModP_name = get(handles.PDomPU,'String');
                guiobj.ModP_value = handles.ModP;
                guiobj.ConsistencyTest_name = get(handles.ConsistResultsLB,'String');
                guiobj.ConsistencyTest_value = handles.ConsistTest;
                guiobj.ResponsePrediction_name = get(handles.PredResultsLB,'String');
                guiobj.ResponsePrediction_value = handles.Prediction;
                guiobj.ParameterOptimization_name = get(handles.OptimResultsLB,'String');
                guiobj.ParameterOptimization_value = handles.ParameterOptimization;

                guiobj.UnitArray = handles.UnitArray;
                guiobj.FreeParameter = handles.ParamAss;
                guiobj.RMcell = handles.MAcell;
            end

        end
        
        function [list,sz] = displayProps(obj)
            list = [];
            sz = '';
        end
    end %public methods
end %classdef
