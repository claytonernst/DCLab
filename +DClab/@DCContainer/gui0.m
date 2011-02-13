function gui(obj)
%how to call gui with object
%builds DC objects from DCContainer
%calls other gui6

ds = DClab.DCDataset(obj);

mop = ds.ModelAndObservationPair;
fp = ds.FreeParameter;

%targetmodels
targets = vertcat(mop.ResponseModel);

%open gui:
gui6(mop,fp,targets)

%ds = DClab.DCDataset(mop,fp);