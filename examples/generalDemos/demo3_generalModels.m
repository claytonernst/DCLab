% DEMO 3  Creating a Dataset whose models have general
% (non-algebraic) dependence on the parameters


% This is a purely artificial example constructed for demonstration
% purposes only. We assemble a dataset based on measurements of two
% quantities Y1 and Y2. Models for each of these quantities are
% functions of three parameters X1, X2, and X3. This demonstration
% differs from demo1.m only in the construction of the model
% assertions. Users who have read demo1 may read this rather quickly,
% focusing only on the construction of the model assertions.
%
% Corresponding to each measured quantity (observable, or subject of
% measurement), we form a 'ModelAndObservationPair', which communicates to
% the software what is known about the quantity -- as expressed both by
% modeling and by recorded observation. The instrument by which a
% model is specified is a 'ResponseModel', while observational
% information is input as an 'ResponseObservation'. In a similar vein,
% known information about possible values of model parameters is
% supplied as 'FreeParameters'. All of this information is
% consolidated into a 'Dataset', which is the core object upon which
% analysis is performed.
%
% This demonstration covers the construction of a Dataset for two
% observations whose observables are modeled by linear functions of
% three parameters.

% 1. Construction of the ModelAndObservationPair that concerns the
% observable Y1. 
%
%    Suppose with the first measurement we observe Y1 to have a value
%    of 5.0 with a 10% level of uncertainty. The following command
%    produces an ResponseObservation object with this information.

RO1 = ResponseObservation(5,0.5)

%    We will use a model for Y1 that takes the form 
%
%                       Y1 = z1(10)
%
%    where z1,z2 satisfy the ordinary differential equation
%    
%    z1_dot = x1*z1+sin(z2), z2_dot = cos(x2*z1) - x3, z1(0) = z2(0) = 0
%
%    where x1,x2,x3 are possible values for the uncertain parameters
%    X1,X2,X3 resp. A useful description of this model should (at
%    minimum) provide this expression and a correspondence between
%    each variable (e.g., x1), and the parameter it represents (e.g.,
%    X1). The Data Collaboration toolbox additionally requires that
%    the range of each variable be specified. For algebraic models
%    (linear or quadratic) this is justified/motivated on physical
%    grounds. Typically such a model arises as an approximation of
%    some process that is only locally (linear/quadratic), so the
%    model should be accompanied by a description of the variable
%    ranges for which it is thought to be valid. When the model is of
%    more general form (such as the present case) supplying a domain
%    of validity (while perhaps relevant) is less natural. Provide
%    this when it makes sense, and use the range [-inf inf]
%    otherwise. In any analyses conducted by the Data Collaboration
%    toolbox software, the model use will be restricted to the domain
%    specified by the parameter assertions.
%
%    We create a model assertion for Y1 with the command
%
RM1 = ResponseModel(@model1)

%    where model1.m is function mfile in this directory. The function
%    conforms to a specific syntax that is elaborated upon on the html
%    documentation on the ModelAssertion constructor. A template
%    functions that serve as the basis for model assertions is
%    provided as ../DClab/filetemplates/dcModel.m. Although you can
%    base ResponseModel on a local functions, we recommend that in
%    your work you use a separate file (as we do in this
%    demonstration) for this purpose. This is so MATLAB can find and
%    evaluate such functions in the event that you load saved sessions
%    of your work.
%
%    Now we create the ModelAndObservationPair. The toolbox references these
%    objects by a unique name (label, identifier) that must be
%    supplied along with the ResponseObservation and ResponseModel
%    already created. We'll use the name 'Y1'. 

mop1 = ModelAndObservationPair(RO1,RM1,'Y1')

% 2. Construction of the ModelAndObservationPair that concerns the
% observable Y2. 
%
%    This proceeds much like the previous case. We will use
%    make-believe data, uncertainty, and model. The only difference of
%    consequence is that is that for illustrative purposes, we assume
%    the model is a function of only two parameters: X1 and X3.
%
%    We assume we measure Y2 to be 11.0 with 5% uncertainty.
%    Additionally we suppose the model for Y2 is the mean value of 
%    z1 over a 5 second time window, where z1 satisfies the
%    ordinary differential equation
%
%         z1_dot = x1*exp(-x2*z1), z1(0) = 0
%
%    Mimicking the previous steps, we create the appropriate 
%    DatasetUnit.

%create the ResponseObservation object
RO2 = ResponseObservation(11,11*0.05) %5 percent uncertainty

%create the model assertion object derived from the file model2.m
%in this directory.
RM2 = ResponseModel(@model2);

mop2 = ModelAndObservationPair(RO2,RM2,'Y2')

% 3. Construction of FreeParameters for X1, X2, and X3. 
%
%    In this segment of the demonstration we make explicit assumptions
%    about the range of each parameter (i.e., its uncertainty). In
%    practice, such information may be obtained through direct
%    measurements on the parameters, or by consideration of prior
%    information regarding the system of study. 
%
%    We construct an assertion for each parameter by supplying the
%    parameter name and a range specifying the uncertainty in its
%    value. For example, the code below asserts that the value of X1
%    lies between 0 and 7.

FP1 = FreeParameter('X1',3,[-2 4])
FP2 = FreeParameter('X2',0,[-5 5])
FP3 = FreeParameter('X3',5.5,[-4.5 4.5])

% 4. Construction of the Dataset
%
%    A dataset consists of a collection of dataset units, together
%    with parameter assertions for all relevant parameters. The object
%    is constructed by concatenating together the dataset units and
%    parameters assertions, and passing these multidimensional
%    entities to the constructor. Notice that we use a vertical
%    concatenation (semicolon delimiter).

MOPArray = [mop1; mop2];
FPArray = [FP1; FP2; FP3];

Dset = DCDataset(MOPArray,FPArray)

% 5. This concludes the demonstration. For information on types of
%    analysis that may be performed on this dataset, see demo4.m

