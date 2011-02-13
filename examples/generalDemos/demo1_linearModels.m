% DEMO 1  Creating a Dataset whose models depend linearly on the
% parameters


% This is a purely artificial example constructed for demonstration
% purposes only. We assemble a dataset based on measurements of two
% quantities Y1 and Y2. Models for each of these quantities are linear
% functions of three parameters X1, X2, and X3.
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

% 1. Construction of the dataset unit that concerns the observable Y1.
%
%    Suppose with the first measurement we observe Y1 to have a value
%    of 5.0 with a 10% level of uncertainty. The following command
%    produces an ResponseObservation object with this information.

RO1 = ResponseObservation(5,0.5)

%    We will use a model for Y1 that takes the form 
%
%                     2 + 3x1 + 5x2 - 2.5x3,
%
%    where x1,x2,x3 are possible values for the uncertain parameters
%    X1,X2,X3 resp. A useful description of this model should (at
%    minimum) provide this algebraic expression and a correspondence
%    between each variable (e.g., x1), and the parameter it represents
%    (e.g., X1). The Data Collaboration toolbox additionally requires
%    that the range of each variable be specified. This is
%    justified/motivated on physical grounds. Typically a linear model
%    arises as an approximation (i.e., linearization) of some process
%    that is only locally linear, so the linear model should be
%    accompanied by a description of the variable ranges for which it
%    is thought to be valid. We create a model assertion for Y1 by
%    first specifying the domain of the linear model. This will take
%    the form of a 3x1 structure array with fields .name and
%    .range. We suppose the linear model was created by assuming each
%    parameter (X1,X2,or X3) takes values between +/- 10. If the
%    linear model is thought to be value for all values of a
%    particular variable, setting the range field to [-inf, inf] would
%    be appropriate.

%initialize a 3x1 structure array with appropriate fields
model_domain = struct('name',cell(3,1),'range',cell(3,1));
%fill out the structure array
model_domain(1).name = 'X1';
model_domain(1).range = [-10 10];
model_domain(2).name = 'X2';
model_domain(2).range = [-10 10];
model_domain(3).name = 'X3';
model_domain(3).range = [-10 10];

%vector representing the algebraic expression (note the required use of a
%column vector)
algForm = [2; 3; 5; -2.5];

%create the assertion object
RM1 = ResponseModel(algForm,model_domain)

%    Notice that the correspondence between variables in the algebraic
%    expression and the parameters is established by order--for
%    example the second component of algForm multiplies the parameter
%    described by the first component of the structure array
%    model_domain, the third multiplies the second, etc.

%    Now we create the DatasetUnit. The toolbox references these
%    objects by a unique name (label, identifier) that must be
%    supplied along with the ResponseObservation and ResponseModel
%    already created. We'll use the name 'Y1'. 

mop1 = ModelAndObservationPair(RO1,RM1,'Y1')

% 2. Construction of the dataset unit that concerns the observable Y2.
%
%    This proceeds much like the previous case. We will use
%    make-believe data, uncertainty, and model. The only difference of
%    consequence is that is that for illustrative purposes, we assume
%    the model is a function of only two parameters: X1 and X3.
%
%    We assume we measure Y2 to be 11.0 with 5% uncertainty.
%    Additionally we suppose the model for Y2 is 2 + 4.1x1 + 0.9x3. 
%    Mimicking the previous steps, we create the appropriate 
%    DatasetUnit.

RO2 = ResponseObservation(11,11*0.05) %5 percent uncertainty

%initialize 2x1 structure
model_domain = struct('name',cell(2,1),'range',cell(2,1));
model_domain(1).name = 'X1';
model_domain(1).range = [-10 10];
model_domain(2).name = 'X3'; %notice the 2nd parameter is X3
model_domain(2).range = [-10 10];

%vector representing the algebraic expression (note the required use 
% of a column vector)
algForm = [2; 4.1; 0.9];

%create the assertion object
RM2 = ResponseModel(algForm,model_domain)

mop2 = ModelAndObservationPair(RO2,RM2,'Y2')

% 3. Construction of FreeParameters for X1, X2, and X3. 
%
%    In the previous steps, we assumed that the linear models were
%    valid for parameters values between +/- 10. In this segment of
%    the demonstration we make explicit assumptions about the range of
%    each parameter (i.e., its uncertainty). In practice, such
%    information may be obtained through direct measurements on the
%    parameters, or by consideration of prior information regarding
%    the system of study. In this demonstration, the ranges that
%    specify the uncertainty in the parameters will be smaller than
%    the ranges over which the linear models were assumed valid. This
%    highlights the possible distinction that exists between the act
%    of hypothesizing a model, and that of estimating the uncertainty
%    in a model parameter (these could even be the work of different
%    research groups). In general, the presence of such differences
%    will depend on the application. We note however, that if the
%    ranges of uncertainty in a parameters (as specified through
%    parameter assertions in this portion of the demo) are larger than
%    those used when constructing the model assertions, an error will
%    be thrown when the dataset is constructed. We seek to prevent a
%    model from being used in conditions (parameter values) where it
%    has not been verified, or worse is known to be invalid. It is to avoid
%    this situation that we require the domain of an asserted model to
%    be made explicit.
%
%    We construct an assertion for each parameter by supplying the
%    parameter name and a range specifying the uncertainty in its
%    value. For example, the code below asserts that the value of X1
%    has a nominal of 3 with uncertainty of -2 below and 4 above (i.e. it
%    lies between 1 and 7).

FP1 = FreeParameter('X1',3,[-2 4])
FP2 = FreeParameter('X2',0,[-5 5])
FP3 = FreeParameter('X3',5.5,[-4.5 4.5])

% 4. Construction of the Dataset
%
%    A dataset consists of a collection of dataset units, together
%    with FreeParameters for all relevant parameters. The object
%    is constructed by concatenating together the ModelAndObservationPairs
%    and FreeParameters, and passing these multidimensional
%    entities to the constructor. Notice that we use a vertical
%    concatenation (semicolon delimiter).

MOPArray = [mop1; mop2];
FPArray = [FP1; FP2; FP3];

Dset = DCDataset(MOPArray,FPArray)

% 5. This concludes the demonstration. For information on types of
%    analysis that may be performed on this dataset, see demo4.m

