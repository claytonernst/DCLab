%% Info for 76 polys
% This file illustrates how the polynomial inequalities
% are accessed, and evaluated.  Andrew Packard, Fall 2010.

%% Load the data
% |D| is a 76-by-1 struct, and |Xfeas| is 102-by-1 |double|
load CAGdata
whos

%% Fields of |D|
% |D| is a |struct| array.  It contains the quadratic forms, and the upper/lower bounds constraints on
% their values.
fieldnames(D)

%% Quadratic forms
% The field |Qform| in |D| represents the quadratic form.  It is a
% symmetric matrix.  Not all rows/columns are given, since each quadratic
% form is only a function of 6-15 of the X's.  Hence the quadratic form
% matrices are sparse and low rank.  For example, the fourth quadratic form
% only depends on 5 of the X coordinates.   Hence the matrix is 6-by-6, and
% the |Pidx| field indicates which X's are involved.
D(4).Qform
D(4).Pidx

%% Create a random Xu, with values between -1 and 1, uniform
Xu = 2*(rand(102,1)-0.5);

%% Check feasability of Xu (randomly chosen in unit cube)
% Most likely, Xu is not feasible on all of the inequalities, but will be feasible
% on some.  Check the first 10.  The field |Interval| is a 1-by-2, which has
% the target value (as an interval) for the quadratic form.  If the value
% lies in the interval defined by [Interval(1) Interval(2)], then x is
% feasible for that given inequality.
for i=1:10
   OneAboveX = [1;Xu(D(i).Pidx)];  % [1;x], only using the relevant coordinates
   ValueAtX = OneAboveX'*D(i).Qform*OneAboveX;  % [1;x]'*M*[1;x]
   [D(i).Interval(1) ValueAtX D(i).Interval(2)]
end
%
% When I run this, it is usually only feasible on 2 or 3 of the first 10.

%% Check that |Xfeas| is indeed feasible
% First, the entries of |Xfeas| need to be within -1 and 1
all(Xfeas>=-1)
all(Xfeas<=1)
%%
% Next, verify that the quadratic forms are all feasible.
feas = false(length(D),1);
for i=1:length(D)
    vec = [1;Xfeas(D(i).Pidx)];
    valTmp = vec'*D(i).Qform*vec;
    feas(i) = D(i).Interval(1)<=valTmp  &  valTmp<=D(i).Interval(2);
end
all(feas)
%%
% Hence the feasible set is nonempty.  Can we get a "decent" convex
% outer-bound of the set?
