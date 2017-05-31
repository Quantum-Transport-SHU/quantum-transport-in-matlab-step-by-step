function [Q,fcnt] = my_quad(funfcn,a,b,tol,trace,varargin)
%QUAD   Numerically evaluate integral, adaptive Simpson quadrature.
%   Q = QUAD(FUN,A,B) tries to approximate the integral of scalar-valued
%   function FUN from A to B to within an error of 1.e-6 using recursive
%   adaptive Simpson quadrature. FUN is a function handle. The function
%   Y=FUN(X) should accept a vector argument X and return a vector result
%   Y, the integrand evaluated at each element of X.
%
%   Q = QUAD(FUN,A,B,TOL) uses an absolute error tolerance of TOL
%   instead of the default, which is 1.e-6.  Larger values of TOL
%   result in fewer function evaluations and faster computation,
%   but less accurate results.  The QUAD function in MATLAB 5.3 used
%   a less reliable algorithm and a default tolerance of 1.e-3.
%
%   Q = QUAD(FUN,A,B,TOL,TRACE) with non-zero TRACE shows the values
%   of [fcnt a b-a Q] during the recursion. Use [] as a placeholder to
%   obtain the default value of TOL.
%
%   [Q,FCNT] = QUAD(...) returns the number of function evaluations.
%
%   Use array operators .*, ./ and .^ in the definition of FUN
%   so that it can be evaluated with a vector argument.
%
%   Notes:
%   QUAD may be most efficient for low accuracies with nonsmooth
%   integrands.
%   QUADL may be more efficient than QUAD at higher accuracies
%   with smooth integrands.
%   QUADGK may be most efficient for oscillatory integrands and any smooth
%   integrand at high accuracies. It supports infinite intervals and can
%   handle moderate singularities at the endpoints. It also supports
%   contour integration along piecewise linear paths.
%   QUADV vectorizes QUAD for array-valued FUN.
%
%   Example:
%      Q = quad(@myfun,0,2);
%   where the file myfun.m defines the function:
%      %-------------------%
%      function y = myfun(x)
%      y = 1./(x.^3-2*x-5);
%      %-------------------%
%
%   or, use a parameter for the constant:
%      Q = quad(@(x)myfun2(x,5),0,2);
%   where the file myfun2.m defines the function:
%      %----------------------%
%      function y = myfun2(x,c)
%      y = 1./(x.^3-2*x-c);
%      %----------------------%
%
%   Class support for inputs A, B, and the output of FUN:
%      float: double, single
%
%   See also QUADV, QUADL, QUADGK, QUAD2D, DBLQUAD, TRIPLEQUAD, TRAPZ, FUNCTION_HANDLE.

%   Based on "adaptsim" by Walter Gander.
%   http://www.inf.ethz.ch/personal/gander
%
%   Reference:
%   [1] W. Gander and W. Gautschi, Adaptive Quadrature - Revisited,
%       BIT Vol. 40, No. 1, March 2000, pp. 84-101.
%
%   Copyright 1984-2009 The MathWorks, Inc.
%   $Revision: 5.26.4.12 $  $Date: 2010/08/23 23:09:11 $

f = fcnchk(funfcn);
if nargin < 4 || isempty(tol), tol = 1.e-6; end;
if nargin < 5 || isempty(trace), trace = 0; end;
if ~isscalar(a) || ~isscalar(b)
    error(message('MATLAB:quad:scalarLimits'));
end

% Initialize with three unequal subintervals.
h = 0.13579*(b-a);
x = [a a+h a+2*h (a+b)/2 b-2*h b-h b];
y = f(x, varargin{:});
fcnt = 7;
if size(y,1) ~= fcnt
    error(message('MATLAB:quad:funNotVectorized'));
end

% Fudge endpoints to avoid infinities.
if ~isfinite(max(max(abs(y(1,:,:)))))
    y(1,:,:) = f(a+eps(superiorfloat(a,b))*(b-a),varargin{:});
    fcnt = fcnt+1;
end
if ~isfinite(max(max(abs(y(7,:,:)))))
    y(7,:,:) = f(b-eps(superiorfloat(a,b))*(b-a),varargin{:});
    fcnt = fcnt+1;
end

% Call the recursive core integrator.
hmin = eps(b-a)/1024;
[Q(1,:,:),fcnt,warn(1)] = ...
    my_quadstep(f,x(1),x(3),y(1,:,:),y(2,:,:),y(3,:,:),tol,trace,fcnt,hmin,varargin{:});
[Q(2,:,:),fcnt,warn(2)] = ...
    my_quadstep(f,x(3),x(5),y(3,:,:),y(4,:,:),y(5,:,:),tol,trace,fcnt,hmin,varargin{:});
[Q(3,:,:),fcnt,warn(3)] = ...
    my_quadstep(f,x(5),x(7),y(5,:,:),y(6,:,:),y(7,:,:),tol,trace,fcnt,hmin,varargin{:});
Q = sum(Q);
warn = max(warn);

switch warn
    case 1
        warning(message('MATLAB:quad:MinStepSize'))
    case 2
        warning(message('MATLAB:quad:MaxFcnCount'))
    case 3
        warning(message('MATLAB:quad:ImproperFcnValue'))
    otherwise
        % No warning.
end

% ------------------------------------------------------------------------

function [Q,fcnt,warn] = my_quadstep(f,a,b,fa,fc,fb,tol,trace,fcnt,hmin,varargin)
%QUADSTEP  Recursive core routine for function QUAD.

maxfcnt = 10000;

% Evaluate integrand twice in interior of subinterval [a,b].
h = b - a;
c = (a + b)/2;
x = [(a + c)/2, (c + b)/2];
y = f(x, varargin{:});
fcnt = fcnt + 2;
fd = y(1,:,:);
fe = y(2,:,:);

% Three point Simpson's rule.
Q1 = (h/6)*(fa + 4*fc + fb);

% Five point double Simpson's rule.
Q2 = (h/12)*(fa + 4*fd + 2*fc + 4*fe + fb);

% One step of Romberg extrapolation.
Q = Q2 + (Q2 - Q1)/15;

if trace
    fprintf('%8.0f %16.10f %18.8e %16.10f\n',fcnt,a,h,max(max(max(Q))));
end

% Check termination criteria.
if ~isfinite(max(max(max(Q))))
    % Infinite or Not-a-Number function value encountered.
    warn = 3;
    return
end
if fcnt > maxfcnt
    % Maximum function count exceeded; singularity likely.
    warn = 2;
    return
end
if max(max(max(abs(Q2 - Q)))) <= tol
    % Accuracy over this subinterval is acceptable.
    warn = 0;
    return
end
if abs(h) < hmin || c == a || c == b
    % Minimum step size reached; singularity possible.
    warn = 1;
    return
end

% Subdivide into two subintervals.
[Qac,fcnt,warnac] = my_quadstep(f,a,c,fa,fd,fc,tol,trace,fcnt,hmin,varargin{:});
[Qcb,fcnt,warncb] = my_quadstep(f,c,b,fc,fe,fb,tol,trace,fcnt,hmin,varargin{:});
Q = Qac + Qcb;
warn = max(warnac,warncb);