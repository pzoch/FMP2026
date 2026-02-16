function [residual, g1, g2, g3] = g_shock_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations.
%                                          Dynare may prepend or append auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g3        [M_.endo_nbr by (M_.endo_nbr)^3] double   Third derivatives matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 6, 1);

%
% Model equations
%

T11 = y(1)^(-params(3));
T19 = y(2)/y(4);
T23 = 1-params(6)+params(1)*T19^(params(1)-1);
T28 = T19^params(1);
T47 = y(4)^(1-params(1));
T48 = y(2)^params(1);
lhs =T11;
rhs =T11*params(2)*T23;
residual(1)= lhs-rhs;
lhs =T11*(1-params(1))*T28;
rhs =y(4)^params(4);
residual(2)= lhs-rhs;
lhs =y(1)+y(3)+y(6);
rhs =y(5);
residual(3)= lhs-rhs;
lhs =y(3);
rhs =y(2)-(1-params(6))*y(2);
residual(4)= lhs-rhs;
lhs =y(6);
rhs =y(6)*params(5)+x(1);
residual(5)= lhs-rhs;
lhs =y(5);
rhs =T47*T48;
residual(6)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(6, 6);

  %
  % Jacobian matrix
  %

T51 = getPowerDeriv(y(1),(-params(3)),1);
T58 = getPowerDeriv(T19,params(1)-1,1);
T63 = getPowerDeriv(T19,params(1),1);
  g1(1,1)=T51-T23*params(2)*T51;
  g1(1,2)=(-(T11*params(2)*params(1)*1/y(4)*T58));
  g1(1,4)=(-(T11*params(2)*params(1)*T58*(-y(2))/(y(4)*y(4))));
  g1(2,1)=T28*(1-params(1))*T51;
  g1(2,2)=T11*(1-params(1))*1/y(4)*T63;
  g1(2,4)=T11*(1-params(1))*T63*(-y(2))/(y(4)*y(4))-getPowerDeriv(y(4),params(4),1);
  g1(3,1)=1;
  g1(3,3)=1;
  g1(3,5)=(-1);
  g1(3,6)=1;
  g1(4,2)=(-(1-(1-params(6))));
  g1(4,3)=1;
  g1(5,6)=1-params(5);
  g1(6,2)=(-(T47*getPowerDeriv(y(2),params(1),1)));
  g1(6,4)=(-(T48*getPowerDeriv(y(4),1-params(1),1)));
  g1(6,5)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],6,36);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],6,216);
end
end
end
end
