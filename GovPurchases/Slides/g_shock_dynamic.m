function [residual, g1, g2, g3] = g_shock_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [nperiods by M_.exo_nbr] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   steady_state  [M_.endo_nbr by 1] double       vector of steady state values
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations.
%                                          Dynare may prepend auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(6, 1);
T11 = y(3)^(-params(3));
T15 = params(2)*y(9)^(-params(3));
T54 = y(6)^(1-params(1));
T60 = 1-params(6)+params(1)*(y(4)/y(10))^(params(1)-1);
T65 = (y(1)/y(6))^params(1);
T71 = y(1)^params(1);
lhs =T11;
rhs =T15*T60;
residual(1)= lhs-rhs;
lhs =T11*(1-params(1))*T65;
rhs =y(6)^params(4);
residual(2)= lhs-rhs;
lhs =y(3)+y(5)+y(8);
rhs =y(7);
residual(3)= lhs-rhs;
lhs =y(5);
rhs =y(4)-(1-params(6))*y(1);
residual(4)= lhs-rhs;
lhs =y(8);
rhs =params(5)*y(2)+x(it_, 1);
residual(5)= lhs-rhs;
lhs =y(7);
rhs =T54*T71;
residual(6)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(6, 11);

  %
  % Jacobian matrix
  %

T74 = getPowerDeriv(y(3),(-params(3)),1);
T82 = getPowerDeriv(y(1)/y(6),params(1),1);
T90 = getPowerDeriv(y(4)/y(10),params(1)-1,1);
  g1(1,3)=T74;
  g1(1,9)=(-(T60*params(2)*getPowerDeriv(y(9),(-params(3)),1)));
  g1(1,4)=(-(T15*params(1)*1/y(10)*T90));
  g1(1,10)=(-(T15*params(1)*T90*(-y(4))/(y(10)*y(10))));
  g1(2,3)=T65*(1-params(1))*T74;
  g1(2,1)=T11*(1-params(1))*1/y(6)*T82;
  g1(2,6)=T11*(1-params(1))*T82*(-y(1))/(y(6)*y(6))-getPowerDeriv(y(6),params(4),1);
  g1(3,3)=1;
  g1(3,5)=1;
  g1(3,7)=(-1);
  g1(3,8)=1;
  g1(4,1)=1-params(6);
  g1(4,4)=(-1);
  g1(4,5)=1;
  g1(5,2)=(-params(5));
  g1(5,8)=1;
  g1(5,11)=(-1);
  g1(6,1)=(-(T54*getPowerDeriv(y(1),params(1),1)));
  g1(6,6)=(-(T71*getPowerDeriv(y(6),1-params(1),1)));
  g1(6,7)=1;

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],6,121);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],6,1331);
end
end
end
end
