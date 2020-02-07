%%% compare the one option small example using bigger or smaller penalty
SysDir = 'C:\Users\Yuting\OneDrive - UCL\Demand Response Business Model\Julia\PSPToyExample\OneOptionExampleADMMNewDatav2\';
DirSmallRHO = strcat(SysDir, 'SmallRHO');
DirBigRHO = strcat(SysDir, 'BigRHO');
Niterrations = 50;
NGenerators = 2; 

%% read the data of small rho
cd (DirSmallRHO)
file_production =  dir('production*.csv');
file_supply_R_X = dir('supplyR_X*.csv');
file_supply_R_Z = dir('supplyR_Z*.csv');
file_supply_R_U = dir('supplyR_U*.csv');
file_totalCost_X = dir('totalCost_X*.csv');
file_totalCost_Z = dir('totalCost_Z*.csv');
file_totalCost_U = dir('totalCost_U*.csv');

production_X_SmallRHO = csvread(file_production(1).name,1,0,[1,0,NGenerators,Niterrations-1]);
supply_R_X_SmallRHO = csvread(file_supply_R_X(1).name,1,0,[1,0,2,Niterrations-1]);
supply_R_Z_SmallRHO = csvread(file_supply_R_Z(1).name,1,0,[1,0,2,Niterrations-1]);
supply_R_U_SmallRHO = csvread(file_supply_R_U(1).name,1,0,[1,0,2,Niterrations-1]);
totalCost_X_SmallRHO = csvread(file_totalCost_X(1).name,1,0,[1,0,2,Niterrations-1]);
totalCost_Z_SmallRHO = csvread(file_totalCost_Z(1).name,1,0,[1,0,2,Niterrations-1]);
totalCost_U_SmallRHO = csvread(file_totalCost_U(1).name,1,0,[1,0,2,Niterrations-1]);

%% read the data of big rho
cd (DirBigRHO)
file_production =  dir('production*.csv');
file_supply_R_X = dir('supplyR_X*.csv');
file_supply_R_Z = dir('supplyR_Z*.csv');
file_supply_R_U = dir('supplyR_U*.csv');
file_totalCost_X = dir('totalCost_X*.csv');
file_totalCost_Z = dir('totalCost_Z*.csv');
file_totalCost_U = dir('totalCost_U*.csv');

production_X_BigRHO = csvread(file_production(1).name,1,0,[1,0,NGenerators,Niterrations-1]);
supply_R_X_BigRHO = csvread(file_supply_R_X(1).name,1,0,[1,0,2,Niterrations-1]);
supply_R_Z_BigRHO = csvread(file_supply_R_Z(1).name,1,0,[1,0,2,Niterrations-1]);
supply_R_U_BigRHO = csvread(file_supply_R_U(1).name,1,0,[1,0,2,Niterrations-1]);
totalCost_X_BigRHO = csvread(file_totalCost_X(1).name,1,0,[1,0,2,Niterrations-1]);
totalCost_Z_BigRHO = csvread(file_totalCost_Z(1).name,1,0,[1,0,2,Niterrations-1]);
totalCost_U_BigRHO = csvread(file_totalCost_U(1).name,1,0,[1,0,2,Niterrations-1]);


%%

% figure
% plot(supply_R_X_BigRHO)
% hold on;
% plot(supply_R_X_SmallRHO)
% title('Supply_R_X')
% 
% figure 
% plot(supply_R_Z_BigRHO)
% hold on;
% plot(supply_R_Z_SmallRHO)
% title('Supply_Z_X')



%% plot the quadratic function
MC1 = 1;
MC2 = 7;
omega = 0.8; % only available in the first period
VB = [6 10];
TotalCostScale = 5;
%omega*(MC1*d1 + MC2*d2) - omega*(d1+d2)*0.5*(VB(1)+VB(2))
% +0.5* RHO*(MC1*d1 + MC2*d2 - totalCost_Z_sol + totalCost_U_sol)^2
% + 0.5*RHO*(d1+d2 - supplyR_Z_sol + supplyR_U_sol)^2
% 0.8d1 + 5.6d2 - 6.4d1-6.4d2
% -5.6d1-0.8d2

iteration = 30;
d1L = 1.9;
d1U = 2;
d2L = 2.05;
d2U = 2.2;
res = 0.0001;

% d1L = 0;
% d1U = 2;
% d2L = 0;
% d2U = 5;
% res = 0.1;

%% big rho
RHO = 100;

% totalCost_Z_sol = totalCost_Z_SmallRHO(:,iteration)/TotalCostScale;
% totalCost_U_sol = totalCost_U_SmallRHO(:,iteration)/TotalCostScale;
% supplyR_Z_sol = supply_R_Z_SmallRHO(:,iteration);
% supplyR_U_sol = supply_R_U_SmallRHO(:,iteration);

totalCost_Z_sol = totalCost_Z_BigRHO(:,iteration)/TotalCostScale;
totalCost_U_sol = totalCost_U_BigRHO(:,iteration)/TotalCostScale;
supplyR_Z_sol = supply_R_Z_BigRHO(:,iteration);
supplyR_U_sol = supply_R_U_BigRHO(:,iteration);

fcn = @(d1,d2) -5.6*d1-0.8*d2 + ...
0.5*RHO*((d1+7*d2)/TotalCostScale - totalCost_Z_sol(1) + totalCost_U_sol(1) ).^2 + ...
0.5*RHO*(-totalCost_Z_sol(2) + totalCost_U_sol(2) ).^2 + ...
0.5*RHO*(d1+d2 - supplyR_Z_sol(1) + supplyR_U_sol(1) ).^2 + ...
0.5*RHO*(-supplyR_Z_sol(2) + supplyR_U_sol(2) ).^2 ;

[X,Y] = meshgrid(d1L:res:d1U, d2L:res:d2U);
figure(1)
meshc(X, Y, fcn(X,Y))
%grid on
%meshc(X, Y)
xlabel('p_{1,1}')
ylabel('p_{2,1}')
xlim([d1L d1U])
ylim([d2L d2U])

Z = fcn(X,Y);
min(min(Z));
[row, col] = find(Z==min(min(Z)));

%% small rho
RHO = 10;

totalCost_Z_sol = totalCost_Z_SmallRHO(:,iteration)/TotalCostScale;
totalCost_U_sol = totalCost_U_SmallRHO(:,iteration)/TotalCostScale;
supplyR_Z_sol = supply_R_Z_SmallRHO(:,iteration);
supplyR_U_sol = supply_R_U_SmallRHO(:,iteration);


fcn = @(d1,d2) -5.6*d1-0.8*d2 + ...
0.5*RHO*((d1+7*d2)/TotalCostScale - totalCost_Z_sol(1) + totalCost_U_sol(1) ).^2 + ...
0.5*RHO*(-totalCost_Z_sol(2) + totalCost_U_sol(2) ).^2 + ...
0.5*RHO*(d1+d2 - supplyR_Z_sol(1) + supplyR_U_sol(1) ).^2 + ...
0.5*RHO*(-supplyR_Z_sol(2) + supplyR_U_sol(2) ).^2 ;

[X,Y] = meshgrid(d1L:res:d1U, d2L:res:d2U);
figure(2)
meshc(X, Y, fcn(X,Y))
%grid on
%meshc(X, Y)
xlabel('p_{1,1}')
ylabel('p_{2,1}')
xlim([d1L d1U])
ylim([d2L d2U])


%% plot the coefficient of d1 and d2

















%% solution of NEXT iteration
d1 = 1.92616;
d2 = 2.13078; 

optimal = -5.6*d1-0.8*d2 +...
    0.5*RHO*((d1+7*d2)/TotalCostScale - totalCost_Z_sol(1) + totalCost_U_sol(1) )^2 + ...
0.5*RHO*(-totalCost_Z_sol(2) + totalCost_U_sol(2) )^2 + ...
0.5*RHO*(d1+d2 - supplyR_Z_sol(1) + supplyR_U_sol(1) )^2 + ...
0.5*RHO*(-supplyR_Z_sol(2) + supplyR_U_sol(2) )^2; 



