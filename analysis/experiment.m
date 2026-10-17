clc, clearvars, close all

net = [2654; 2438; 2024; 1831; 1690];
wat_cm = [1; 2; 4; 5; 6];

figure(1);
hold on;
scatter(wat_cm, net, 'filled');
title("Original data");
f = fit(wat_cm, net, 'exp1');
plot(f, wat_cm, net);
disp(f);
xlabel("water thickness [cm]");
ylabel("NET measurements");
wat_coeff = 0.09158;
 
%alternate way
% net = I0*exp(-c*wat_cm) => c = -log(net/I0)/wat_cm
c12 = -log(net(2)/net(1));
c13 = -log(net(3)/net(1))/3;
c14 = -log(net(4)/net(1))/4;
c15 = -log(net(5)/net(1))/5;
c23 = -log(net(3)/net(2))/2;
c24 = -log(net(4)/net(2))/3;
c25 = -log(net(5)/net(2))/4;
c34 = -log(net(4)/net(3));
c35 = -log(net(5)/net(3))/2;
c45 = -log(net(5)/net(4));

alt_wat_coeff = (c12+c13+c14+c15+c23+c24+c25+c34+c35+c45) / 10;

%computation of the actual water coefficient from book
%supposing there is a linear regression model of 1st degree that fits the
%coeff - gamma energy data
a = 600;
a_coeff = 0.0896;
b = 800;
b_coeff = 0.0786;
linear_coeff = (a_coeff - b_coeff) / (a - b);
act_wat_coeff = a_coeff + linear_coeff*(662-a); % first-draft interpolation, superseded below
act_wat_coeff = 0.0896; % reference value at 662 keV, as derived in the report



fprintf("Deviation from actual value: %d\n", wat_coeff-act_wat_coeff);


