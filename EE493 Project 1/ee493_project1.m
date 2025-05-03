% topic:    EE493 Project 1 Script
% author:   Austin Driggs
% created:  2025-01-25


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD AND FORMAT I-V DATA

% data has to be in the working directory
load pnjunction_2025; 
V = vi(:, 1);
I = vi(:, 2);

% get positive I data set for log calculations
log_calculations = I > 0;
V_log = V(log_calculations);
I_log = I(log_calculations);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DETERMINE RANGE OF EXPONENTIAL OPERATION

% define exponential range based on trial and error!!!
exp_range = (V_log >= 0.4) & (V_log <= 0.8);

% refit
V_exp = V_log(exp_range);
I_exp = I_log(exp_range);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT RAW DATA

% plot 1: linear-linear plot
figure;
plot(V, I, 'o', 'MarkerSize', 5);
xlabel('Voltage (Volts)');
ylabel('Current (Amps)');
title('Plot 1: Linear-Linear Plot: Voltage vs Current');
grid on;
yl = ylim;
hold on;
plot([V_exp(1) V_exp(1)], yl, 'k--');
plot([V_exp(end) V_exp(end)], yl, 'k--');

% plot 2: log-linear plot
figure;
semilogy(V_log, I_log, 'o', 'MarkerSize', 5);
xlabel('Voltage (Volts)');
ylabel('Log(Current) (Amps)');
title('Plot 2: Log-Linear Plot: Voltage vs Log(Current)');
grid on;
yl = ylim;
hold on;
plot([V_exp(1) V_exp(1)], yl, 'k--');
plot([V_exp(end) V_exp(end)], yl, 'k--');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CURVE FITTING AND PARAMETER CALCULATIONS

% fit a straight line to log(I) vs V in exp range
A = polyfit(V_exp, log(I_exp), 1);
slope = A(1);
intercept = A(2);

% calculate parameters
I_0 = exp(intercept);
q = 1.602e-19;
k = 1.38e-23;
T = 300;
n = (q / (k * T)) / slope;
disp(['I_0 = ', num2str(I_0), ' Amps']);
disp(['n = ', num2str(n)]);

% voltage-to-current equation
I_fit = @(V) I_0 * exp((q * V) / (n * k * T)); % anon function
I_fitted = I_fit(V); % call for all V values

% diode equation
I_diode = @(V) I_0 * (exp(V / (25.86e-3 * n)) - 1); % for plot 5


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% CALCULATE PERCENT ERROR

% percent error within the defined exponential range
I_exp_fitted = I_fit(V_exp); 
percent_err = abs((I_exp - I_exp_fitted) ./ I_exp) * 100; 
mean_percent_err = mean(percent_err); 
disp(['percent error in exp range: ', num2str(mean_percent_err), ' %']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FITTED CURVE WITH DATA

% plot 3: data and fitted curve on a linear-linear plot with highlighted exponential range
figure;
plot(V, I, 'o', 'MarkerSize', 5); % original
hold on;
plot(V_exp, I_exp, 'ro', 'MarkerSize', 5); % highlited exponential
plot(V, I_fitted, 'r-', 'LineWidth', 2); % fitted curve
xlabel('Voltage (Volts)');
ylabel('Current (Amps)');
title('Plot 3: Linear-Linear Plot with Fitted Curve');
grid on;
legend(['I_0 = ', num2str(I_0, '%.2e'), ' A, n = ', num2str(n, '%.2f')], 'Location', 'Best');
yl = ylim;
plot([V_exp(1) V_exp(1)], yl, 'k--');
plot([V_exp(end) V_exp(end)], yl, 'k--');

% plot 4: highlight exponential range on the log-linear plot
figure;
semilogy(V_log, I_log, 'o', 'MarkerSize', 5); % original
hold on;
semilogy(V_exp, I_exp, 'ro', 'MarkerSize', 5); % highlited exponential
semilogy(V_log, I_fit(V_log), 'r-', 'LineWidth', 2);  % fitted curve
xlabel('Voltage (Volts)');
ylabel('Log(Current) (Amps)');
title('Plot 4: Log-Linear Plot with Highlighted Exponential Range');
grid on;
legend(['I_0 = ', num2str(I_0, '%.2e'), ' A, n = ', num2str(n, '%.2f')], 'Location', 'Best');
yl = ylim;
plot([V_exp(1) V_exp(1)], yl, 'k--');
plot([V_exp(end) V_exp(end)], yl, 'k--');

% plot 5: data and fitted curve on a linear-linear plot with diode equation
figure;
plot(V, I, 'o', 'MarkerSize', 5); % original
hold on;
plot(V_exp, I_exp, 'ro', 'MarkerSize', 5); % highlited exponential
plot(V, I_diode(V), 'g-', 'LineWidth', 2); % diode equation
xlabel('Voltage (Volts)');
ylabel('Current (Amps)');
title('Plot 5: Linear-Linear Plot with Diode Equation');
grid on;
legend(['I_0 = ', num2str(I_0, '%.2e'), ' A, n = ', num2str(n, '%.2f')], 'Location', 'Best');
yl = ylim;
plot([V_exp(1) V_exp(1)], yl, 'k--');
plot([V_exp(end) V_exp(end)], yl, 'k--');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE OUTPUTS

% print -dpng 'plot1_linear-linear.png';
% print -dpng 'plot2_log-linear.png';
% print -dpng 'plot3_exp-range.png';
% print -dpng 'plot4_fitted-curve.png';
