% topic:    EE493 Project 3 Script
% author:   Austin Driggs and Caleb Edwards
% created:  2025-01-29

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD AND FORMAT I-V DATA

% data has to be in the working directory
load mosfet_2025.mat;

% extract data into variables
nfet_gs = nfet_gatesweep; % [Vg, Id]
nfet_ds = nfet_drainsweeps; % [Vd, Id (for different Vg)]
pfet_gs = pfet_gatesweep;
pfet_ds = pfet_drainsweeps;
Vth = 0.6530; % added in after the calculation ran
%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 1: N-FET CHARACTERIZATION

%% Drain Sweep Plot

% Extract parameters through curve fitting for nFET
Ut = 25.86e-3;
idx = (nfet_gs(:,1) >= 0) & (nfet_gs(:,1) <= 0.8); % exponential range 0.1 to 0.8 V
p = polyfit(nfet_gs(idx,1), log(nfet_gs(idx,2)), 1);
kappa = p(1) * Ut;
I0 = exp(p(2));

% print extracted values
fprintf('\nnFET Parameters:\n');
fprintf('Kappa: %.4f\n', kappa);
fprintf('I0: %.4e A\n', I0);

% Drain sweep analysis
figure;
hold on;

% Define the saturation region where drain voltage is higher than threshold
saturation_region = nfet_ds(:,1) > Vth;  % Identify the saturation region

% Plot original drain sweeps
for i = 2:size(nfet_ds,2)
    plot(nfet_ds(:,1), nfet_ds(:,i));  % Original drain sweep curve
end
xline(0.65, 'k--', 'LineWidth', 1); % end of exponential region at 0.9V
annotation('textbox', [0.3, 0.6, 0.1, 0.1], 'String', 'Vth = 0.65V' , ...
           'EdgeColor', 'k', 'BackgroundColor', 'w', 'FontSize', 10);
xlabel('Drain Voltage (V)');
ylabel('Drain Current (A)');
title('nFET Drain Sweeps');
grid on;
hold off;


%% Log Gate Sweep
% gate sweep analysis - log
figure;
semilogy(nfet_gs(:,1), nfet_gs(:,2), 'o'); 
xlabel('Gate Voltage (V)');
ylabel('Drain Current (A)');
title('nFET Gate Sweep (Log Plot)');
grid on;

% Mark exponential region with vertical dotted lines (0.1V to 0.9V)
hold on;
xline(0.1, 'k--', 'LineWidth', 1); % start of exponential region at 0.1V
xline(0.8, 'k--', 'LineWidth', 1); % end of exponential region at 0.9V
hold off;

annotation('textbox', [0.4, 0.6, 0.2, 0.15], 'String', ...
    sprintf('Exponential Range = 0.1V - 0.8V\nI0 = %.2e\nKappa = %.4f', I0, kappa), ...
    'EdgeColor', 'k', 'BackgroundColor', 'w', 'FontSize', 10, 'FitBoxToText', 'on');

% Fit the exponential range (0.1V to 0.5V)
exp_range_idx = nfet_gs(:,1) >= 0.1 & nfet_gs(:,1) <= 0.8;
p_exp = polyfit(nfet_gs(exp_range_idx, 1), log(nfet_gs(exp_range_idx, 2)), 1);
I_exp_fit = exp(p_exp(2)) * exp(p_exp(1) * nfet_gs(exp_range_idx, 1));

% Plot the exponential fit
hold on;
plot(nfet_gs(exp_range_idx, 1), I_exp_fit, 'r-', 'LineWidth', 2, 'DisplayName', 'Exp Fit');
hold off;



%% Linear Gate Sweep
% gate sweep analysis - linear
figure;
plot(nfet_gs(:,1), nfet_gs(:,2), 'o'); 
xlabel('Gate Voltage (V)');
ylabel('Drain Current (A)');
title('nFET Gate Sweep (Linear Plot)');
grid on;

% Mark exponential region with vertical dotted lines (0.1V to 0.5V)
hold on;
xline(0.1, 'k--', 'LineWidth', 1); % start of exponential region at 0.1V
xline(0.8, 'k--', 'LineWidth', 1); % end of exponential region at 0.5V
hold off;

annotation('textbox', [0.4, 0.6, 0.2, 0.15], 'String', ...
    sprintf('Exponential Range = 0.1V - 0.8V\nI0 = %.2e\nKappa = %.4f', I0, kappa), ...
    'EdgeColor', 'k', 'BackgroundColor', 'w', 'FontSize', 10, 'FitBoxToText', 'on');

% Fit the exponential range (0.1V to 0.5V)
p_exp_lin = polyfit(nfet_gs(exp_range_idx, 1), nfet_gs(exp_range_idx, 2), 1);
I_exp_fit_lin = polyval(p_exp_lin, nfet_gs(exp_range_idx, 1));

% Plot the exponential fit
hold on;
plot(nfet_gs(exp_range_idx, 1), I_exp_fit_lin, 'r-', 'LineWidth', 2, 'DisplayName', 'Exp Fit');
hold off;




%% Calculations For Threshold Voltage

% Take the Linear Drain Current and Gate Voltage Graph 
% Take the square root of the drain current and extend that line to the
% x-intercept which should be our Vth


%Square Root of Drain Current vs Gate Voltage (Vgs > end of exponential range 0.8)

% Filter data for Vgs > 0.8V
valid_idx = nfet_gs(:,1) > 0.7; 
Vgs_filtered = nfet_gs(valid_idx, 1);
sqrt_Id_filtered = sqrt(nfet_gs(valid_idx, 2));

% Plot the square root of drain current vs gate voltage
figure;
plot(Vgs_filtered, sqrt_Id_filtered, 'o', 'DisplayName', 'Filtered Data');
xlabel('Gate Voltage (V)');
ylabel('Square Root of Drain Current (A^{1/2})');
title('nFET Gate Sweep (Filtered for V_{GS} > 0.1V)');
grid on;
hold on;

% Fit a line to the filtered data
p_sqrt_fit = polyfit(Vgs_filtered, sqrt_Id_filtered, 1);
sqrt_fit_line = polyval(p_sqrt_fit, Vgs_filtered);

% Plot the fitted line
plot(Vgs_filtered, sqrt_fit_line, 'r-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');

% Find the x-intercept (threshold voltage)
Vth = -p_sqrt_fit(2) / p_sqrt_fit(1); % x-intercept formula = about 0.65V

% Extend the fit line to x-axis
Vgs_extend = linspace(min(Vgs_filtered), Vth, 100);
sqrt_fit_extend = polyval(p_sqrt_fit, Vgs_extend);
plot(Vgs_extend, sqrt_fit_extend, 'b--', 'LineWidth', 1, 'DisplayName', 'Extended Fit');

% Mark the threshold voltage on the plot
xline(Vth, 'k--', 'LineWidth', 1, 'Label', sprintf('V_{th} = %.3f V', Vth));

legend;
hold off;

% Display threshold voltage
fprintf('\nThreshold Voltage (Vth): %.4f V\n', Vth);



%% Gate Voltages Computation for Each Drain Current (I think we should use the Voltage by reference instead below)

% Load drain sweep data
nfet_ds = nfet_drainsweeps; % Assuming nfet_drainsweeps contains [Vd, Id1, Id2, ..., Idn] data

% Extract drain currents for each gate voltage (assumed to be in columns 2 to 7)
ID_all = nfet_ds(:, 2:end); 

% Preallocate array to store calculated gate voltages
VG_all = zeros(1, size(ID_all, 2));

% Loop through each column of currents to calculate the gate voltage
for i = 1:size(ID_all, 2)
    % Get the average drain current for the ith column
    ID_avg = mean(ID_all(:, i)); % Alternatively, select a specific current if needed
    
    % Compute the gate voltage for this average current
    VG_all(i) = Vth + kappa * Ut * log(exp(sqrt(ID_avg / I0)) - 1);
end

% Display the calculated gate voltages
fprintf('\nCalculated Gate Voltages (V) for Each Column of Drain Currents:\n')
disp(VG_all)

%% Gate Voltages By Reference ( When Drain Voltage is 3.3V)
% 3.30000000000000	2.50648400000000e-09	4.76503700000000e-09	9.60619400000000e-09	0.000150169800000000	0.000303766800000000	0.000506068300000000
% Above is Column 331 of the drainsweep which is the drain currents when
% the drain voltage is 3.3V (Same in both drainsweep and gatesweep now)
% If we look at the currents when VD is 3.3V then we can reference them on
% the Gatesweep Datasheet to find the VG for each current value 
% Starting with column 2 the gate voltages are approximately:
% 0.428, 0.435, 0.476, 1.375, 1.65, 1.87
% Using this data we can see columns 2-4 are subthreshold and 5-6 are above
disp('Referenced Gate Voltages for Each Current')
disp('0.415, 0.448, 0.578, 1.4, 1.7, 2')





%% Subthreshold nFet Drainsweep Plot
figure;
hold on;

% Define the saturation region where drain voltage is higher than threshold
saturation_region = nfet_ds(:,1) > Vth;  % Identify the saturation region

% Initialize array to store Isat values
Isat_values = zeros(1, 3);
colors = lines(3); % MATLAB's "lines" colormap for distinct colors

% Plot subthreshold drain sweeps and find Isat at Vth
for i = 2:4  % Limit i to the range [2, 4]
    color_idx = i - 1;
    plot(nfet_ds(:,1), nfet_ds(:,i), 'Color', colors(color_idx, :), 'LineWidth', 1.5);  % Original drain sweep curve
    hold on; % Ensure multiple plots are retained
    
    % Find Isat as the current at Vth
    [~, idx] = min(abs(nfet_ds(:,1) - Vth)); % Find index closest to Vth
    Isat_values(i-1) = nfet_ds(idx, i);
end

xline(Vth, 'k--', 'LineWidth', 1, 'Label', sprintf('Vth = %.2f V', Vth));
xlabel('Drain Voltage (V)');
ylabel('Drain Current (A)');
title('Subthreshold nFET Drain Sweeps');
grid on;
hold off;

% Create a legend with Isat values
legend_entries = cell(1, length(Isat_values) + 1);
for i = 1:length(Isat_values)
    legend_entries{i} = sprintf('Case %d: I_{sat} = %.2e A', i, Isat_values(i));
end
legend_entries{end} = sprintf('Vth = %.2f V', Vth);
legend(legend_entries, 'Location', 'northoutside');

% Display Isat values
fprintf('\nSaturation Current (Isat) for each subthreshold case:\n');
for i = 1:3
    fprintf('Case %d: %.4e A\n', i, Isat_values(i));
end




%% Subthreshold Drainsweep Plot Traced To Find Early Voltage
figure;
hold on;

% Store VA values for each current
VA_values = [];

% Track min/max X for visibility adjustments
min_x = 0;
max_x = max(nfet_ds(:,1));

% Define a color map for consistency
colors = lines(3); % MATLAB's "lines" colormap generates distinguishable colors

% Plot subthreshold drain sweeps
fprintf('\nNegative Early Voltage Values From Tracing\n');
for i = 2:4  % Limit i to the range [2, 4]
    color_idx = i - 1; % Index to match colors
    
    % Plot the original drain sweep curve
    plot(nfet_ds(:,1), nfet_ds(:,i), 'Color', colors(color_idx, :), 'LineWidth', 1.5);  
    hold on; 
    
    % Extract data points after Vth (0.65V) for curve tracing
    idx = nfet_ds(:,1) > Vth;  % Identify the region where Vd > Vth
    Vd_sub = nfet_ds(idx,1);   % Extract corresponding Vd values
    Id_sub = nfet_ds(idx,i);   % Extract corresponding Id values
    
    % Perform linear regression (fit a line) on the selected region
    coeffs = polyfit(Vd_sub, Id_sub, 1); % Fit a first-degree polynomial (y = mx + b)
    slope = coeffs(1);
    intercept = coeffs(2);
    
    % Extend the fitted line backward to find x-intercept (-VA)
    VA = -intercept / slope; % Solve for x when y = 0
    VA_values = [VA_values; VA];  % Store VA for printing later
    min_x = min(min_x, VA); % Track minimum X value for full graph visibility
    
    % Generate points for the extended trend line
    Vd_ext = linspace(VA, max(Vd_sub), 100); % Extend from VA to max Vd
    Id_ext = polyval(coeffs, Vd_ext); % Compute corresponding Id values
    
    % Plot the extended fitted line in matching color
    plot(Vd_ext, Id_ext, '--', 'Color', colors(color_idx, :), 'LineWidth', 1.5); 
    
    % Print VA value to console
    fprintf('-V_A for I_d%d: %.2f V\n', i, VA);
end

% Mark the transition voltage with a vertical dashed line
xline(Vth, 'k--', 'LineWidth', 1);

% Annotation for threshold voltage
annotation('textbox', [0.3, 0.6, 0.1, 0.1], 'String', 'Vth = 0.65V' , ...
           'EdgeColor', 'k', 'BackgroundColor', 'w', 'FontSize', 10);

% Adjust axes for full visibility
xlim([min_x, max_x]);  % Ensure all curves and VA lines are visible

% Make x-axis bold
ax = gca;
ax.XAxis.LineWidth = 2; % Thicker x-axis

% Axis labels and title
xlabel('Drain Voltage (V)');
ylabel('Drain Current (A)');
title('Subthreshold nFET Drain Sweeps with Curve Tracing');
grid on;
hold off;

% Find the average of Early voltages
Vearly = mean(VA_values);
fprintf('\nEarly Voltage Average: %.4f V\n', Vearly);

% Create a legend with Early Voltage values
legend_entries = cell(1, length(VA_values) + 1);
for i = 1:length(VA_values)
    legend_entries{i} = sprintf('Case %d: V_A = %.2f V', i, VA_values(i));
end
legend_entries{end} = sprintf('Avg V_A: %.2f V', Vearly);
legend(legend_entries, 'Location', 'northeastoutside');





%% Calculate Isat for each subthreshold case
% Isat_values = zeros(1, 3);
% for i = 2:4  % Columns corresponding to subthreshold cases
%     saturation_region = nfet_ds(:,1) > Vth; % Identify saturation region
%     Isat_values(i-1) = mean(nfet_ds(saturation_region, i));
% end
% 
% fprintf('\nSaturation Current (Isat) for each subthreshold case:\n');
% for i = 1:3
%     fprintf('Case %d: %.4e A\n', i, Isat_values(i));
% end


%% Drain Sweeps – Transition from Ohmic to Saturation
% Provide a plot that indicates the transition from ohmic to saturation for each subthreshold case
figure;
hold on;
colors = lines(3); % Use distinct colors for each case
transition_points = zeros(1, 3);

for i = 2:4  % Columns corresponding to subthreshold cases
    dId_dVd = diff(nfet_ds(:,i)) ./ diff(nfet_ds(:,1)); % Compute derivative
    transition_idx = find(dId_dVd < 0.1 & nfet_ds(1:end-1,1) > 0.1, 1); % Identify transition point
    
    if ~isempty(transition_idx)
        transition_voltage = nfet_ds(transition_idx,1);
        transition_points(i-1) = transition_voltage;
        plot(nfet_ds(:,1), nfet_ds(:,i), 'Color', colors(i-1, :), 'DisplayName', sprintf('Case %d', i-1));
        xline(transition_voltage, '--', 'Color', colors(i-1, :), 'LineWidth', 1.5, 'Label', sprintf('%.2f V', transition_voltage));
    end
end

xlabel('Drain Voltage (V)');
ylabel('Drain Current (A)');
title('Ohmic to Saturation Transition for Subthreshold Cases');
legend;
grid on;
hold off;

% Discussion Output
fprintf('\nOhmic-to-Saturation Transition Voltages:\n');
for i = 1:3
    if transition_points(i) > 0
        fprintf('Case %d: %.4f V\n', i, transition_points(i));
    else
        fprintf('Case %d: No clear transition detected\n', i);
    end
end

fprintf('\nDiscussion:\n');
fprintf('The transition locations were determined by finding the first point where the derivative dId/dVd \n');
fprintf('drops below 0.1, indicating the shift from the ohmic region (linear increase) to saturation. \n');
fprintf('These transition voltages should theoretically be close to Vd = Vg - Vth, and we will compare \n');
fprintf('them to this expected value in further analysis.');





%% Make a table summarizing all drain-sweep characteristics
% this can probably just be done in the slideshow



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 2: P-FET CHARACTERIZATION

% Calculate V_GS (Gate-Source Voltage)
Vgs = pfet_gs(:,1) - 3.3;  % Assuming source voltage is 3.3V

%% Plot V_GS vs Drain Current
figure;
semilogy(Vgs, abs(pfet_gs(:,2)), 'o');
xlabel('Gate-Source Voltage (V)');
ylabel('|Drain Current| (A)');
title('pFET Gate-Source Voltage vs Drain Current (Log Plot)');
grid on;
% Set x-axis limits from -3.3V to 0V
xlim([-3.3 0]);
% Mark exponential region with vertical dotted lines (-0.8V to 0V)
hold on;
xline(-0.8, 'k--', 'LineWidth', 1); % start of exponential region at -0.8V
xline(0, 'k--', 'LineWidth', 1); % end of exponential region at 0V
hold off;

% Add labels to the exponential region box
annotation('textbox', [0.3, 0.6, 0.1, 0.1], 'String', 'Exponential Range = -0.8V - 0V', ...
           'EdgeColor', 'k', 'BackgroundColor', 'w', 'FontSize', 10);

%% Extract Parameters through Curve Fitting (Exponential Range -0.8V to 0V)
exp_range_idx = Vgs >= -0.8 & Vgs <= 0;  % Exponential range index for Vgs
p = polyfit(Vgs(exp_range_idx), log(abs(pfet_gs(exp_range_idx,2))), 1);

% Extract Kappa and I0
kappa_p = abs(p(1) * Ut);  % Take the absolute value of kappa
I0_p = exp(p(2));

% Print extracted values with absolute kappa
fprintf('\npFET Parameters:\n');
fprintf('Kappa: %.4f\n', kappa_p);
fprintf('I0: %.4e A\n', I0_p);

%% Gate Sweep Analysis with Log-Linear Plot and Curve Fit
figure;
semilogy(Vgs, abs(pfet_gs(:,2)), 'o');
xlabel('Gate-Source Voltage (V)');
ylabel('|Drain Current| (A)');
title('pFET Gate-Source Voltage vs Drain Current (Log Plot)');
grid on;
% Mark exponential region with vertical dotted lines (-0.8V to 0V)
hold on;
xline(-0.8, 'k--', 'LineWidth', 1); % start of exponential region at -0.8V
xline(0, 'k--', 'LineWidth', 1); % end of exponential region at 0V
hold off;

% Add labels to the exponential region box with I0 and Kappa
annotation('textbox', [0.3, 0.5, 0.25, 0.15], 'String', sprintf('Exponential Range = -0.8V - 0V\nI0 = %.2e\nKappa = %.4f', I0_p, kappa_p), ...
           'EdgeColor', 'k', 'BackgroundColor', 'w', 'FontSize', 10, 'FitBoxToText', 'on');

% Fit the exponential range (-0.8V to 0V)
I_exp_fit = exp(p(2)) * exp(p(1) * Vgs(exp_range_idx));

% Plot the exponential fit
hold on;
plot(Vgs(exp_range_idx), I_exp_fit, 'r-', 'LineWidth', 2, 'DisplayName', 'Exp Fit');
hold off;

%% Gate Sweep Analysis - Linear Plot
figure;
plot(Vgs, pfet_gs(:,2), 'o');
xlabel('Gate-Source Voltage (V)');
ylabel('Drain Current (A)');
title('pFET Gate-Source Voltage vs Drain Current (Linear Plot)');
grid on;
% Mark exponential region with vertical dotted lines (-0.8V to 0V)
hold on;
xline(-0.8, 'k--', 'LineWidth', 1); % start of exponential region at -0.8V
xline(0, 'k--', 'LineWidth', 1); % end of exponential region at 0V
hold off;

% Add labels to the exponential region box with I0 and Kappa
annotation('textbox', [0.3, 0.5, 0.25, 0.15], 'String', sprintf('Exponential Range = -0.8V - 0V\nI0 = %.2e\nKappa = %.4f', I0_p, kappa_p), ...
           'EdgeColor', 'k', 'BackgroundColor', 'w', 'FontSize', 10, 'FitBoxToText', 'on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
