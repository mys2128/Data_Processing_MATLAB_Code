%BMEN E3810: Biomedical Engineering Laboratory I (INSTRON Data Analysis)
%Written by: Mia Raneri, Sofia Garcia, Andre Perez, & Miranda Salazar
%Fall 2024
%Initialize the Workspace:

clear; clc; close all;

%% ORGANIZE DATA
%Define Data Parameters:
materials = {'Silicone_Thin_Instron', 'Silicone_Thick_Instron', ...
    'Silicone_Thin_Minstron', 'Silicone_Thick_Minstron', ...
    'Skin_Raw', 'Skin_Treated'};
teamnames = {'Wed01',  'Wed02',  'Wed03',  'Wed04', ...
    'Thurs01','Thurs02','Thurs03', 'Thurs04', ...
    'Thurs05','Thurs06','Thurs07'};

filename = 'I.2 - Biomechanical Testing - Class Data.xlsx';

%Import Raw Data:
for k = 1:length(materials)
    rawdata.(materials{k}) = readcell(filename, 'Sheet', materials{k});
    rawdata.(materials{k})(cellfun(@(x) any(ismissing(x)), rawdata.(materials{k}))) = {NaN};
end

%Organize Data:
for k = 1:length(materials)
    for n = 1:length(teamnames)
        data.(materials{k}).L0(1,n)   = (rawdata.(materials{k})(3,4*n-1));
        data.(materials{k}).Width(1,n) = (rawdata.(materials{k})(4,4*n-1));
        data.(materials{k}).Thickness(1,n)  = (rawdata.(materials{k})(5,4*n-1));
        data.(materials{k}).Area(1,n)   = (rawdata.(materials{k})(6,4*n-1));
        
        data.(materials{k}).Pos{1,n}    = (rawdata.(materials{k})(9:end,4*n-2));
        data.(materials{k}).Force{1,n}  = (rawdata.(materials{k})(9:end,4*n-1));
    end
end

%some data points were weird
data.(materials{3}).Pos{1,2}{1} = 0.4;
data.(materials{3}).Pos{1,10}{27} = 5.25;
data.(materials{4}).Pos{1,2}{39} = 7.75;

%create arrays of all teams forces and positions:
AllForces = {};
AllPos = {};
for k = 1:length(materials)
    for n = 1:length(teamnames)
        Force = cell2mat(data.(materials{k}).Force{1,n});
        Pos  = cell2mat(data.(materials{k}).Pos{1,n});
       
        % Zero force/position data
        Force = Force - Force(1);
        Pos = Pos - Pos(1);

        %some forces were negative
        neg_thin_minstron = [2,4,5,8,9, 11];
        neg_thick_minstron = [2,4,5,7,8,9,11];
        if k == 3
            if ismember(n, neg_thin_minstron)
                Force = -1 * Force;
            end
        end

        if k == 4
            if ismember(n,neg_thick_minstron)
                Force = -1 * Force;
            end
        end

        AllForces{k,n} = Force;
        AllPos{k,n} = Pos;
    end
    
end


%Plots for Load vs. Displacement
silicone = {"Thin Instron","Thick Instron","Thin Minstron","Thick Minstron"};
chicken = {"Chicken Raw","Chicken Treated"};
colors4= {"r", "b", "g", "m"};
colors2 = {"m" , "b"};

for k = 1:4
    for n = 1:11
        figure(1);
        h(k) = plot(AllPos{k,n}, AllForces{k,n}, 'Color', colors4{k});
        hold on
        
        title(strcat("Load_vs_Displacement_(Silicone)"), 'Interpreter', 'None')
        legend(h, silicone, 'location','best')
        xlabel('Position Displacement (mm)');
        ylabel('Force (N)')
    end
end
hold off

for k = 5:6
    for n = 1:11
        figure(2);
        j(k-4) = plot(AllPos{k,n}, AllForces{k,n}, 'Color', colors2{k-4});
        hold on

        title(strcat("Load_vs_Displacement_(Chicken)"), 'Interpreter', 'None')
        legend(j, chicken, 'location','best')
        xlabel('Position Displacement (mm)');
        ylabel('Force (N)')
    end
end
hold off

%% STIFFNESS ANALYSIS

Stiffness_allmaterials = {};
for k = 1:6  
    for n = 1:length(teamnames)
        Force = AllForces{k,n};
        Pos = AllPos{k,n};
        
        % Determine if and where NaN entries exist in force/position data
        idx = find(isnan(Force), 1, 'first');
        if isempty(idx) == 0
            Force = Force(1:idx-1);
            Pos = Pos(1:idx-1);
        end

        %focus only on first half of data for better results
        midpoint = round(length(Force)/2);
        Force_subset = Force(1:midpoint);
        Pos_subset = Pos(1:midpoint);

        % Crop off all values outside of linear range
    
        % Take first and second derivatives of force with respect to position
        first_deriv = gradient(Force,Pos);      % First derivative
        second_deriv = diff(first_deriv);     % Second derivative
        
        std_dev_dd = std(second_deriv);
        
        % True/False evaluation: 1 if the absolute value of the second derivative 
        % is less than one standard deviation, 0 otherwise
        linear_region = abs(second_deriv) < 2*std_dev_dd;
        linear_region = [false; linear_region; false]; %ensures it starts and ends with 0
        % Find the longest consecutive region of ones (linear region)
        difference = diff(linear_region); %helps locate areas where linear regions are starting/ending
        start_region = find(difference == 1); %index of locations where regions start
        end_region = find(difference == -1); %index of locations where regions end
        end_region = end_region - 1 ; %want the index of the last 1, not the first 0
    
        %now we want to find the longest interval of consecutive 1s:
        region_length = end_region - start_region; %length 
        [~, longest_region_idx] = max(region_length);
        start_idx = start_region(longest_region_idx);
        end_idx = end_region(longest_region_idx);
       
        Force_linear = Force(start_idx:end_idx);
        Pos_linear = Pos(start_idx:end_idx);
    
        % Calculate average slope of force/position data to determine stiffness
        slope_linear_region = first_deriv(start_idx:end_idx);
        stiffness = mean(slope_linear_region)/1000; %divide by 1000 to convert mm to m
        Stiffness_allmaterials{k,n} = stiffness;

    end
   
end

thin_stiffness_instron = [Stiffness_allmaterials{1,:}];
thick_stiffness_instron = [Stiffness_allmaterials{2,:}];
thin_stiffness_minstron = [Stiffness_allmaterials{3,:}];
thick_stiffness_minstron = [Stiffness_allmaterials{4,:}];
rawchicken_k = [Stiffness_allmaterials{5,:}];
treatedchicken_k = [Stiffness_allmaterials{6,:}];

mean_k_thin_instron = mean(thin_stiffness_instron, 'omitnan');
std_k_thin_instron = std(thin_stiffness_instron, 'omitnan');

mean_k_thin_minstron = mean(thin_stiffness_minstron, 'omitnan');
std_k_thin_minstron = std(thin_stiffness_minstron, 'omitnan');

mean_k_thick_instron = mean(thick_stiffness_instron,'omitnan');
std_k_thick_instron = std(thick_stiffness_instron,'omitnan');

mean_k_thick_minstron = mean(thick_stiffness_minstron,'omitnan');
std_k_thick_minstron = std(thick_stiffness_minstron,'omitnan');

mean_k_rawchicken = mean(rawchicken_k, 'omitnan');
std_k_rawchicken = std(rawchicken_k, 'omitnan');

mean_k_treatedchicken = mean(treatedchicken_k, 'omitnan');
std_k_treatedchicken = std(treatedchicken_k, 'omitnan');

%% STRESS & STRAIN ANALYSIS

%again, master array of all stress and strains (like force & position)
AllStresses = {};
AllStrains = {};
E_allmaterials = {};
ultimate_tensile_stress = {};
extensibility = {};

for k = 1:6
    for n = 1:11
        Force = AllForces{k,n};
        Pos = AllPos{k,n};
      
        Pos = Pos / 1000; %convert to m
        A = cell2mat(data.(materials{k}).Area(1,n));
        area = A * 10^-6; %convert mm^2 to m^2
        L0 = cell2mat(data.(materials{k}).L0(1,n));
        length = L0/1000; % mm to m
        
        Strain = Pos / length;  % Compute Strain (strain = displacement / original length)
        Stress = Force / area; % Compute Stress (stress = force / cross-sectional area)
        Stress = Stress * 10^-6; %convert from Pa to MPa

        AllStresses{k,n} = Stress;
        AllStrains{k,n} = Strain;

        % Take first and second derivatives of stress with strain to position
        first_deriv = gradient(Stress,Strain);      % First derivative
        second_deriv = diff(first_deriv);     % Second derivative
        
        std_dev_dd = std(second_deriv);
        
        % True/False evaluation: 1 if the absolute value of the second derivative 
        % is less than one standard deviation, 0 otherwise
        linear_region = abs(second_deriv) < 2*std_dev_dd;
        linear_region = [false; linear_region; false]; %ensures it starts and ends with 0
        % Find the longest consecutive region of ones (linear region)
        difference = diff(linear_region); %helps locate areas where linear regions are starting/ending
        start_region = find(difference == 1); %index of locations where regions start
        end_region = find(difference == -1); %index of locations where regions end
        end_region = end_region - 1 ; %want the index of the last 1, not the first 0
    
        %now we want to find the longest interval of consecutive 1s:
        region_length = end_region - start_region; %length 
        [~, longest_region_idx] = max(region_length);
        start_idx = start_region(longest_region_idx);
        end_idx = end_region(longest_region_idx);
    
        idx2 = find(Stress == max(Stress), 1, 'first'); %this is max Stress

        ultimate_tensile_stress{k,n} = Stress(idx2);
        extensibility{k,n} = Strain(idx2) * 100;

        Stress_linear = Stress(start_idx:end_idx);
        Strain_linear = Strain(start_idx:end_idx);
    
        % Calculate average slope of stress/strain data to determine E
        slope_linear_region = first_deriv(start_idx:end_idx);     
        E = mean(slope_linear_region);
        E_allmaterials{k,n} = E;
    end
   
end
  
thin_E_instron = [E_allmaterials{1,:}];
thin_E_minstron = [E_allmaterials{3,:}];
thick_E_instron = [E_allmaterials{2,:}];
thick_E_minstron = [E_allmaterials{4,:}];
rawchicken_E = [E_allmaterials{5,:}];
treatedchicken_E = [E_allmaterials{6,:}];

mean_E_thin_instron = mean(thin_E_instron, 'omitnan');
std_E_thin_instron = std(thin_E_instron, 'omitnan');

mean_E_thin_minstron = mean(thin_E_minstron, 'omitnan');
std_E_thin_minstron = std(thin_E_minstron, 'omitnan');

mean_E_thick_instron = mean(thick_E_instron, 'omitnan');
std_E_thick_instron = std(thick_E_instron, 'omitnan');

mean_E_thick_minstron = mean(thick_E_minstron, 'omitnan');
std_E_thick_minstron = std(thick_E_minstron, 'omitnan');

mean_E_rawchicken = mean(rawchicken_E, 'omitnan');
std_E_rawchicken = std(rawchicken_E,'omitnan');

mean_E_treatedchicken = mean(treatedchicken_E, 'omitnan');
std_E_treatedchicken = std(treatedchicken_E, 'omitnan');

%plots for Stress vs. Strain:
for k = 1:4
    for n = 1:11
        figure(3);
        h(k) = plot(AllStrains{k,n}, AllStresses{k,n}, 'Color', colors4{k});
        hold on
        
        title(strcat("Stress_vs_Strain_(Silicone)"), 'Interpreter', 'None')
        legend(h, silicone, 'location','best')
        xlabel('Strain');
        ylabel('Stress (Pa)');
    end
end
hold off

for k = 5:6
    for n = 1:11
        figure(4);
        j(k-4) = plot(AllStrains{k,n}, AllStresses{k,n} * 10^-6, 'Color', colors2{k-4});
        hold on

        title(strcat("Stress_vs_Strain_(Chicken)"), 'Interpreter', 'None')
        legend(j, chicken, 'location','best')
        xlabel('Strain');
        ylabel('Stress (MPa)');
   end
end

%% Bar Plots for Material Properties
groupnames = {'Silicone Thin Instron', 'Silicone Thick Instron', 'Silicone Thin Minstron', 'Silicone Thick Minstron','Raw Chicken', 'Treated Chicken'};

%Young's Modulus (E):
Emeans = [mean_E_thin_instron; mean_E_thick_instron; mean_E_thin_minstron; mean_E_thick_minstron; mean_E_rawchicken; mean_E_treatedchicken];
Estds = [std_E_thin_instron; std_E_thick_instron; std_E_thin_minstron; std_E_thick_minstron; std_E_rawchicken; std_E_treatedchicken];

figure(5); 
bar(Emeans);
hold on;

x = 1:6;
errorbar(x, Emeans, Estds, 'k', 'linestyle', 'none'); 

set(gca, 'XTickLabel', groupnames); 
title('Youngs Modulus (E) by Material');
xlabel('Material');
ylabel('Youngs Modulus, E (MPa)');

for i = 1:6
    meanLabel = sprintf('Mean: %.2f', Emeans(i));
    stdLabel = sprintf('SD: %.2f', Estds(i));
    
    text(i, Emeans(i) + Estds(i) + 1.0, meanLabel, 'HorizontalAlignment', 'center');
    text(i, Emeans(i) + Estds(i) + 0.5, stdLabel, 'HorizontalAlignment', 'center');
end
hold off

% Stiffness:
kmeans = [mean_k_thin_instron; mean_k_thick_instron; mean_k_thin_minstron; mean_k_thick_minstron; mean_k_rawchicken; mean_k_treatedchicken];
kstds = [std_k_thin_instron; std_k_thick_instron; std_k_thin_minstron; std_k_thick_minstron; std_k_rawchicken; std_k_treatedchicken];

figure(6);
bar(kmeans);
hold on

x = 1:6;
errorbar(x, kmeans, kstds, 'k', 'linestyle', 'none'); 

set(gca, 'XTickLabel', groupnames); 
title('Stiffness (k) by Material');
xlabel('Material');
ylabel('Stiffness, k (N/m)');

for i = 1:6
    meanLabel = sprintf('Mean: %.5f', kmeans(i));
    stdLabel = sprintf('SD: %.5f', kstds(i));
    
    text(i, kmeans(i) + kstds(i) - 0.0002, meanLabel, 'HorizontalAlignment', 'center');
    text(i, kmeans(i) + kstds(i) - 0.0005, stdLabel, 'HorizontalAlignment', 'center');
end
hold off

%Extensibility:
thin_ex_instron = [extensibility{1,:}];
thick_ex_instron = [extensibility{2,:}];
thin_ex_minstron = [extensibility{3,:}];
thick_ex_minstron = [extensibility{4,:}];
raw_ex = [extensibility{5,:}];
treated_ex = [extensibility{6,:}];

exmeans = [mean(thin_ex_instron); mean(thick_ex_instron); mean(thin_ex_minstron); mean(thick_ex_minstron); mean(raw_ex); mean(treated_ex)];
exstds = [std(thin_ex_instron); std(thick_ex_instron); std(thin_ex_minstron); std(thick_ex_minstron); std(raw_ex); std(treated_ex)];

figure(7); 
bar(exmeans);
hold on;

x = 1:6;
errorbar(x, exmeans, exstds, 'k', 'linestyle', 'none'); 

set(gca, 'XTickLabel', groupnames); 
title('Extensibility (ζ) by Material');
xlabel('Material');
ylabel('Extensibility, ζ (%)');

for i = 1:6
    meanLabel = sprintf('Mean: %.2f', exmeans(i));
    stdLabel = sprintf('SD: %.2f', exstds(i));
    
    text(i, exmeans(i) + exstds(i) + 1.6, meanLabel, 'HorizontalAlignment', 'center');
    text(i, exmeans(i) + exstds(i) + 0.6, stdLabel, 'HorizontalAlignment', 'center');
end
hold off


%Ultimate Tensile Strength
thin_UT_instron = [ultimate_tensile_stress{1,:}]; 
thick_UT_instron = [ultimate_tensile_stress{2,:}]; 
thin_UT_minstron = [ultimate_tensile_stress{3,:}];
thick_UT_minstron = [ultimate_tensile_stress{4,:}];
raw_UT = [ultimate_tensile_stress{5,:}];
treated_UT = [ultimate_tensile_stress{6,:}];

UTmeans = [mean(thin_UT_instron); mean(thick_UT_instron); mean(thin_UT_minstron); mean(thick_UT_minstron); mean(raw_UT); mean(treated_UT)];
UTstds = [std(thin_UT_instron); std(thick_UT_instron); std(thin_UT_minstron); std(thick_UT_minstron); std(raw_UT); std(treated_UT)];

figure(8);
bar(UTmeans);
hold on

x = 1:6;
errorbar(x, UTmeans, UTstds, 'k', 'linestyle', 'none'); 

set(gca, 'XTickLabel', groupnames);  
title('Ultimate Tensile Strength (σUT) by Material');
xlabel('Material');
ylabel('Ultimate Tensile Strength, σUT (MPa)');

for i = 1:6
    meanLabel = sprintf('Mean: %.2f', UTmeans(i));
    stdLabel = sprintf('SD: %.2f', UTstds(i));
    
    text(i, UTmeans(i) + UTstds(i) + 0.4, meanLabel, 'HorizontalAlignment', 'center');
    text(i, UTmeans(i) + UTstds(i) + 0.2, stdLabel, 'HorizontalAlignment', 'center');
end
hold off

%% Section 2 Methods: 2.1 Sample Preparation

%calculating the sample dimensions for all teams data
% Data for raw skin
raw_skin_data = [
    35.89, 48, 45, 41, 30.8, 42, 49, 47.07, 50, 35.38, 45;   % Original Length (mm)
    13.03, 7, 10, 10.03, 8.6, 13.63, 10, 8.49, 7.59, 5.689, 10.36;   % Width (mm)
    0.27, 0.64, 0.55, 1.02, 0.595, 0.7, 0.62, 0.9, 0.39, 0.8636, 0.776;   % Thickness (mm)
    3.5181, 4.48, 5.5, 10.2306, 5.117, 9.541, 6.2, 7.641, 2.9601, 4.9130204, 8.03936  % Cross-sectional Area (mm^2)
];

% Data for treated skin
treated_skin_data = [
    32.76, 47.5, 46, 32, 29.5, 28, 48, 36.15, 41.5, 27.33, 27;   % Original Length (mm)
    8.94, 6, 4.81, 7.37, 10.9, 14.3, 9, 4.28, 5.65, 4.37, 7.14;   % Width (mm)
    0.22, 0.55, 0.65, 0.67, 0.83, 1.2, 0.46, 0.63, 0.39, 0.69, 1.08;   % Thickness (mm)
    1.9668, 3.3, 3.1265, 4.9379, 9.047, 17.16, 4.14, 2.6964, 2.2035, 3.0153, 7.7112  % Cross-sectional Area (mm^2)
];

% Data for silicone thin Instron
silicon_thin_instron_data = [
    42.79, 47, 37, 38, 57.7, 40, 35, 37, 35, 40, 36;   % Original Length (mm)
    12.7, 12, 12.21, 7, 41.06, 14, 12, 12.96, 12.45, 13.5, 11.4554;   % Width (mm)
    1.83, 1.85, 1.87, 2, 1.81, 2, 1.88, 1.74, 1.85, 2, 1.7;   % Thickness (mm)
    23.241, 22.2, 22.8327, 24.1425, 74.3186, 28, 22.56, 22.5504, 23.0325, 27, 19.47418  % Cross-sectional Area (mm^2)
];

% Data for silicone thin Minstron
silicon_thin_minstron_data = [
    54.2, 55, 52, 57, 57.7, 43, 39, 36.84, 35, 50, 36;   % Original Length (mm)
    12.7, 12.4, 12.21, 13.05, 41.06, 14, 12, 13.27, 12.45, 13.5, 11.4554;   % Width (mm)
    1.83, 1.82, 1.87, 1.85, 1.81, 2, 1.88, 1.85, 1.85, 2, 1.7;   % Thickness (mm)
    23.241, 22.568, 22.8327, 24.1425, 74.3186, 28, 22.6, 24.5495, 23.0325, 27, 19.47418  % Cross-sectional Area (mm^2)
];

% Data for silicone thick Instron
silicon_thick_instron_data = [
    38.63, 41, 43, 35, 54.11, 39, 36, 36, 37, 40, 44;   % Original Length (mm)
    11.61, 12, 13.53, 13, 40.24, 13, 13, 11.41, 12.31, 14, 11.56;   % Width (mm)
    3.29, 3.21, 3.27, 3, 3.27, 3.13, 3.28, 3.07, 3.29, 3.35, 3.12;   % Thickness (mm)
    38.1969, 38.52, 44.2431, 39, 131.5848, 40.69, 42.64, 35.0287, 40.4999, 46.9, 36.0672  % Cross-sectional Area (mm^2)
];

% Data for silicone thick Minstron
silicon_thick_minstron_data = [
    55.05, 57, 53, 59, 54.11, 42, 37, 39.35, 37, 40, 44;   % Original Length (mm)
    11.61, 14.3, 13.53, 14.7, 40.24, 13, 12, 12.45, 12.31, 14, 11.56;   % Width (mm)
    3.29, 3.28, 3.27, 3.22, 3.27, 3.13, 3.28, 3.15, 3.29, 3.35, 3.12;   % Thickness (mm)
    38.1969, 46.904, 44.2431, 47.334, 131.5848, 40.69, 39.4, 39.2175, 40.4999, 46.9, 36.0672  % Cross-sectional Area (mm^2)
];

% Define the headers for display
headers = {'Original Length (mm)', 'Width (mm)', 'Thickness (mm)', 'Cross-sectional Area (mm^2)'};

% Create a function to calculate mean ± std for each material
function result = calculate_mean_std(data)
    data_mean = mean(data, 2);
    data_std = std(data, 0, 2);
    result = strcat(num2str(data_mean, '%.2f'), ' ± ', num2str(data_std, '%.2f'));
end

% Calculate and display results for each material
disp('Raw Skin Results:');
disp(array2table(calculate_mean_std(raw_skin_data)', 'VariableNames', headers));

disp('Treated Skin Results:');
disp(array2table(calculate_mean_std(treated_skin_data)', 'VariableNames', headers));

disp('Silicone Thin Instron Results:');
disp(array2table(calculate_mean_std(silicon_thin_instron_data)', 'VariableNames', headers));

disp('Silicone Thin Minstron Results:');
disp(array2table(calculate_mean_std(silicon_thin_minstron_data)', 'VariableNames', headers));

disp('Silicone Thick Instron Results:');
disp(array2table(calculate_mean_std(silicon_thick_instron_data)', 'VariableNames', headers));

disp('Silicone Thick Minstron Results:');
disp(array2table(calculate_mean_std(silicon_thick_minstron_data)', 'VariableNames', headers));