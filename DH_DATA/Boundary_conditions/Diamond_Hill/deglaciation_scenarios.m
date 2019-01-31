% Use this script to create deglaciation scenarios for
% DH_deglaciation_scenarios.mat. It makes use of the SLM toolbox, which can
% be found here: https://www.mathworks.com/matlabcentral/fileexchange/24443-slm-shape-language-modeling
close all

%% First load up the geochronological data from Diamond Hill, and (for now) only 
%use the ones from the Darwin Glacier side
load(['/Users/trevorhillebrand/Documents/Antarctica/Darwin-Hatherton/',...
    'Modeling/Koutnik model/DH_code/DH_DATA/Geochronology data/cosmo_data.mat']);

erratics_names = {'14-HAT-012-DH';'14-HAT-014-DH';'14-HAT-015-DH';'14-HAT-016-DH';...
    '14-HAT-017-DH';'14-HAT-022-DH';'14-HAT-024-DH';'14-HAT-032-DH'};

erratics_ind = ismember(data.parsed_output.DH.erratics.SampleName, erratics_names);

erratics_age = data.parsed_output.DH.erratics.t10St(erratics_ind);
erratics_err = data.parsed_output.DH.erratics.ext10St(erratics_ind);
erratics_height = data.parsed_output.DH.erratics.HeightAboveIceMargin(erratics_ind);


bedrock_names = {'14-HAT-006-DH';'14-HAT-026-DH';'14-HAT-033-DH'};

bedrock_ind = ismember(data.parsed_output.DH.bedrock.SampleName, bedrock_names);

bedrock_age = data.parsed_output.DH.bedrock.t14St(bedrock_ind);
bedrock_err = data.parsed_output.DH.bedrock.ext14St(bedrock_ind);
bedrock_height = data.parsed_output.DH.bedrock.HeightAboveIceMargin(bedrock_ind);
%% ??????????????????????????????????????????????????????????????????????????

%  Here, you have two nodes to play with to create thinning curves. Unless
%  you tell it to do otherwise, the model will use the C-14 saturated
%  (14-HAT-006-DH) sample to define the LGM thickness. The two nodes should be placed in
%  between 006-DH and 026-DH to create a thinning curve that fits the data,
%  but explores the largely unknown space in between. If you only want one 
%  such point, set the other to [0, 0] These points, along
%  with the name of the thinning curve (e.g., Rapid9ka) are defined in
%  Deglaciation_scenarios.csv. The following cells loop through those
%  scenarios, construct the curves, and save them.

scenarios = readtable('Deglaciation_scenarios.csv');
save_on = 'off'; %set to 'on' if you want to automatically save every scenario
%% ??????????????????????????????????????????????????????????????????????????
for jj = 1:length(scenarios.ScenarioName)


%% add the points you want to use to define the thinning curve
pt_1_age = scenarios.StopAge(jj);
pt_1_height = scenarios.StopHeight(jj);

pt_2_age = scenarios.StartAge(jj);
pt_2_height = scenarios.StartHeight(jj);

pt_ages = [pt_1_age; pt_2_age];
pt_height = [pt_1_height; pt_2_height];

%%
%concatenate erratics, bedrock, and your new points
all_ages = [erratics_age; bedrock_age; pt_1_age; pt_2_age];
all_height = [erratics_height; bedrock_height; pt_1_height; pt_2_height];
% all_err = [erratics_err; bedrock_height];


figure(jj); clf
herrorbar(erratics_age, erratics_height, erratics_err, 'o')
hold on
herrorbar(bedrock_age, bedrock_height, bedrock_err, '<')
plot(pt_ages, pt_height, '*')

%%
%sort by age
[all_ages, I] = sort(all_ages);  %I is an index vector
all_height = all_height(I);
% all_err = all_err(I);

if all_ages(1) ~= 0
    all_ages = [0; all_ages];
    all_height = [0; all_height];
end
%make the time vector to define points on
time = 0:100:pt_2_age;

%% Try a piecewise linear fit
d_knot_linear = 1.5e3;
% 
% slm_linear = slmengine((all_ages(1:end-1)),all_height(1:end-1),'degree',...
%     'linear','knots',0:d_knot_linear:pt_2_age+d_knot_linear, 'plot','on');
% S_at_GL_linear = [slmeval(time, slm_linear), max(all_height)]; 
% 
% % if the 2 kyr window leads to a large (> 0.01) negative slopes, try a 1 kyr window
% while any(slmeval(time, slm_linear, 1) < -0.01) 
% d_knot_linear = 0.75* d_knot_linear;
%     slm_linear = slmengine((all_ages(1:end-1)),all_height(1:end-1),'degree',...
%     'linear','knots',0:d_knot_linear:pt_2_age+d_knot_linear, 'plot','on');
%     S_at_GL_linear = [slmeval(time, slm_linear), max(all_height)]; 
%     
% end

slm_linear = slmengine(all_ages([1:10, 12:end-1]),all_height([1:10, 12:end-1]),'degree',...
    'linear','knots',all_ages([1:10, 12:end-1]), 'plot','on');
S_at_GL_linear = [slmeval(time, slm_linear), max(all_height)]; 

%% fit a 2nd-degree polynomial, excluding the saturated sample (this doesn't always work well)
[P,S] = polyfit(all_ages(1:end-1), all_height(1:end-1),2);
d_knot_quad = 1e3;
S_at_GL_quad = [polyval(P, time), max(all_height)];

%% Or a cubic fit

d_knot_cubic = 2e3; %spacing between knots. Might have to be adjusted to avoid large gradients

slm_cubic = slmengine((all_ages(1:end-1)),all_height(1:end-1),'degree',...
    'cubic','knots',0:d_knot_cubic:pt_2_age + d_knot_cubic, 'plot','off');
S_at_GL_cubic = [slmeval(time, slm_cubic), max(all_height)];

%the while-loop changes the knot spacing to prevent any large gradients
while any(slmeval(time, slm_cubic, 1) < -0.02) 

d_knot_cubic = 1.25*d_knot_cubic;

slm_cubic = slmengine((all_ages(1:end-1)),all_height(1:end-1),'degree',...
    'cubic','knots',0:d_knot_cubic:pt_2_age+d_knot_cubic, 'plot','off');
S_at_GL_cubic = [slmeval(time, slm_cubic), max(all_height)];
end

%% or a step function
d_knot_step = 1e3;

slm_step = slmengine((all_ages(1:end-1)),all_height(1:end-1),'degree',...
    'constant','knots',[all_ages([1:5,7:end-1]); pt_2_age+d_knot_step], 'plot','off');
S_at_GL_step = [slmeval(time, slm_step), max(all_height)];


%% Now plot these up
figure(jj)
hold on
plot([time, 30e3], S_at_GL_linear, 'k', 'linewidth', 1.5);
plot([time, 30e3], S_at_GL_cubic, 'k:', 'linewidth', 3);
plot([time, 30e3], S_at_GL_step, 'k--', 'linewidth', 1.5);
plot([time, 30e3], S_at_GL_quad, 'k-.', 'linewidth', 1.5);

xlim([0 2e4]); ylim([-50 500])
title(scenarios.ScenarioName{jj})
% and store them
eval([(scenarios.ScenarioName{jj}), '.linear = S_at_GL_linear;'])
eval([(scenarios.ScenarioName{jj}), '.cubic = S_at_GL_cubic;'])
eval([(scenarios.ScenarioName{jj}), '.step = S_at_GL_step;'])
eval([(scenarios.ScenarioName{jj}), '.quad = S_at_GL_quad;'])
eval([(scenarios.ScenarioName{jj}), '.time = -[time, 30e3];'])
end

%% save them if you want
if strcmp(save_on, 'on') && strfind(pwd, 'DH_DATA/Boundary_conditions/Diamond_Hill')
  
    for jj = 1:length(scenarios.ScenarioName)
        
        eval(['save ', (scenarios.ScenarioName{jj}), '.mat ', (scenarios.ScenarioName{jj}),';'])
    
    end
end

