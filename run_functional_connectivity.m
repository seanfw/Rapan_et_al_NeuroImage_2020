% % Isolate left hemisphere myelin data
% !wb_command -cifti-separate /Users/seanfw/Dropbox/PRIME-DE_Oxford/data/cifti.func_pp_sm3.dtseries.nii COLUMN -metric CORTEX_LEFT /Users/seanfw/Dropbox/PRIME-DE_Oxford/data/oxford_lh_resting.func.gii
% % Isolate right hemisphere myelin data
% !wb_command -cifti-separate /Users/seanfw/Dropbox/PRIME-DE_Oxford/data/cifti.func_pp_sm3.dtseries.nii COLUMN -metric CORTEX_RIGHT /Users/seanfw/Dropbox/PRIME-DE_Oxford/data/oxford_rh_resting.func.gii

%% load in the data
func_lh = gifti('/Users/seanfw/Dropbox/PRIME-DE_Oxford/data/oxford_lh_resting.func.gii');
julich_motor_atlas = gifti('motor/julich_sections_2d_filled_10k.label.gii');
lh_inflated_10k = gifti('/Users/seanfw/Dropbox/PRIME-DE_Oxford/cifti_10k_1.5mm/MacaqueYerkes19.L.inflated.10k_fs_LR.surf.gii')
kennedy_atlas_91 = gifti('motor/kennedy_atlas_91_10k.label.gii');
areaList_Donahue = kennedy_atlas_91.labels.name(2:end)';

%%
num_regions = max(julich_motor_atlas.cdata);
num_vertices = size(func_lh.cdata,1);
num_timepoints_all = size(func_lh.cdata,2);
num_monkeys = 19;
num_timepoints_per_monkey = num_timepoints_all/num_monkeys;

%% Calculate the principal component timecourses and functional connectivity maps for each motor/premotor area

mean_fc_all_regions = nan(num_vertices,num_regions);
pc1_fc_all_regions = nan(num_vertices,num_regions);
pc1_timeseries_all_regions = nan(num_timepoints_all,num_regions);

for current_area = 1:num_regions
    vertices_in_parcel = find(julich_motor_atlas.cdata==current_area);
    vertices_in_parcel_timeseries = func_lh.cdata(vertices_in_parcel,:);
    vertices_in_parcel_timeseries_demeaned = nan(size(vertices_in_parcel_timeseries));
    for current_monkey = 1:num_monkeys
        % demean current monkeys data (seems like it is *almost* demeaned
        % already, so shouldn't make a big difference
        current_monkey_timepoints = (current_monkey-1)*num_timepoints_per_monkey + (1:num_timepoints_per_monkey);
        vertices_in_parcel_timeseries_demeaned(:,current_monkey_timepoints) = vertices_in_parcel_timeseries(:,current_monkey_timepoints) - repmat(mean(vertices_in_parcel_timeseries(:,current_monkey_timepoints),2),1,num_timepoints_per_monkey);

    end
    
    % Calculate functional connectivity for mean ROI timeseries
    mean_roi_timeseries = mean(func_lh.cdata(vertices_in_parcel,:));
%     mean_fc_all_regions(:,current_area) = corr(func_lh.cdata',mean_roi_timeseries');
   
    [coeff,score,~,~,explained] = pca(vertices_in_parcel_timeseries_demeaned');
    pc1_timeseries_all_regions(:,current_area) = score(:,1);
    pc1_fc_all_regions(:,current_area) = corr(func_lh.cdata',pc1_timeseries_all_regions(:,current_area));
    sprintf('variance explained by PC1 - region %d - %.02f', current_area , explained(1))


end

% mean_fc_all_regions_z = atanh(mean_fc_all_regions);
pc1_fc_all_regions_z = atanh(pc1_fc_all_regions);
% 
% % Save functional connectivity as gifti
% mean_func_conn_motor = gifti(mean_fc_all_regions_z);
% save(mean_func_conn_motor,'mean_func_conn_motor.func.gifti','Base64Binary');
pc1_func_conn_motor = gifti(pc1_fc_all_regions_z);
save(pc1_func_conn_motor,'pc1_func_conn_motor.func.gifti','Base64Binary');

%% Calculate the principal component timecourses and functional connectivity maps for each area of the Kennedy atlas

num_regions_Kennedy = max(kennedy_atlas_91.cdata);

mean_fc_all_regions_Kennedy = nan(num_vertices,num_regions_Kennedy);
pc1_fc_all_regions_Kennedy = nan(num_vertices,num_regions_Kennedy);
pc1_timeseries_all_regions_Kennedy = nan(num_timepoints_all,num_regions);

vertices_in_motor_atlas = find(julich_motor_atlas.cdata);
percent_non_motor_vertices = nan(num_regions_Kennedy,1);


for current_area = 1:num_regions_Kennedy
    vertices_in_parcel = find(kennedy_atlas_91.cdata==current_area);
    non_motor_vertices_in_parcel = vertices_in_parcel(~ismember(vertices_in_parcel,vertices_in_motor_atlas));
    percent_non_motor_vertices(current_area) = 100*length(non_motor_vertices_in_parcel)/length(vertices_in_parcel);
    vertices_in_parcel_timeseries = func_lh.cdata(non_motor_vertices_in_parcel,:);
    vertices_in_parcel_timeseries_Kennedy_demeaned = nan(size(vertices_in_parcel_timeseries));

    for current_monkey = 1:num_monkeys
        % demean current monkeys data (seems like it is *almost* demeaned
        % already, so shouldn't make a big difference
        current_monkey_timepoints = (current_monkey-1)*num_timepoints_per_monkey + (1:num_timepoints_per_monkey);
        vertices_in_parcel_timeseries_Kennedy_demeaned(:,current_monkey_timepoints) = vertices_in_parcel_timeseries(:,current_monkey_timepoints) - repmat(mean(vertices_in_parcel_timeseries(:,current_monkey_timepoints),2),1,num_timepoints_per_monkey);

    end
    

    % Calculate functional connectivity for first principal component
    % timeseries of each area
    [coeff,score,~,~,explained] = pca(vertices_in_parcel_timeseries_Kennedy_demeaned');
    pc1_timeseries_all_regions_Kennedy(:,current_area) = score(:,1);
    sprintf('variance explained by PC1 - region %d - %.02f', current_area , explained(1))
%     pc1_fc_all_regions_Kennedy(:,current_area) = corr(func_lh.cdata',score(:,1));
end

[overlaps_sorted overlap_ind] = sort(percent_non_motor_vertices);

regions_ordered_by_overlap = areaList_Donahue(overlap_ind)
% mean_fc_all_regions_z_Kennedy = atanh(mean_fc_all_regions_Kennedy);
% pc1_fc_all_regions_z_Kennedy  = atanh(pc1_fc_all_regions_Kennedy);

%% Correlate Lucija's regions with Kennedy regions
Julich_Lyon_FC = corr(pc1_timeseries_all_regions,pc1_timeseries_all_regions_Kennedy);
Julich_Julich_FC = corr(pc1_timeseries_all_regions);
Lyon_Lyon_FC = corr(pc1_timeseries_all_regions_Kennedy);

%% Create connectivity fingerprints
Lyon_regions_of_interest = {'1','2','3','5','9','10','24c','32','45A','46d','46v','7A','7B','7m','8l','8m','8r','9/46d','9/46v','LIP','MIP','ProM','VIP'}
num_Lyon_regions_of_interest = length(Lyon_regions_of_interest);
Lyon_regions_of_interest_index = nan(num_Lyon_regions_of_interest,1);
for current_region = 1:num_Lyon_regions_of_interest
    Lyon_regions_of_interest_index(current_region) =  find(strcmp(Lyon_regions_of_interest{current_region},areaList_Donahue));
end

axlims = repmat([-0.1;0.9],1,num_Lyon_regions_of_interest);
Julich_Lyon_FC_ROIs = Julich_Lyon_FC(:,Lyon_regions_of_interest_index);
Julich_Lyon_FC_ROIs_z = atanh(Julich_Lyon_FC_ROIs);

%%
% import Julich area names
fid = fopen('motor/motor_data_order.txt');
Julich_area_list = textscan(fid,'%s%s%s');
fclose(fid);
Julich_area_list = Julich_area_list{1};
Julich_area_list = strrep(Julich_area_list,'area_','');


%% New for paper revision - investigate strength of connectivity between areas F2v, F3 & 4

ind_F2v = find(strcmp('F2v',Julich_area_list))
ind_F3 = find(strcmp('F3',Julich_area_list))
ind_4m = find(strcmp('4m',Julich_area_list))
ind_4a = find(strcmp('4a',Julich_area_list))
ind_4p = find(strcmp('4p',Julich_area_list))

ind_F1_Lyon = find(strcmp('F1',areaList_Donahue))
ind_F2_Lyon = find(strcmp('F2',areaList_Donahue))
ind_F3_Lyon = find(strcmp('F3',areaList_Donahue))
ind_F4_Lyon = find(strcmp('F4',areaList_Donahue))
ind_F5_Lyon = find(strcmp('F5',areaList_Donahue))
ind_F6_Lyon = find(strcmp('F6',areaList_Donahue))
ind_F7_Lyon = find(strcmp('F7',areaList_Donahue))


%
% Get the list of all the Julich areas followed by all the Lyon areas
Julich_Lyon_concatenated_area_list = [Julich_area_list; areaList_Donahue];
% Remove Lyon motor regions from area list
Julich_Lyon_concatenated_area_list(num_regions + ind_F1_Lyon: num_regions + ind_F7_Lyon) = []


% Concatenated the functional connectivity between Julich areas with the FC
% between Julich and Lyon areas. Now each row in a Julich area, and the
% columns are Julich (first 16) and Lyon (last 91) areas.
Julich_Julich_Lyon_FC_temp = [Julich_Julich_FC, Julich_Lyon_FC];
Julich_Julich_Lyon_FC = Julich_Julich_Lyon_FC_temp;
% Remove connectivity of Lyon motor regions (in practice, this is only a
% few voxels that lie outside the Julich motor atlas)
Julich_Julich_Lyon_FC(:,num_regions + ind_F1_Lyon: num_regions + ind_F7_Lyon) = [];


%%
conn_F2v = Julich_Julich_Lyon_FC(ind_F2v,:)
[conn_F2v_sorted,conn_F2v_sorted_ind] = sort(conn_F2v,'descend')

conn_F2v_sorted_area_list = Julich_Lyon_concatenated_area_list(conn_F2v_sorted_ind)

% visualise F2v connectivity with other areas;

fHand = figure('units','normalized','outerposition',[0 0 0.5 1]);
set(gcf,'color','w');
aHand = axes('parent', fHand);
hold(aHand, 'on')

for current_region = 2:length(conn_F2v_sorted)
    if current_region == find(strcmp('4a',conn_F2v_sorted_area_list)) || current_region == find(strcmp('4m',conn_F2v_sorted_area_list)) || current_region == find(strcmp('4p',conn_F2v_sorted_area_list));
         barh(current_region, conn_F2v_sorted(current_region),'r');
    else
        barh(current_region, conn_F2v_sorted(current_region),'k');
    end
end
set(gca,'Ydir','reverse')
%title('Best predictors of acute plasticity')
set(gca,'ytick',2:length(conn_F2v_sorted))
set(gca,'yticklabel',conn_F2v_sorted_area_list(2:end))
xlabel('Correlation Coefficient')
caxis([-0.1 0.5])
%c=colorbar;
%ylabel(c,{'hippocampal functional connectivity, pre-lesion';'(Pearsons r)'},'FontSize',20);
set(gca, 'FontSize', 10)
title('F2v functional connectivity');
saveas(fHand,'motor/functional_conn_images/F2v_conn_barh_chart.png');

%%
% repeat for F3
conn_F3 = Julich_Julich_Lyon_FC(ind_F3,:);
[conn_F3_sorted,conn_F3_sorted_ind] = sort(conn_F3,'descend')

conn_F3_sorted_area_list = Julich_Lyon_concatenated_area_list(conn_F3_sorted_ind);

% visualise F3 connectivity with other areas;

fHand = figure('units','normalized','outerposition',[0 0 0.5 1]);
set(gcf,'color','w');
aHand = axes('parent', fHand);
hold(aHand, 'on')

for current_region = 2:length(conn_F3_sorted)
    if current_region == find(strcmp('4a',conn_F3_sorted_area_list)) || current_region == find(strcmp('4m',conn_F3_sorted_area_list)) || current_region == find(strcmp('4p',conn_F3_sorted_area_list));
         barh(current_region, conn_F3_sorted(current_region),'r');
    else
        barh(current_region, conn_F3_sorted(current_region),'k');
    end
end
set(gca,'Ydir','reverse')
%title('Best predictors of acute plasticity')
set(gca,'ytick',2:length(conn_F3_sorted))
set(gca,'yticklabel',conn_F3_sorted_area_list(2:end))
xlabel('Correlation Coefficient')
caxis([-0.1 0.5])
%c=colorbar;
%ylabel(c,{'hippocampal functional connectivity, pre-lesion';'(Pearsons r)'},'FontSize',20);
set(gca, 'FontSize', 10)

title('F3 functional connectivity');


saveas(fHand,'motor/functional_conn_images/F3_conn_barh_chart.png');

%% Create revised surface plots with colorbrewer2 colormap - new for revision

% Colours from Colorbrewer2 - with extra yellow to highlight region of
% interest
eleven_class_RdBu_Ye = [103,0,31;178,24,43;214,96,77;244,165,130;253,219,199;247,247,247;209,229,240;146,197,222;67,147,195;33,102,172;0,0,0;]./255


customcmap = flipud(eleven_class_RdBu_Ye);

for current_region = 1:num_regions
    pc1_func_conn_motor_temp = pc1_func_conn_motor;
    vertices_in_parcel = find(julich_motor_atlas.cdata==current_region);
    % show seed region in black
    pc1_func_conn_motor_temp.cdata(vertices_in_parcel,current_region) = -1;
    % show midline in white
    pc1_func_conn_motor_temp.cdata(isnan(pc1_func_conn_motor_temp.cdata))=0;
    
    myfig = figure('units','normalized','outerposition',[0.5 0.4 0.5 0.6])
    set(gcf,'color','w');
    colormap(customcmap)
    mysurf = plot(lh_inflated_10k,pc1_func_conn_motor_temp,current_region);
    direction1 = [0 1 0];
    direction2 = [0 0 1];
    direction3 = [1 0 0];
    rotate(mysurf,direction1,60);
    rotate(mysurf,direction2,90);
    rotate(mysurf,direction3,15);
    material([0.5,0.5,0.15])
    h = colorbar;
    set(h,'location','south')
    h.FontSize = 20;
    caxis([-0.25,0.25])
    ylabel(h, 'functional connectivity (z)')
    saveas(myfig,sprintf('motor/functional_conn_images/func_conn_%s_lateral.png',Julich_area_list{current_region}));

    rotate(mysurf,direction1,180);
    rotate(mysurf,direction3,30);
    light('Position',[0 1 0])



    saveas(myfig,sprintf('motor/functional_conn_images/func_conn_%s_medial.png',Julich_area_list{current_region}));


end



%% plot

close all
figure
subplot(3,2,1)
spider_plot(Julich_Lyon_FC_ROIs_z(1:3,:),'AxesLabels',Lyon_regions_of_interest,'AxesLimits',axlims, 'AxesDisplay', 'one');
legend(Julich_area_list{1:3}, 'Location', 'northeastoutside');

subplot(3,2,2)
axlims2 = repmat([-0.1;0.4],1,num_Lyon_regions_of_interest);
spider_plot(Julich_Lyon_FC_ROIs_z([4:6,13],:),'AxesLabels',Lyon_regions_of_interest,'AxesLimits',axlims2, 'AxesDisplay', 'one');
legend(Julich_area_list{[4:6,13]}, 'Location', 'northeastoutside');

subplot(3,2,3)
axlims3 = repmat([-0.1;0.6],1,num_Lyon_regions_of_interest);
spider_plot(Julich_Lyon_FC_ROIs_z(7:9,:),'AxesLabels',Lyon_regions_of_interest,'AxesLimits',axlims3, 'AxesDisplay', 'one');
legend(Julich_area_list{7:9}, 'Location', 'northeastoutside');

subplot(3,2,4)
spider_plot(Julich_Lyon_FC_ROIs_z(10:12,:),'AxesLabels',Lyon_regions_of_interest,'AxesLimits',axlims, 'AxesDisplay', 'one');
legend(Julich_area_list{10:12}, 'Location', 'northeastoutside');

subplot(3,2,5)
spider_plot(Julich_Lyon_FC_ROIs_z(14:16,:),'AxesLabels',Lyon_regions_of_interest,'AxesLimits',axlims3, 'AxesDisplay', 'one');
legend(Julich_area_list{14:16}, 'Location', 'northeastoutside');

%% Save as Excel for Lucija
% Data for FC fingerprints
filename = 'motor/func_conn_data/connectivity_fingerprint_data.xlsx';
writematrix(Julich_Lyon_FC_ROIs_z,filename,'Sheet',1,'Range','B2:X17')
writecell(Julich_area_list,filename,'Sheet',1,'Range','A2');
writecell(Lyon_regions_of_interest,filename,'Sheet',1,'Range','B1');


%% Whole cortex data 
Julich_Julich_Lyon_FC_z = atanh(Julich_Julich_Lyon_FC);

filename = 'motor/func_conn_data/connectivity_data_whole_cortex_Julich_Lyon_r.xlsx';
writematrix(Julich_Julich_Lyon_FC,filename,'Sheet',1,'Range','B2:DD17')
writecell(Julich_area_list,filename,'Sheet',1,'Range','A2');
writecell(Julich_Lyon_concatenated_area_list',filename,'Sheet',1,'Range','B1');

filename = 'motor/func_conn_data/connectivity_data_whole_cortex_Julich_Lyon_z.xlsx';
writematrix(Julich_Julich_Lyon_FC_z,filename,'Sheet',1,'Range','B2:DD17')
writecell(Julich_area_list,filename,'Sheet',1,'Range','A2');
writecell(Julich_Lyon_concatenated_area_list',filename,'Sheet',1,'Range','B1');

save Julich_Julich_Lyon_FC_16x107.mat Julich_Julich_Lyon_FC Julich_Julich_Lyon_FC_z Julich_area_list Julich_Lyon_concatenated_area_list


%% plot 

close all;

figure
plot(lh_inflated_10k,mean_func_conn_motor,1)
colorbar

figure
plot(lh_inflated_10k,mean_func_conn_motor,2)
colorbar


%% Now create individual monkey/area func conn maps

Lyon_area_list = areaList_Donahue; 

Julich_Lyon_FC_indiv = nan(num_regions,num_regions_Kennedy,num_monkeys);

pc1_timeseries_indiv_monkeys = nan(num_timepoints_per_monkey,num_regions,num_monkeys);
pc1_timeseries_indiv_monkeys_Kennedy = nan(num_timepoints_per_monkey,num_regions_Kennedy,num_monkeys);

vertices_in_motor_atlas = find(julich_motor_atlas.cdata);
percent_non_motor_vertices = nan(num_regions_Kennedy,1);


for current_monkey = 1:num_monkeys
    
    current_monkey_timepoints = (1:num_timepoints_per_monkey) + num_timepoints_per_monkey*(current_monkey-1);
    
    % 16 Julich areas
    for current_area = 1:num_regions
        sprintf('monkey %d - area %d, %s',current_monkey,current_area, Julich_area_list{current_area})
        vertices_in_parcel = find(julich_motor_atlas.cdata==current_area);
        vertices_in_parcel_timeseries = func_lh.cdata(vertices_in_parcel,current_monkey_timepoints);

       [coeff,score,~,~,explained] = pca(vertices_in_parcel_timeseries');
       pc1_timeseries_indiv_monkeys(:,current_area,current_monkey) = score(:,1);


        % Filling in the maps for all the areas in monkey 1 first, then the
        % maps for all the areas in monkey 2 etc
    %     fc_indiv_monkey_regions(:,current_area + num_regions*(current_monkey-1)) = corr(func_lh.cdata(:,current_monkey_timepoints)',mean_roi_timeseries');
    end
    
    % 91 Lyon areas
    for current_area = 1:num_regions_Kennedy
        sprintf('monkey %d - Kennedy area %d, %s',current_monkey,current_area, areaList_Donahue{current_area})
        vertices_in_parcel = find(kennedy_atlas_91.cdata==current_area);
        non_motor_vertices_in_parcel = vertices_in_parcel(~ismember(vertices_in_parcel,vertices_in_motor_atlas));
        percent_non_motor_vertices(current_area) = 100*length(non_motor_vertices_in_parcel)/length(vertices_in_parcel);

        % only extract data from vertices not overlapping with the
        % motor/premotor atlas
        vertices_in_parcel_timeseries = func_lh.cdata(non_motor_vertices_in_parcel,current_monkey_timepoints);

       [coeff,score,~,~,explained] = pca(vertices_in_parcel_timeseries');
       pc1_timeseries_indiv_monkeys_Kennedy(:,current_area,current_monkey) = score(:,1);
    end
    
    Julich_Lyon_FC_indiv(:,:,current_monkey) = corr(pc1_timeseries_indiv_monkeys(:,:,current_monkey),pc1_timeseries_indiv_monkeys_Kennedy(:,:,current_monkey));
    
end



Julich_Lyon_FC_indiv_z = atanh(Julich_Lyon_FC_indiv);

% save Julich_Lyon_FC_indiv_z_16x91.mat Julich_Lyon_FC_indiv_z Julich_area_list Lyon_area_list
% 
% func_conn_indiv_motor = gifti(fc_indiv_monkey_regions_z);
% %%
% save(func_conn_indiv_motor,'motor/func_conn_data/func_conn_indiv_motor.func.gifti','Base64Binary');

%% Create connectivity fingerprints
Lyon_regions_of_interest = {'1','2','3','5','9','10','24c','32','45A','46d','46v','7A','7B','7m','8l','8m','8r','9/46d','9/46v','LIP','MIP','ProM','VIP'}
Lyon_regions_of_interest = {'F1','F2','F3','F4','F5','F6','F7'}

num_Lyon_regions_of_interest = length(Lyon_regions_of_interest);
Lyon_regions_of_interest_index = nan(num_Lyon_regions_of_interest,1);
for current_region = 1:num_Lyon_regions_of_interest
    Lyon_regions_of_interest_index(current_region) =  find(strcmp(Lyon_regions_of_interest{current_region},areaList_Donahue));
end

axlims = repmat([-0.1;0.9],1,num_Lyon_regions_of_interest);
Julich_Lyon_FC_indiv_ROIs_z = Julich_Lyon_FC_indiv_z(:,Lyon_regions_of_interest_index,:);

% save Julich_Lyon_FC_indiv_z_16x23.mat Julich_Lyon_FC_indiv_ROIs_z Lyon_regions_of_interest Julich_area_list

%% Create mask

% masking out the midline (which is all nans)
lh_cortex_mask_10k = double(~isnan(fc_indiv_monkey_regions_z(:,1)));

lh_cortex_mask_10k_gifti = gifti(lh_cortex_mask_10k);
save(lh_cortex_mask_10k_gifti,'motor/func_conn_data/lh_cortex_mask_10k.func.gifti','Base64Binary');

