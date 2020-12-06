close all; clear all; clc;
%% merge borders
% Find the border files
motor_borders_list = dir('../motor/motor_borders/area*')

motor_borders_list_cell = struct2cell(motor_borders_list);
motor_borders_list_cell = motor_borders_list_cell(1,:);

julich_borders_list = [motor_borders_list_cell];

% How many border files do we have ?
num_motor_areas = length(motor_borders_list);

% where are the border files stored?
path_to_motor_borders = '../motor/motor_borders/';

%% Let's put together the motor borders (drawn by Lucija)
% This just creates the motor.border file by putting together the first two regions
merge_motor = sprintf('wb_command -border-merge ../motor/motor.border -border  %s -border  %s ', fullfile([path_to_motor_borders motor_borders_list(1).name]), fullfile([path_to_motor_borders motor_borders_list(2).name]));
[u v]=system(merge_motor);

% This adds on each of the other regions to the motor.border file one at a time. 
for current_border = 3:num_motor_areas
    sprintf('adding %s',motor_borders_list(current_border).name)
    merge_motor = sprintf('wb_command -border-merge ../motor/motor.border -border  ../motor/motor.border -border  %s ', fullfile([path_to_motor_borders motor_borders_list(current_border).name]));
    [u v]=system(merge_motor);

end

%% change motor borders to ROIs
!echo "filling in borders to create ROIs"
!wb_command -border-to-rois ../MacaqueYerkes19_v1.2_976nz/MNINonLinear/fsaverage_LR32k/MacaqueYerkes19_v1.2.L.inflated.32k_fs_LR.surf.gii ../motor/motor.border E:/workbench/test/motor/motor_3d.func.gii

% Let's change the mask to a ROI
!echo "creating a ROI from the mask border file "
!wb_command -border-to-rois ../MacaqueYerkes19_v1.2_976nz/MNINonLinear/fsaverage_LR32k/MacaqueYerkes19_v1.2.L.inflated.32k_fs_LR.surf.gii ../motor/motor_mask.border ../motor/non_motor_mask.func.gii -inverse

%% load in 3d ROI file
full_cortex_ROIs = gifti('../motor/motor_3d.func.gii');

%% identify and remove overlaps
overlapping_vertices = find(sum(full_cortex_ROIs.cdata,2)>1);
num_overlapping_vertices = length(overlapping_vertices);

for current_vertex = 1:num_overlapping_vertices
    
    overlapping_ROIs = find(full_cortex_ROIs.cdata(overlapping_vertices(current_vertex),:));
    num_overlapping_ROIs = length(overlapping_ROIs);
    max_overlapping_ROIs = max(overlapping_ROIs);
    full_cortex_ROIs.cdata(overlapping_vertices(current_vertex),overlapping_ROIs) = 0;
    full_cortex_ROIs.cdata(overlapping_vertices(current_vertex),max_overlapping_ROIs) = 1;

end

%% merge 3d file into a single 2d surface

num_ROIs = size(full_cortex_ROIs.cdata,2);
num_vertices = size(full_cortex_ROIs.cdata,1);
ROI_idx = 1:num_ROIs;
ROI_idx_mat = repmat(ROI_idx,num_vertices,1);
full_cortex_ROIs_map_3d = full_cortex_ROIs.cdata.*ROI_idx_mat;
full_cortex_ROIs_map_2d = sum(full_cortex_ROIs_map_3d,2);

%% write new gifti file

temp = gifti('../motor/macaque_myelin.func.gii');
full_cortex_ROIs_2d = temp;
full_cortex_ROIs_2d.cdata = full_cortex_ROIs_map_2d;

save(full_cortex_ROIs_2d,'../motor/motor_2d.func.gifti','Base64Binary');
%% plot 
% load in inflated left hemisphere surface
l_inflated = gifti('../MacaqueYerkes19_v1.2_976nz/MNINonLinear/fsaverage_LR32k/MacaqueYerkes19_v1.2.L.inflated.32k_fs_LR.surf.gii');

figure
plot(l_inflated,full_cortex_ROIs_2d);

%% Fill in the gaps in the motor atlas 

% Identify vertices that have a value 0 but lie inside the motor mask.
% This should identify border vertices between the brain sections.

% Find all the vertices without a parcel number - this should be the
% unfilled border areas and unsampled cortex.
border_plus_unsampled = find(full_cortex_ROIs_2d.cdata == 0);

% Let's load in the frontoparietooccipital mask
mask = gifti('../motor/non_motor_mask.func.gii');
mask_ind = find(mask.cdata);

% Find all the vertices without a parcel number excluding the unsampled cortex.
julich_border_ind = setdiff(border_plus_unsampled,mask_ind);

julich_borders = temp;
julich_borders.cdata = zeros(num_vertices,1);
julich_borders.cdata(julich_border_ind,1) = 1;

% Visualise borders
close all;
figure
plot(l_inflated,julich_borders);

figure
plot(l_inflated,mask);
% save(julich_borders,'../motor/julich_borders.func.gii','Base64Binary');
%%

% Let's find the distances for each border vertex to the other vertices

l_midthickness = gifti('../MacaqueYerkes19_v1.2_976nz/MNINonLinear/fsaverage_LR32k/MacaqueYerkes19_v1.2.L.midthickness.32k_fs_LR.surf.gii');

num_border_vertices = length(julich_border_ind);

full_cortex_ROIs_2d_filled = full_cortex_ROIs_2d;

for current_vertex = 1:num_border_vertices
   sprintf('Filling in border vertices - %f percent complete', 100*round(current_vertex/num_border_vertices,2))
   % assign each vertex to its closest sample    
   command = ['wb_command -surface-geodesic-distance ../MacaqueYerkes19_v1.2_976nz/MNINonLinear/fsaverage_LR32k/MacaqueYerkes19_v1.2.L.midthickness.32k_fs_LR.surf.gii ' sprintf('%d',julich_border_ind(current_vertex)-1) ' geo_dist_current_vertex.func.gii'];
   system(command);
   geo_dist_current_vertex = gifti('geo_dist_current_vertex.func.gii');

   [sorted_dist,dist_ind] = sort(geo_dist_current_vertex.cdata);

   % find the closest vertex with a non-zero parcel value
   closest_nonzero_vertex = dist_ind(min(find(full_cortex_ROIs_2d.cdata(dist_ind))));
   % now let's add this border vertex to the closest parcel
   full_cortex_ROIs_2d_filled.cdata(julich_border_ind(current_vertex),1) = full_cortex_ROIs_2d.cdata(closest_nonzero_vertex,1);

end

% save 
save(full_cortex_ROIs_2d_filled,'../julich_sections_2d_filled.func.gii','Base64Binary');

%% Now let's place receptor density in the appropriate parcels (AMPA)
load AMPA_density.mat
AMPA_density_surf = temp;
AMPA_density_surf.cdata = zeros(size(AMPA_density_surf.cdata));

for current_area = 1:num_motor_areas;
    
    sprintf('filling in AMPA density for %s',motor_data_order{current_area});
    vertices_in_area = find(full_cortex_ROIs_2d_filled.cdata==current_area);
    AMPA_density_surf.cdata(vertices_in_area) = AMPA_density(current_area);    
    
end

save(AMPA_density_surf,'../motor/AMPA_density.func.gii','Base64Binary');

%% Now let's place receptor density in the appropriate parcels (kain)
load kain_density.mat
kain_density_surf = temp;
kain_density_surf.cdata = zeros(size(kain_density_surf.cdata));

for current_area = 1:num_motor_areas;
    
    sprintf('filling in kain density for %s',motor_data_order{current_area});
    vertices_in_area = find(full_cortex_ROIs_2d_filled.cdata==current_area);
    kain_density_surf.cdata(vertices_in_area) = kain_density(current_area);    
    
end

save(kain_density_surf,'../motor/kain_density.func.gii','Base64Binary');

%% Now let's place receptor density in the appropriate parcels (NMDA)
load NMDA_density.mat
NMDA_density_surf = temp;
NMDA_density_surf.cdata = zeros(size(NMDA_density_surf.cdata));

for current_area = 1:num_motor_areas;
    
    sprintf('filling in NMDA density for %s',motor_data_order{current_area});
    vertices_in_area = find(full_cortex_ROIs_2d_filled.cdata==current_area);
    NMDA_density_surf.cdata(vertices_in_area) = NMDA_density(current_area);    
    
end

save(NMDA_density_surf,'../motor/NMDA_density.func.gii','Base64Binary');

%% Now let's place receptor density in the appropriate parcels (GABAa)
load GABAa_density.mat
GABAa_density_surf = temp;
GABAa_density_surf.cdata = zeros(size(GABAa_density_surf.cdata));

for current_area = 1:num_motor_areas;
    
    sprintf('filling in GABAa density for %s',motor_data_order{current_area});
    vertices_in_area = find(full_cortex_ROIs_2d_filled.cdata==current_area);
    GABAa_density_surf.cdata(vertices_in_area) = GABAa_density(current_area);    
    
end

save(GABAa_density_surf,'../motor/GABAa_density.func.gii','Base64Binary');

%% Now let's place receptor density in the appropriate parcels (GABAb)
load GABAb_density.mat
GABAb_density_surf = temp;
GABAb_density_surf.cdata = zeros(size(GABAb_density_surf.cdata));

for current_area = 1:num_motor_areas;
    
    sprintf('filling in GABAb density for %s',motor_data_order{current_area});
    vertices_in_area = find(full_cortex_ROIs_2d_filled.cdata==current_area);
    GABAb_density_surf.cdata(vertices_in_area) = GABAb_density(current_area);    
    
end

save(GABAb_density_surf,'../motor/GABAb_density.func.gii','Base64Binary');

%% Now let's place receptor density in the appropriate parcels (BZ)
load BZ_density.mat
BZ_density_surf = temp;
BZ_density_surf.cdata = zeros(size(BZ_density_surf.cdata));

for current_area = 1:num_motor_areas;
    
    sprintf('filling in BZ density for %s',motor_data_order{current_area});
    vertices_in_area = find(full_cortex_ROIs_2d_filled.cdata==current_area);
    BZ_density_surf.cdata(vertices_in_area) = BZ_density(current_area);    
    
end

save(BZ_density_surf,'..t/motor/BZ_density.func.gii','Base64Binary');

%% Now let's place receptor density in the appropriate parcels (M1)
load M1_density.mat
M1_density_surf = temp;
M1_density_surf.cdata = zeros(size(M1_density_surf.cdata));

for current_area = 1:num_motor_areas;
    
    sprintf('filling in M1 density for %s',motor_data_order{current_area});
    vertices_in_area = find(full_cortex_ROIs_2d_filled.cdata==current_area);
    M1_density_surf.cdata(vertices_in_area) = M1_density(current_area);    
    
end

save(M1_density_surf,'../motor/M1_density.func.gii','Base64Binary');

%% Now let's place receptor density in the appropriate parcels (M2)
load M2_density.mat
M2_density_surf = temp;
M2_density_surf.cdata = zeros(size(M2_density_surf.cdata));

for current_area = 1:num_motor_areas;
    
    sprintf('filling in M2 density for %s',motor_data_order{current_area});
    vertices_in_area = find(full_cortex_ROIs_2d_filled.cdata==current_area);
    M2_density_surf.cdata(vertices_in_area) = M2_density(current_area);    
    
end

save(M2_density_surf,'../motor/M2_density.func.gii','Base64Binary');

%% Now let's place receptor density in the appropriate parcels (M3)
load M3_density.mat
M3_density_surf = temp;
M3_density_surf.cdata = zeros(size(M3_density_surf.cdata));

for current_area = 1:num_motor_areas;
    
    sprintf('filling in M3 density for %s',motor_data_order{current_area});
    vertices_in_area = find(full_cortex_ROIs_2d_filled.cdata==current_area);
    M3_density_surf.cdata(vertices_in_area) = M3_density(current_area);    
    
end

save(M3_density_surf,'../motor/M3_density.func.gii','Base64Binary');

%% Now let's place receptor density in the appropriate parcels (alpha1)
load alpha1_density.mat
alpha1_density_surf = temp;
alpha1_density_surf.cdata = zeros(size(alpha1_density_surf.cdata));

for current_area = 1:num_motor_areas;
    
    sprintf('filling in alpha1 density for %s',motor_data_order{current_area});
    vertices_in_area = find(full_cortex_ROIs_2d_filled.cdata==current_area);
    alpha1_density_surf.cdata(vertices_in_area) = alpha1_density(current_area);    
    
end

save(alpha1_density_surf,'../motor/alpha1_density.func.gii','Base64Binary');

%% Now let's place receptor density in the appropriate parcels (alpha2)
load alpha2_density.mat
alpha2_density_surf = temp;
alpha2_density_surf.cdata = zeros(size(alpha2_density_surf.cdata));

for current_area = 1:num_motor_areas;
    
    sprintf('filling in alpha2 density for %s',motor_data_order{current_area});
    vertices_in_area = find(full_cortex_ROIs_2d_filled.cdata==current_area);
    alpha2_density_surf.cdata(vertices_in_area) = alpha2_density(current_area);    
    
end

save(alpha2_density_surf,'../motor/alpha2_density.func.gii','Base64Binary');

%% Now let's place receptor density in the appropriate parcels (ht1a)
load ht1a_density.mat
ht1a_density_surf = temp;
ht1a_density_surf.cdata = zeros(size(ht1a_density_surf.cdata));

for current_area = 1:num_motor_areas;
    
    sprintf('filling in ht1a density for %s',motor_data_order{current_area});
    vertices_in_area = find(full_cortex_ROIs_2d_filled.cdata==current_area);
    ht1a_density_surf.cdata(vertices_in_area) = ht1a_density(current_area);    
    
end

save(ht1a_density_surf,'../motor/ht1a_density.func.gii','Base64Binary');

%% Now let's place receptor density in the appropriate parcels (ht2)
load ht2_density.mat
ht2_density_surf = temp;
ht2_density_surf.cdata = zeros(size(ht2_density_surf.cdata));

for current_area = 1:num_motor_areas;
    
    sprintf('filling in ht2 density for %s',motor_data_order{current_area});
    vertices_in_area = find(full_cortex_ROIs_2d_filled.cdata==current_area);
    ht2_density_surf.cdata(vertices_in_area) = ht2_density(current_area);    
    
end

save(ht2_density_surf,'../motor/ht2_density.func.gii','Base64Binary');
