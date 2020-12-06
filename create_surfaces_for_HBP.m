 
% load in Julich atlas
julich_atlas = gifti('submission_to_BALSA/scene/julich_sections_2d_filled.label.gii');

% load in inflated left hemisphere
lh_inflated = gifti('/Users/seanfw/Dropbox/SFW_Wang_projects/Active_projects/Receptor_gradients/walsh-receptor-gradients/MacaqueYerkes19_v1.2_976nz/MNINonLinear/fsaverage_LR32k/MacaqueYerkes19_v1.2.L.inflated.32k_fs_LR.surf.gii');
 
%% Get each receptor name
receptor_list_vlong = dir('submission_to_BALSA/scene/receptor_density/*.gii');
receptor_list_long = {receptor_list_vlong.name}';
receptor_list = cellfun(@(x) x(1:end-17), receptor_list_long, 'un', 0)
clear receptor_list_vlong receptor_list_long ;

%% import Julich area names
fid = fopen('motor/motor_data_order.txt');
Julich_area_list = textscan(fid,'%s%s%s');
fclose(fid);
Julich_area_list = Julich_area_list{1};
Julich_area_list = strrep(Julich_area_list,'area_','');
clear fid



%% Create individual surface files for each brain area in the atlas

num_regions = length(Julich_area_list);
num_receptors = length(receptor_list);
template_surface = gifti('motor/macaque_myelin.func.gii');
template_surface.cdata = zeros(length(template_surface.cdata),1);

mkdir submission_to_HBP/atlas_func_files/
mkdir submission_to_HBP/atlas_label_tables/
mkdir submission_to_HBP/atlas_label_surfaces/

for current_region = 1:num_regions
    
    % create the gifti 'func' file that holds the data for the atlas
    template_surface.cdata = zeros(length(template_surface.cdata),1);
    sprintf('creating the surface file for area %s',Julich_area_list{current_region})
    vertices_in_network = find(julich_atlas.cdata==current_region); % careful - data goes from 0-16, names in the label are listed from 1-17
    
    template_surface.cdata(vertices_in_network) = 1;
    
    save(template_surface,sprintf('submission_to_HBP/atlas_func_files/area_%s.func.gii',Julich_area_list{current_region}),'Base64Binary');

    % create the table with the label name information
    current_label_name = julich_atlas.labels.name(current_region+1);
    current_label_colour = round(julich_atlas.labels.rgba(current_region+1,:)*255);

    fileID = fopen(sprintf('submission_to_HBP/atlas_label_tables/area_%s.txt',Julich_area_list{current_region}),'w');
    fprintf(fileID,'%s \n',sprintf('area_%s',Julich_area_list{current_region}));
    fprintf(fileID,'%d %d %d %d %d\n',[1 current_label_colour]);
    fclose(fileID);
    
    
    % create the label surface file 
    create_label_file = sprintf('wb_command -metric-label-import submission_to_HBP/atlas_func_files/area_%s.func.gii submission_to_HBP/atlas_label_tables/area_%s.txt submission_to_HBP/atlas_label_surfaces/area_%s.label.gii' ,Julich_area_list{current_region},Julich_area_list{current_region},Julich_area_list{current_region});
    [u v]=system(create_label_file);

end


%% Create individual surface files for each receptor in each brain area in the atlas

num_regions = length(Julich_area_list);
num_receptors = length(receptor_list);
template_surface = gifti('motor/macaque_myelin.func.gii');
template_surface.cdata = zeros(length(template_surface.cdata),1);

for current_receptor = 1:num_receptors
     
    mkdir(sprintf('submission_to_HBP/%s_receptor_surfaces',receptor_list{current_receptor}))

    for current_region = 1:num_regions
        
        current_receptor_surface = gifti(sprintf('submission_to_BALSA/scene/receptor_density/%s_density.func.gii',receptor_list{current_receptor}));
        template_surface.cdata = zeros(length(template_surface.cdata),1);
        sprintf('creating the surface file for receptor %s in area %s ',receptor_list{current_receptor} ,Julich_area_list{current_region})
        vertices_in_network = find(julich_atlas.cdata==current_region); % careful - data goes from 0-16, names in the label are listed from 1-17

        template_surface.cdata(vertices_in_network) = current_receptor_surface.cdata(vertices_in_network);

        save(template_surface,sprintf('submission_to_HBP/%s_receptor_surfaces/%s_area_%s.func.gii',receptor_list{current_receptor},receptor_list{current_receptor},Julich_area_list{current_region}),'Base64Binary');


    end

end

%% Remove temporary files 

rmdir('submission_to_HBP/atlas_func_files/','s')

rmdir('submission_to_HBP/atlas_label_tables/','s')

