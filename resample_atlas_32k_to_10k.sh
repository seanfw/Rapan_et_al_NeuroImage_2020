# echo "downsampling Lucija's atlas from 32k to 10k"
# # Downsample Lucija's map from 32k to 10k to match Ting's connectivity data
# wb_command -label-resample ../julich_sections_2d_filled.label.gii /Users/seanfw/Dropbox/SFW_Wang_projects/Active_projects/Receptor_gradients/walsh-receptor-gradients/CMI_gradients/Macaque_32k/MacaqueYerkes19.L.sphere.32k_fs_LR.surf.gii /Users/seanfw/Dropbox/SFW_Wang_projects/Active_projects/Receptor_gradients/walsh-receptor-gradients/CMI_gradients/Macaque_10k/MacaqueYerkes19.L.sphere.10k_fs_LR.surf.gii ADAP_BARY_AREA ../julich_sections_2d_filled_10k.label.gii -area-surfs /Users/seanfw/Dropbox/SFW_Wang_projects/Active_projects/Receptor_gradients/walsh-receptor-gradients/CMI_gradients/Macaque_32k/MacaqueYerkes19.L.midthickness.32k_fs_LR.surf.gii /Users/seanfw/Dropbox/SFW_Wang_projects/Active_projects/Receptor_gradients/walsh-receptor-gradients/CMI_gradients/Macaque_10k/MacaqueYerkes19.L.midthickness.10k_fs_LR.surf.gii

# echo "downsampling Lucija's borders from 32k to 10k"
# # Downsample Lucija's map from 32k to 10k to match Ting's connectivity data
# borders_list=$(ls ../motor_borders_32k/*.border)
# 
# # echo $borders_list
# 
# 
# for current_border in $borders_list
# do
#     current_border_short=$(echo $current_border | cut -d'/' -f 3)
#     echo "downsampling $current_border_short"
# 
#      wb_command -border-resample $current_border /Users/seanfw/Dropbox/SFW_Wang_projects/Active_projects/Receptor_gradients/walsh-receptor-gradients/CMI_gradients/Macaque_32k/MacaqueYerkes19.L.sphere.32k_fs_LR.surf.gii /Users/seanfw/Dropbox/SFW_Wang_projects/Active_projects/Receptor_gradients/walsh-receptor-gradients/CMI_gradients/Macaque_10k/MacaqueYerkes19.L.sphere.10k_fs_LR.surf.gii ../motor_borders_10k/$current_border_short
#     
#     
# done

echo "downsampling Kennedy's atlas from 32k to 10k"
# Downsample Kennedy's map from 32k to 10k to match Ting's connectivity data
wb_command -label-resample /Users/seanfw/Dropbox/SFW_Wang_projects/Research_areas/Gradients/Palomero_Gallagher/walsh-macaque-receptor-gradients/atlases/kennedy_atlas_91.label.gii /Users/seanfw/Dropbox/SFW_Wang_projects/Active_projects/Receptor_gradients/walsh-receptor-gradients/CMI_gradients/Macaque_32k/MacaqueYerkes19.L.sphere.32k_fs_LR.surf.gii /Users/seanfw/Dropbox/SFW_Wang_projects/Active_projects/Receptor_gradients/walsh-receptor-gradients/CMI_gradients/Macaque_10k/MacaqueYerkes19.L.sphere.10k_fs_LR.surf.gii ADAP_BARY_AREA ../kennedy_atlas_91_10k.label.gii -area-surfs /Users/seanfw/Dropbox/SFW_Wang_projects/Active_projects/Receptor_gradients/walsh-receptor-gradients/CMI_gradients/Macaque_32k/MacaqueYerkes19.L.midthickness.32k_fs_LR.surf.gii /Users/seanfw/Dropbox/SFW_Wang_projects/Active_projects/Receptor_gradients/walsh-receptor-gradients/CMI_gradients/Macaque_10k/MacaqueYerkes19.L.midthickness.10k_fs_LR.surf.gii

