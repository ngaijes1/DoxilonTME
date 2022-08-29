%% Histology_CellDeath_CellType_final.m Adapted from Benjamin Kingston by Jessica Ngai (Warren Chan lab)
matlab_folder = pwd;

%Select folder contining 3 channel images
Main_folder = uigetdir();
cd(Main_folder);

%% Define cell type and key terms for directories
files = dir('*MOD*');
cell_marker = 1; 
%Options are: 
% 1 = F480
% 3 = Ly6g

save_dir = strcat(Main_folder,'\Results with Cell Type F480 800 thresh'); %change cell type's bkg threshold

    if exist(save_dir, 'dir')~=7
        mkdir(save_dir);
    end

cd(Main_folder)
%% CHANGE FOR NAME OF SECONDARY ONLY FILES
sec_files = dir('*Secondary Only*'); 

if cell_marker == 1   %F480
    img_file = [dir('*F4_80-Cy5*.tiff');dir('*F4_80-Cy5_01*.tiff')]; 
        
else cell_marker == 3   %Ly6G 
    img_file = [dir('*Ly6G-Cy5*.tiff');dir('*Ly6G-Cy5_01*.tiff')];

end

for  i = 1:size(img_file,1)  %loop the size of array of img_file
tic                          
 cd(Main_folder)             
 
     [~,shortfile] = fileparts(img_file(i).name);   
display(['Analyzing ' shortfile])                   

whole_img = loadtiff(img_file(i).name);             

dapi = whole_img(:,:,1);                
cell_type = whole_img(:,:,3);           

dapi_brightspots = stretch(dapi)<253;
dapi_brightspots = uint16(dapi_brightspots);
dapi = dapi.*dapi_brightspots;

dapi_thresh_lvl = adaptthresh(dapi,0.4);         
dapi_thresh = imbinarize(dapi,dapi_thresh_lvl);  

%% Create the mask or edge of tissue using nuclei
tissue_boundary_thresh = dapi>700;      %accept only positions of dapi channel greater than 700 (intensity)
level_erode = strel('sphere',3);        
level_dilate = strel('sphere',3);       
tissue_boundary_erode = imerode(tissue_boundary_thresh,level_erode);    
tissue_boundary_dilate = imdilate(tissue_boundary_erode,level_dilate);  

dt_tissue = bwdist(tissue_boundary_dilate); 
dt_tissue_boundary = dt_tissue<10;          
dt_tissue_boundary_fill = imfill(dt_tissue_boundary,'holes');   
dt_tissue_boundary_fill_rmove_small = bwareaopen(dt_tissue_boundary_fill,10000);  %remove all connected components that have less than 10000 pix

tissue_boundary = dt_tissue_boundary_fill_rmove_small;  

dapi_thresh_crop = dapi_thresh.*tissue_boundary;        %element-wise multiplication of matrix to apply tissue boundary onto dapi matrix; 
%% Dilating the DAPI for nuclei mask
threshnuc = dapi_thresh_crop>0;         
    dt_threshnuc = dt(threshnuc);       
    seeds = maxima(dt_threshnuc,2,0);   
    seeds2 = dilation(seeds,4.5)>0;     
    image_out = waterseed(seeds2,max(dt_threshnuc)-dt_threshnuc,1,0,0); 
    threshnuc(image_out) = false;       

tissue_boundary_rmove_edge_dt = dt(tissue_boundary);            
tissue_boundary_rmove_edge = tissue_boundary_rmove_edge_dt>70;  %keep anything above 70 intensity, save as logical matrix
dapi_thresh_crop_noedge = tissue_boundary_rmove_edge.*threshnuc; 

nuclei_seg = uint8(dapi_thresh_crop_noedge)>0;      
SmallNucleiRemove = bwareaopen(nuclei_seg, 5);      
LargeNuclei = bwareaopen(SmallNucleiRemove, 275,4); 
Nuclei_cleaned = SmallNucleiRemove - LargeNuclei;
nuclei_seg = Nuclei_cleaned;

nuc_label = label(nuclei_seg>0); 
nuc_dilated_label = dip_growregions(nuc_label,[],[],1,1,'low_first'); 
nuc_label = double(nuc_dilated_label); 

%% Assigning cell type 
all_cells_stats = regionprops(nuc_label,cell_type,'MeanIntensity','PixelIdxList'); 
all_cells_stats_struct = all_cells_stats;

max_cell_in_tissue = max(nuc_label(:));   
min_cell_in_tissue = min(nuc_label(:))+1; 

new_all_nuc =  uint32(nuc_label);     

    for b = min_cell_in_tissue:max_cell_in_tissue
            new_all_nuc(all_cells_stats_struct(b).PixelIdxList) = all_cells_stats_struct(b).MeanIntensity;
    end
    
    all_nuc = new_all_nuc; 
    
%% Calculate background from secondary only controls
cd(Main_folder)
means_table_cell = zeros(2,2); 

for  j = 1:size(sec_files,1)
  
    secondary_img = loadtiff(sec_files(j).name);    
    secondary_img_cell = secondary_img(:,:,3);      
    crop_secondary_img_cell = secondary_img_cell>800;  %binary great than 800 for f480, 800 for ly6g
    crop_secondary_img_done_cell = secondary_img_cell.*uint16(crop_secondary_img_cell); 
    mean_background_cell = mean(crop_secondary_img_done_cell(:));   
    std_background_cell = std(double(crop_secondary_img_done_cell(:))); 
 
    means_table_cell(j,1) = mean_background_cell;  
    means_table_cell(j,2) = std_background_cell;   
        
end
mean_background_cell = mean(means_table_cell(:,1)); 
std_background_cell = mean(means_table_cell(:,2));  

%% Cell Type Postive cut-offs
background_sing_thresh_cell_3std = mean_background_cell+(3*std_background_cell); 
background_sing_thresh_cell_2std = mean_background_cell+(2*std_background_cell);

%WITHIN 2 STDEV
test_cell_std2 = all_nuc>background_sing_thresh_cell_2std; %compare if nuc intensity is > than 2 std threshold bkg; logical array
label_test_cell_2std = bwlabel(test_cell_std2);            
num_mac_cell_2std = max(label_test_cell_2std(:));          
num_total = max(nuc_label(:));                   

percent_cell_2std = num_mac_cell_2std/num_total*100;  

%WITHIN 3 STDEV
test_cell_std3 = all_nuc>background_sing_thresh_cell_3std; %compare if nuc intensity is > than 3 std threshold bkg; logical array 
    label_test_cell_3std = bwlabel(test_cell_std3);        
num_mac_cell_3std = max(label_test_cell_3std(:));          
num_total = max(nuc_label(:));                  

percent_cell_3std = num_mac_cell_3std/num_total*100;  


%% Combine Results Tables 

summary_stats = [background_sing_thresh_cell_2std percent_cell_2std background_sing_thresh_cell_3std percent_cell_3std num_mac_cell_2std num_mac_cell_3std num_total];
sumary_stats_table = array2table(summary_stats);
sumary_stats_table.Properties.VariableNames = {'TWOSTD_Signal_Cutoff' 'TWOSTD_Percent_cell_type_positive' 'THREESTD_Signal_Cut_off' 'THREESTD_Percent_cell_type_positive' 'NUM_CELLS_CELLTYPE_twostd' 'NUM_CELLS_CELLTYPE_threestd' 'TOTAL_CELL_NUMBER'};

%% Write Results
cd(save_dir)

summary_stats_table_name = strcat(shortfile,'_Results_Summary_stats.csv');
writetable(sumary_stats_table,summary_stats_table_name);

Cell_type_int_name = strcat(shortfile,'_cell_type_intensity.tiff');
clear options; 
            options.overwrite = true;
            options.compress = 'lzw';
            saveastiff(uint16(all_nuc), Cell_type_int_name, options); 
            
segmented_nuclei_name = strcat(shortfile,'_segmented_nuclei.tiff');
clear options;
            options.overwrite = true;
            options.compress = 'lzw';
            saveastiff(uint16(nuclei_seg), segmented_nuclei_name, options); 
            
clear sumary_stats_table all_nuc all_nuc_dist_vess nuclei_seg vessel_seg
toc
end
