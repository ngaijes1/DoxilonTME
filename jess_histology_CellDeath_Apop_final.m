%% Histology_CellDeath_Apop_final.m Adapted from Benjamin Kingston by Jessica Ngai (Warren Chan lab)
matlab_folder = pwd;

%Select folder contining 3 channel images
Main_folder = uigetdir();
cd(Main_folder);

%% Define cell type and key terms for directories
files = dir('*MOD*');
cell_marker = 3; 
%Options are: 
% 1 = F480
% 3 = Ly6g

save_dir = strcat(Main_folder,'\Results with Apoptosis 200 thresh + Ly6G 800 thresh'); %change cell type's bkg threshold

    if exist(save_dir, 'dir')~=7
        mkdir(save_dir);
    end

%%%%%%%%%%%%%%%%%%%%%    
cd(Main_folder)
%% CHANGE FOR NAME OF SECONDARY ONLY FILES
sec_files = dir('*Secondary Only*'); 

if cell_marker == 1   %F480
    img_file = [dir('*F4_80-Cy5_01*.tiff');dir('*F4_80-Cy5*.tiff')]; 
   
else cell_marker == 3   %Ly6G 
    img_file = [dir('*Ly6G-Cy5*.tiff');dir('*Ly6G-Cy5_01*.tiff')];
    
end

for  i = 1:size(img_file,1)  %loop the size of array of img_file
tic                          
 cd(Main_folder)             
 
     [~,shortfile] = fileparts(img_file(i).name);   %get parts of file name, just name as string array
display(['Analyzing ' shortfile])                   

whole_img = loadtiff(img_file(i).name);             %loads tiff file to temp var

dead = whole_img(:,:,2);                %save dead stain data/channel to temp dead variable
cell_type = whole_img(:,:,3);           %save cell type data/channel to temp cell_type variable

dead_brightspots = stretch(dead)<253;
dead_brightspots = uint16(dead_brightspots);
dead = dead.*dead_brightspots;

dead_thresh_lvl = adaptthresh(dead,0.4);         
dead_thresh = imbinarize(dead,dead_thresh_lvl);  

%% Create the mask or edge of tissue using Caspase3
tissue_boundary_thresh = dead>200;      %accept only positions of caspase3 channel greater than 200 (intensity)
level_erode = strel('sphere',3);        
level_dilate = strel('sphere',3);       
tissue_boundary_erode = imerode(tissue_boundary_thresh,level_erode);    
tissue_boundary_dilate = imdilate(tissue_boundary_erode,level_dilate);  

dt_tissue = bwdist(tissue_boundary_dilate); 
dt_tissue_boundary = dt_tissue<10;          
dt_tissue_boundary_fill = imfill(dt_tissue_boundary,'holes');   
dt_tissue_boundary_fill_rmove_small = bwareaopen(dt_tissue_boundary_fill,10000);  

tissue_boundary = dt_tissue_boundary_fill_rmove_small;  
dead_thresh_crop = dead_thresh.*tissue_boundary;        %element-wise multiplication of matrix to apply tissue boundary onto dapi matrix
%% Dilating the CASPASE for dead cell mask
threshnuc = dead_thresh_crop>0;         
    dt_threshnuc = dt(threshnuc);       
    seeds = maxima(dt_threshnuc,2,0);   
    seeds2 = dilation(seeds,4.5)>0;     
    image_out = waterseed(seeds2,max(dt_threshnuc)-dt_threshnuc,1,0,0); 
    threshnuc(image_out) = false;       
                                        
tissue_boundary_rmove_edge_dt = dt(tissue_boundary);            
tissue_boundary_rmove_edge = tissue_boundary_rmove_edge_dt>70;  
dapi_thresh_crop_noedge = tissue_boundary_rmove_edge.*threshnuc; 

dead_seg = uint8(dapi_thresh_crop_noedge)>0;      
SmallDeadRemove = bwareaopen(dead_seg, 5);      
LargeDead = bwareaopen(SmallDeadRemove, 275,4); 
Dead_cleaned = SmallDeadRemove - LargeDead;
dead_seg = Dead_cleaned;

dead_label = label(dead_seg>0); 
nuc_dilated_label = dip_growregions(dead_label,[],[],1,1,'low_first'); 
dead_label = double(nuc_dilated_label); 

%% Assigning cell type 
all_cells_stats = regionprops(dead_label,cell_type,'MeanIntensity','PixelIdxList'); 
all_cells_stats_struct = all_cells_stats; 

max_cell_in_tissue = max(dead_label(:));   
min_cell_in_tissue = min(dead_label(:))+1; 

new_all_dead =  uint32(dead_label);     

    for b = min_cell_in_tissue:max_cell_in_tissue
            new_all_dead(all_cells_stats_struct(b).PixelIdxList) = all_cells_stats_struct(b).MeanIntensity;            
    end
    
    all_dead = new_all_dead; 
    
%% Calculate background from secondary only controls
cd(Main_folder)
means_table_dead = zeros(2,2); 

for  j = 1:size(sec_files,1)

    secondary_img = loadtiff(sec_files(j).name);    
    secondary_img_dead = secondary_img(:,:,3);      
    crop_secondary_img_dead = secondary_img_dead>800;  %CHANGE (800 for CC3+F480)  (800 for CC3+Ly6g)
    crop_secondary_img_done_dead = secondary_img_dead.*uint16(crop_secondary_img_dead); %applying the mask on actual img
    mean_background_dead = mean(crop_secondary_img_done_dead(:));   
    std_background_dead = std(double(crop_secondary_img_done_dead(:))); 
 
    means_table_dead(j,1) = mean_background_dead;  
    means_table_dead(j,2) = std_background_dead;   
       
end

mean_background_dead = mean(means_table_dead(:,1)); 
std_background_dead = mean(means_table_dead(:,2));  


%% Dead Cells_postive cut-offs
background_sing_thresh_3std = mean_background_dead+(3*std_background_dead); 
background_sing_thresh_2std = mean_background_dead+(2*std_background_dead);

%WITHIN 2 STDEV
test_dead_std2 = all_dead>background_sing_thresh_2std; %compare if nuc intensity is > than 2 std threshold bkg; logical array
label_test_dead_2std = bwlabel(test_dead_std2);           
num_mac_dead_2std = max(label_test_dead_2std(:));         
num_total = max(dead_label(:));                           

percent_dead_2std = num_mac_dead_2std/num_total*100;  

%WITHIN 3 STDEV
test_dead_std3 = all_dead>background_sing_thresh_3std; %compare if nuc intensity is > than 3 std threshold bkg; logical array 
    label_test_dead_3std = bwlabel(test_dead_std3);       
num_mac_dead_3std = max(label_test_dead_3std(:));          
num_total = max(dead_label(:));                   

percent_dead_3std = num_mac_dead_3std/num_total*100;  

%% Combine Results Tables 

summary_stats = [background_sing_thresh_2std percent_dead_2std background_sing_thresh_3std percent_dead_3std num_mac_dead_2std num_mac_dead_3std num_total];
sumary_stats_table = array2table(summary_stats);
sumary_stats_table.Properties.VariableNames = {'TWOSTD_Signal_Cutoff' 'TWOSTD_Percent_dead_positive' 'THREESTD_Signal_Cut_off' 'THREESTD_Percent_dead_positive' 'NUM_CELLS_Dead_twostd' 'NUM_CELLS_Dead_threestd' 'TOTAL_Dead_NUMBER'};

%% Write Results
cd(save_dir)

summary_stats_table_name = strcat(shortfile,'_Results_Summary_stats.csv');
writetable(sumary_stats_table,summary_stats_table_name);

Cell_type_int_name = strcat(shortfile,'_cell_type_intensity.tiff');
clear options; 
            options.overwrite = true;
            options.compress = 'lzw';
            saveastiff(uint16(all_dead), Cell_type_int_name, options); 

segmented_nuclei_name = strcat(shortfile,'_segmented_dead.tiff');
clear options;
            options.overwrite = true;
            options.compress = 'lzw';
            saveastiff(uint16(dead_seg), segmented_nuclei_name, options); 
                   
clear sumary_stats_table all_nuc all_nuc_dist_vess nuclei_seg vessel_seg
toc
end
