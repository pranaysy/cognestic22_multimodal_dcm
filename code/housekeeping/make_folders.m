% Path to root directory
rootdir = '/imaging/henson/Wakeman/cognestic22_multimodal_dcm';

% Define levels of folder hierarchy 
level1 = {'fits'};
level2 = {'batch_gui', 'batch_script', 'direct_script'};
level3 = {'fmri', 'meg'};
level4 = {'templates'};
level5 = {'DCMs', 'GCMs'};
level6 = {'Full', 'Full_vs_Self', 'Families'}; % For GCMs only

% Loop over each and make folders
for l1=1:numel(level1)
    for l2=1:numel(level2)
        for l3=1:numel(level3)
            for l4=1:numel(level4)
                
                % Make level4 directory
                l4_dir = fullfile(rootdir, level1{l1}, level2{l2}, level3{l3}, level4{l4});
                
                % DCMs
                folderpath = fullfile(l4_dir, level5{1});
                mkdir(folderpath);
                
                % Create empty file
                if exist(fullfile(folderpath, '.empty'), 'file') == 0
                    f = fopen(fullfile(folderpath, '.empty'), 'w' );  
                    fclose(f);
                end
                
                % GCMs
                for l6=1:numel(level6)
                    
                    folderpath = fullfile(l4_dir, level5{2}, level6{l6});
                    mkdir(folderpath);
                    
                    % Create empty file
                    if exist(fullfile(folderpath, '.empty'), 'file') == 0
                        f = fopen(fullfile(folderpath, '.empty'), 'w' );  
                        fclose(f);
                    end
                end
                    
            end
        end
    end
end


% Fancy method because me no like for loops
% Construct grid with combinations of names
% [l1, l2, l3, l4] = ndgrid(1:numel(level1), 1:numel(level2), 1:numel(level3), 1:numel(level4));
% namegrid = vertcat(level1(l1(:))', level2(l2(:)), level3(l3(:)), level4(l4(:))');

