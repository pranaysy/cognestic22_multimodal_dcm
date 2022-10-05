% Path to root directory
rootdir = '/imaging/henson/Wakeman/cognestic22_multimodal_dcm';

% Fits folder
fitsdir = fullfile(rootdir, 'fits');

% Populate all .mat files
filelist = dir(fullfile(fitsdir, '**', '*.mat'));

% Get default recycling state
state = recycle;
recycle('on') % Permanently delete!

% Delete everything (commented for safety by default!)
for f=1:length(filelist)
    filename = fullfile(filelist(f).folder, filelist(f).name);
    %delete(filename) % Uncomment to actually delete *.mat files in the 'fits' folder
end

% Reset recycling behaviour back to default
recycle(state)
