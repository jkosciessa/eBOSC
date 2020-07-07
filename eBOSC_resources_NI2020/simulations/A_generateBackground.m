% This script generates the scale-free background that is used for the
% simulations of standard BOSC and eBOSC.

% Note that the function 'f_alpha_gaussian' from the CNOISE toolbox 
% has to be added to the path. The toolbox is available at
% https://people.sc.fsu.edu/~jburkardt/m_src/cnoise/cnoise.html

% For reproducibility, we use a fixed seed.

pn.backgroundData = ""; % set output path

% initialize background matrix
bckgrnd      = zeros(1000,5000);
bckgrnd_filt = zeros(1000,5000);

% seed
randn('seed',20160118);

% loop 1000 repetitions
for k = 1:1000

    display(num2str(k))

    % generate 1/f background
    bckgrnd(k,:) = f_alpha_gaussian(5000,1,1);

    % bandpass filter signal (consecutive low + high-pass filter)
    
    [B,A]  = butter(4,70/(250/2),'low'); 
    bckgrnd_filt(k,:) = filtfilt(B,A,bckgrnd(k,:)); clear A B
    
    [B,A]  = butter(4,.5/(250/2),'high'); 
    bckgrnd_filt(k,:) = filtfilt(B,A,bckgrnd_filt(k,:)); clear A B

end; clear k

% save background
save([pn.out 'background'],'bckgrnd','bckgrnd_filt')
