% Determination of subspace detection threshold using Song et al., 2014
clear
% An example containing: matrix A of event cluster and svd results
load('svd_marmara_G01.mat')

%% Get the subspace detection dimension d from the fractional energy of 
% singular vectors in U (use singular values in S)
% Compute average fractional energy (eq. 23)
clear fc d

for ii = 1:size(A,2)
    fc(ii,1) = trace(S(1:ii,1:ii)'*S(1:ii,1:ii)) / trace(S'*S);
end

figure; plot(1:size(A,2),fc,'-*'); 
xlabel('Dimension of representation d'); ylabel('f_c')
d = 4;

%% Compute the effective dimension of the embedding space Nhat
% Needs to compute the standard deviation of the sample correlation
% coefficient c_hat between noise data and event signal (eq. 27)
addpath('MatSAC') % for rdSac.m
clear c_hat c_std Nhat
NumTries = 200; % Number of noise parts randomly chosen to be compared with events

% Path to continuous data
pathcont = ['/Volumes/Enaffaire/Marmesonet2011/Data/' stan '/' stan opt '/']; 
filesac = dir([pathcont '*.sac']);

% Do correlation between events and random pieces of data
kk = 1;
for ii = 1:size(A,2)
    disp(['Event ' num2str(ii) '/' num2str(size(A,2))])
    for jj = 1:NumTries
        % Get random file
        fid = randi(length(filesac));
        [data,hd] = rdSac([pathcont filesac(fid).name]);
        
        % Get random piece of data
        sid = randi(length(data) - size(A,1));
        data = data(sid:sid+size(A,1)-1);
        %data = randn(601,1);
        % Process data first, as the event data
        data = data - mean(data); % Remove the mean
        data = filtfilt(b,a,data); % Filter as events and detection later
        data = data / max(abs(data)); % Normalize
                
        c_hat(kk,1) = (data(:)'*A(:,ii))/sqrt((A(:,ii)'*A(:,ii))*(data(:)'*data(:)));
        kk = kk+1;
        clear fid data hd sid cctmp iv
    end
end

figure; histogram(c_hat,100); xlabel('Correlation values'); ylabel('#')

c_std = std(c_hat);
Nhat = 1 + (c_std)^(-2); %~797 (no filter), ~99 (with filter)
r = -1:0.01:1;
fr = (1/beta(0.5,(Nhat/2)-1)) * (1 - r.^2).^((Nhat/2)-2); % Eq. 41
hold on; plot(r,60*fr,'--g','LineWidth',2)

clear kk ii jj pathcont filesac r NumTries

%% Chose the correlation detection threshold and get the detection threshold
% Done by setting 1 false alarm per 5400 correlation samples:

disp('Zoom on the bottom of the histogram and press enter')
pause;
disp('choose the limit of the correlation')
[pc,~] = ginput(1);
plot(pc,0,'ko','MarkerSize',8,'MarkerFaceColor','r')
plot(-pc,0,'ko','MarkerSize',8,'MarkerFaceColor','r')
disp(['The correlation detection threshold pc is ' num2str(abs(pc))])
gamma_c = pc^2;

%% Calculate the false alarm rate Pf with eq. 42 (or eq. 33)

Pf = 1 - fcdf((gamma_c / (1 - gamma_c))*((Nhat-1)/1),1,Nhat-1);
disp(['The false alarm rate Pf is ' num2str(Pf)])
