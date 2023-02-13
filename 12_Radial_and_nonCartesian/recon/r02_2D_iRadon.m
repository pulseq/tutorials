% very basic inverse Radon radial reconstruction
%
% it loads Matlab .mat files with the rawdata in the format 
%     adclen x channels x readouts
% it also seeks an accompanying .seq file with the same name to interpret
%     the data

%% Load the latest file from the specified directory
%path='/data/Dropbox/ismrm2021pulseq_liveDemo/dataLive/Vienna_7T_Siemens'; % directory to be scanned for data files

path='~/20211025-AMR/data';

pattern='*.seq';
D=dir([path filesep pattern]);
[~,I]=sort([D(:).datenum]);
data_file_path=[path filesep D(I(18)).name]; % use end-1 to reconstruct the second-last data set, etc...
                                             % or replace I(end-0) with I(1) to process the first dataset, I(2) for the second, etc...
  
% basic path and filename without the extension
[p,n,e] = fileparts(data_file_path);
basic_file_path=fullfile(p,n);

% load data
data_file_path=[basic_file_path '.mat'];
fprintf(['loading `' data_file_path '´ ...\n']);
data_unsorted = load(data_file_path);
if isstruct(data_unsorted)
    fn=fieldnames(data_unsorted);
    assert(length(fn)==1); % we only expect a single variable
    data_unsorted=data_unsorted.(fn{1});
end
% keep basic filename without the extension
[p,n,e] = fileparts(data_file_path);
basic_file_path=fullfile(p,n);

%% Load sequence from file 
seq = mr.Sequence();              % Create a new sequence object
seq_file_path = [basic_file_path '.seq'];
seq.read(seq_file_path,'detectRFuse'); % detectRFuse is an important option for SE sequences
fprintf(['loaded sequence `' seq.getDefinition('Name') '´\n']);

%% raw data preparation

if ndims(data_unsorted)<3 && size(data_unsorted,1)==1
    % OCRA data may need some fixing
    [~, ~, eventCount]=seq.duration();
    readouts=eventCount(6);
    data_unsorted=reshape(data_unsorted,[length(data_unsorted)/readouts,1,readouts]);
elseif ndims(data_unsorted)==2 && size(data_unsorted,2)>1
    data_unsorted=reshape(data_unsorted,[size(data_unsorted,1),1,size(data_unsorted,2)]);
end

[adc_len,channels,readouts]=size(data_unsorted);
rawdata = permute(data_unsorted, [1,3,2]); % channels last

%% Analyze the nominal trajectory

[ktraj_adc_nom,t_adc] = seq.calculateKspacePP('trajectory_delay',0); 

% detect slice dimension
max_abs_ktraj_adc=max(abs(ktraj_adc_nom'));
[~, slcDim]=min(max_abs_ktraj_adc);
encDim=find([1 2 3]~=slcDim);

ktraj_adc_nom = reshape(ktraj_adc_nom, [3, adc_len, size(ktraj_adc_nom,2)/adc_len]);

prg_angle=unwrap(squeeze(atan2(ktraj_adc_nom(encDim(2),2,:)-ktraj_adc_nom(encDim(2),1,:),...
                               ktraj_adc_nom(encDim(1),2,:)-ktraj_adc_nom(encDim(1),1,:))));
nproj=length(prg_angle);

%% from k-space to projections (1D FFTs)

data_fft1=ifftshift(ifft(ifftshift(rawdata,1)),1);

figure; imab(abs(squeeze(data_fft1))); title('sinogramm view')

%% crop the data to remove oversampling and adapt to the iradon
target_matrix_size=256;
shift=0; % 12 was found experimentally for the first Benjamin's data set
cropLeft=(adc_len-target_matrix_size)/2+shift;
cropRight=(adc_len-target_matrix_size)/2-shift;
data_fft1c=data_fft1(2+cropLeft:end-cropRight,:,:);

% visualize the matching of positive and negative directions 
p1=1;
[~,p2]=min(abs(mod(prg_angle-prg_angle(p1)+2*pi,2*pi)-pi)); % look for the projection closest to the opposite to p1 
figure;plot(abs(data_fft1c(1:end,p1,1)));hold on;plot(abs(data_fft1c(end:-1:1,p2,1))); title('comparing opposite projections');

%% the actuall iRadon transform
theta=270-prg_angle/pi*180;
for c=1:channels
    % the classical (absolute value) transform
    irad_a=iradon(abs(data_fft1c(:,:,c)),theta,'linear','Hann');
    %irad_a=iradon(abs(data_fft1c(:,:,c)),theta);
    %irad_a=iradon(abs(data_fft1c(:,:,c)),theta,'linear','Shepp-Logan');
    if (c==1)
        irad_abs=zeros([size(irad_a) channels]);
        irad_cmpx=zeros([size(irad_a) channels]);
    end
    irad_abs(:,:,c)=irad_a;
    % MR-specific complex-valued transform
    irad_r=iradon(real(data_fft1c(:,:,c)),theta,'linear','Hann');
    irad_i=iradon(imag(data_fft1c(:,:,c)),theta,'linear','Hann');
    irad_cmpx(:,:,c)=irad_r + 1i*irad_i;
end

figure;imab(abs(irad_abs));colormap('gray'); title('abs iRadon recon')
saveas(gcf,[basic_file_path '_image_iradon'],'png');

figure;imab(abs(irad_cmpx));colormap('gray'); title('complex iRadon recon')
%axis('equal');

%% Sum of squares combination
if channels>1
    sos=abs(sum(irad_cmpx.^2,ndims(irad_cmpx)).^(1/2));
    sos=sos./max(sos(:));
    figure;imab(sos);colormap('gray'); title('SOS complex iRadon recon')
    %imwrite(sos, ['img_combined.png'])
end
