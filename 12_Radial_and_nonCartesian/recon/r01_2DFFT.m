% Reconstruction of 2D Cartesian Pulseq data
% provides an example on how data reordering can be detected from the MR
% sequence with almost no additional prior knowledge
%
% it loads Matlab .mat files with the rawdata in the format 
%     adclen x channels x readouts
% it also seeks an accompanying .seq file with the same name to interpret
%     the data

%% Load the latest file from the specified directory
path='/data/Dropbox/ismrm2021pulseq_liveDemo/dataLive/Vienna_7T_Siemens'; % directory to be scanned for data files

pattern='*.mat';
D=dir([path filesep pattern]);
[~,I]=sort([D(:).datenum]);
data_file_path=[path filesep D(I(end-0)).name]; % use end-1 to reconstruct the second-last data set, etc...
                                                % or replace I(end-0) with I(1) to process the first dataset, I(2) for the second, etc...
% load data
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
end

[adc_len,channels,readouts]=size(data_unsorted);
% the incoming data order is [kx coils acquisitions]
data_coils_last = permute(data_unsorted, [1, 3, 2]);

%% Plot and analyze the trajectory data (ktraj_adc)

[ktraj_adc, ktraj, t_excitation, t_refocusing, t_adc] = seq.calculateKspace();
figure; plot(ktraj(1,:),ktraj(2,:),'b',...
             ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % a 2D plot
axis('equal'); title('2D k-space trajectory');

% try to detect the data ordering
k_extent=max(abs(ktraj_adc),[],2);
k_scale=max(k_extent);
k_threshold=k_scale/5000;

% detect unused dimensions and delete them
if any(k_extent<k_threshold)
    ktraj_adc(k_extent<k_threshold,:)=[]; % delete rows
    k_extent(k_extent<k_threshold)=[];
end

% detect dK, k-space reordering and repetitions (or slices, etc)
kt_sorted=sort(ktraj_adc,2);
dk_all=kt_sorted(:,2:end)-kt_sorted(:,1:(end-1));
dk_all(dk_all<k_threshold)=NaN;
dk_min=min(dk_all,[],2);
dk_max=max(dk_all,[],2);
dk_all(dk_all-dk_min(:,ones(1,size(dk_all,2)))>k_threshold)=NaN;
dk_all_cnt=sum(isfinite(dk_all),2);
dk_all(~isfinite(dk_all))=0;
dk=sum(dk_all,2)./dk_all_cnt;
[~,k0_ind]=min(sum(ktraj_adc.^2,1));
kindex=round((ktraj_adc-ktraj_adc(:,k0_ind*ones(1,size(ktraj_adc,2))))./dk(:,ones(1,size(ktraj_adc,2))));
kindex_min=min(kindex,[],2);
kindex_mat=kindex-kindex_min(:,ones(1,size(ktraj_adc,2)))+1;
kindex_end=max(kindex_mat,[],2);
sampler=zeros(kindex_end');
repeat=zeros(1,size(ktraj_adc,2));
for i=1:size(kindex_mat,2)
    if (size(kindex_mat,1)==3)
        ind=sub2ind(kindex_end,kindex_mat(1,i),kindex_mat(2,i),kindex_mat(3,i));
    else
        ind=sub2ind(kindex_end,kindex_mat(1,i),kindex_mat(2,i)); 
    end
    repeat(i)=sampler(ind);
    sampler(ind)=repeat(i)+1;
end
if (max(repeat(:))>0)
    kindex=[kindex;(repeat+1)];
    kindex_mat=[kindex_mat;(repeat+1)];
    kindex_end=max(kindex_mat,[],2);
end
%figure; plot(kindex(1,:),kindex(2,:),'.-');

%% sort the k-space data into the data matrix
data_coils_last=reshape(data_coils_last, [adc_len*readouts, channels]);

data=zeros([kindex_end' channels]);
if (size(kindex,1)==3)
    for i=1:size(kindex,2)
        data(kindex_mat(1,i),kindex_mat(2,i),kindex_mat(3,i),:)=data_coils_last(i,:);
    end
else
    for i=1:size(kindex,2)
        data(kindex_mat(1,i),kindex_mat(2,i),:)=data_coils_last(i,:);
    end
end

if size(kindex,1)==3
    nImages=size(data,3);
else
    nImages=1;
    data=reshape(data, [size(data,1) size(data,2) 1 size(data,3)]); % we need a dummy images/slices dimension
end

%figure; imab(data);

%% Reconstruct coil images

images = zeros(size(data));
%figure;

for ii = 1:channels
    %images(:,:,:,ii) = fliplr(rot90(fftshift(fft2(fftshift(data(:,:,:,ii))))));
    images(:,:,:,ii) = fftshift(fft2(fftshift(data(end:-1:1,:,:,ii))));
end

%% Image display with optional sum of squares combination
figure;
if channels>1
    sos=abs(sum(images.*conj(images),ndims(images))).^(1/2);
    imab(sos); title('reconstructed image(s), sum-of-squares');
    %sos=sos./max(sos(:));    
    %imwrite(sos, ['img_combined.png']
else
    imab(abs(images));title('reconstructed image(s)');
end
colormap('gray');
saveas(gcf,[basic_file_path '_image_2dfft'],'png');
