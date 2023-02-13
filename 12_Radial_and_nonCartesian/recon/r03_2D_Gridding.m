% very basic and crude non-Cartesian recon using griddata()
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
elseif ndims(data_unsorted)==2 && size(data_unsorted,2)>1
    data_unsorted=reshape(data_unsorted,[size(data_unsorted,1),1,size(data_unsorted,2)]);
end

[adc_len,channels,readouts]=size(data_unsorted);
rawdata = double(permute(data_unsorted, [1,3,2])); % channels last
rawdata = reshape(rawdata, [adc_len*readouts,channels]);

%% reconstruct the trajectory 

traj_recon_delay=[0 0 0]*1e-6; % adjust this parameter to potentially improve resolution & geometric accuracy. 
                       % It can be calibrated by inverting the spiral revolution dimension and making 
                       % two images match. for our Prisma and a particular trajectory we found 1.75e-6
                       % it is also possisible to provide a vector of 3 delays (varying per axis)

[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP('trajectory_delay',traj_recon_delay); 

% detect slice dimension
max_abs_ktraj_adc=max(abs(ktraj_adc'));
[~, slcDim]=min(max_abs_ktraj_adc);
encDim=find([1 2 3]~=slcDim);

% figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
% 
figure; plot(ktraj(encDim(1),:),ktraj(encDim(2),:),'b',...
             ktraj_adc(encDim(1),:),ktraj_adc(encDim(2),:),'r.'); % a 2D plot
axis('equal'); title('2D k-space trajectory');

%% Define FOV and resolution and simple off-resonance frequency correction 

%fov=30e-3; Nx=128; Ny=Nx; % OCRA
fov=256e-3; Nx=256; Ny=Nx; % whole-body scanners
deltak=1/fov;
os=2; % oversampling factor (we oversample both in image and k-space)
offresonance=0; % global off-resonance in Hz

%% rudimentary off-resonance correction
nex=length(t_excitation);
t_adc_ex=t_adc;
if nex>1
    for e=2:nex
        i1=find(t_adc>t_excitation(e),1);
        if e<nex
            i2=max(find(t_adc<t_excitation(e+1)));
        else
            i2=length(t_adc);
        end
        t_adc_ex(i1:i2)=t_adc_ex(i1:i2)-t_excitation(e);
    end
end

for c=1:channels
    rawdata(:,c) = rawdata(:,c) .* exp(-1i*2*pi*t_adc_ex'*offresonance);
end

%% here we expect Nx, Ny, deltak to be set already
% and rawdata ktraj_adc loaded (and having the same dimensions)

kxm=round(os*os*Nx/2);
kym=round(os*os*Ny/2);

[kyy,kxx] = meshgrid(-kxm:(kxm-1), -kym:(kym-1));
kyy=-kyy*deltak/os; % we swap the order ind invert one sign to account for Matlab's strange column/line convention
kxx=kxx*deltak/os;

kgd=zeros([size(kxx) channels]);
for c=1:channels
    kgd(:,:,c)=griddata(ktraj_adc(encDim(1),:),ktraj_adc(encDim(2),:),rawdata(:,c),kxx,kyy,'cubic'); 
end
kgd(isnan(kgd))=0;

figure;imagesc(log(abs(kgd(:,:,1))));axis('square');title('k-space data after gridding');

igd=ifftshift(ifft2(ifftshift(kgd())));

Nxo=round(Nx*os);
Nyo=round(Ny*os);
Nxs=round((size(igd,1)-Nxo)/2);
Nys=round((size(igd,2)-Nyo)/2);
igdc = igd((Nxs+1):(Nxs+Nxo),(Nys+1):(Nys+Nyo),:);
if slcDim==1
    igdc=rot90(igdc,-1); % this makes sagittal images look more natural
end
figure;imab(abs(igdc));colormap('gray');
%axis('equal');

%% Sum of squares combination
if channels>1
    sos=abs(sum(igdc.^2,ndims(igdc)).^(1/2));
    sos=sos./max(sos(:));
    figure;imab(sos);colormap('gray');
    %imwrite(sos, ['img_combined.png'])
end

saveas(gcf,[basic_file_path '_image_2dgrd'],'png');
