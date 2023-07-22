% very basic and crude non-Cartesian recon using griddata()
%
% it loads Matlab .mat files with the rawdata in the format 
%     adclen x channels x readouts
% it also seeks an accompanzing .seq file with the same name to interpret
%     the data

%% Load the latest file from the specified directory
%path='../IceNIH_RawSend/'; % directory to be scanned for data files
%path='/data/Dropbox/ismrm2021pulseq_liveDemo/dataLive/Vienna_7T_Siemens'; % directory to be scanned for data files
%path='/data/Dropbox/ismrm2021pulseq_liveDemo/dataPrerecorded/Vienna_7T_Siemens'
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

[ktraj_adc_nom,t_adc] = seq.calculateKspacePP('trajectory_delay',0e-6); 

% detect slice dimension
max_abs_ktraj_adc=max(abs(ktraj_adc_nom'));
[~, slcDim]=min(max_abs_ktraj_adc);
encDim=find([1 2 3]~=slcDim);

ktraj_adc_nom = reshape(ktraj_adc_nom, [3, adc_len, size(ktraj_adc_nom,2)/adc_len]);

prg_angle=squeeze(atan2(ktraj_adc_nom(encDim(2),2,:)-ktraj_adc_nom(encDim(2),1,:),ktraj_adc_nom(encDim(2),2,:)-ktraj_adc_nom(encDim(2),1,:)));
nproj=length(prg_angle);

i_pure2=find(abs(ktraj_adc_nom(encDim(2),1,:))<eps);
i_pure1=find(abs(ktraj_adc_nom(encDim(1),1,:))<eps);
assert(length(i_pure2)==length(i_pure1)); % the code below assumes this

delta180=nproj/length(i_pure2);

%%

data_fft1=ifftshift(ifft(ifftshift(rawdata,1)),1);

ip=1:nproj;
cmplx_diff=data_fft1(2:end,ip,:).*conj(data_fft1(end:-1:2,mod(ip-1+delta180,nproj)+1,:));
cmplx_diff_no_channels=sum(cmplx_diff,3);

figure; imagesc(angle(cmplx_diff_no_channels));

%% pick the pure X and pure Y differences, plot them and estimate the slope
if length(i_pure2)<2
    [~,p2]=min(abs(mod(prg_angle-prg_angle(1)+2*pi,2*pi)-pi)); % look for the projection closest to the opposite to the 1st 
    i_pure1=[1 p2];
    i_pure2=i_pure1; % this is a hack!
end

thresh=0.05;
cmplx_diff_pure_axes=cmplx_diff_no_channels(:,[i_pure2 i_pure1]);
mpa=max(abs(cmplx_diff_pure_axes));
for i=1:size(cmplx_diff_pure_axes,2)
    cmplx_diff_pure_axes(abs(cmplx_diff_pure_axes(:,i))<thresh*mpa(i),i)=0;
end

% plot
figure; plot(angle(cmplx_diff_pure_axes(:,[1 3])));

% just get the mean slope of the phase
msop=angle(sum(cmplx_diff_pure_axes(2:end,:).*conj(cmplx_diff_pure_axes(1:end-1,:))));
delays_imageSpace=-msop/2/pi/2*(t_adc(2)-t_adc(1))*adc_len; 

if all(i_pure1==i_pure2)    
    fprintf('this is not a full sweep scan so we can only roughly estimate the delay as: %g us\n', round(delays_imageSpace(1)*1e9)/1e3);
    fprintf('you can try to re-run recon04_2D_Gridding setting traj_recon_delay to the above delay value\n');
    return
end

%% a more accurate function
delays=[0 0 0];
delays(encDim)=calc_delays(rawdata, ktraj_adc_nom(encDim,:),i_pure2,i_pure1,t_adc(2)-t_adc(1));

% improve iteratively & verify
n_it=4;
d=[0 0 0];
for n=1:n_it
    d=sum(delays,1);
    ktraj_adc = seq.calculateKspacePP('trajectory_delay',d);
    d(encDim)=calc_delays(rawdata, ktraj_adc(encDim,:),i_pure2,i_pure1,t_adc(2)-t_adc(1));
    delays = [delays; d];
end
d=sum(delays,1);

fprintf('found delays: [%g %g %g] us\n', round(d*1e9)/1e3);
fprintf('now re-run recon04_2D_Gridding setting traj_recon_delay to the above delay values\n');

% END 

function delays=calc_delays(rawdata, ktraj_adc,i_pureX,i_pureY,dt)
    adc_len=size(rawdata,1);
    n_chan=size(rawdata,3);
    n_sel=length(i_pureX)+length(i_pureY);
    data_selected=rawdata(:,[i_pureX i_pureY],:);
    data_sel_res=zeros(adc_len,n_sel,n_chan);
    dk=zeros(1,n_sel);abs(diff(ktraj_adc(1,(i_pureX(1)-1)*adc_len+(1:2))));
    for i=1:n_sel
        if i<=length(i_pureX)
            kr_selected=ktraj_adc(1,(i_pureX(i)-1)*adc_len+(1:adc_len));
        else
            kr_selected=ktraj_adc(2,(i_pureY(i-length(i_pureX))-1)*adc_len+(1:adc_len));
        end
        [~,ikmin]=min(abs(kr_selected));
        kmin=kr_selected(ikmin);
        dk(i)=diff(kr_selected(ikmin:ikmin+1));
        % in casa of large delays we can get out of the gradient and need to
        % fix this part not to confuse the interpolation function 
        d=(kr_selected(2:end)-kr_selected(1:end-1))/dk(i);
        d(d<=1e-10)=1e-10; % replace zeros or negative steps with very small positive numbers
        kr_selected=dk(i)*cumsum([0 d]);
        kr_selected=kr_selected-kr_selected(ikmin)+kmin;
        % done fixing...
        for c=1:n_chan
            data_sel_res(:,i,c)=interp1(kr_selected,data_selected(:,i,c),(-adc_len/2:(adc_len/2-1))*abs(dk(i)),'pchip',0);
        end
    end

    data_sel_fft1=ifftshift(ifft(ifftshift(data_sel_res,1)),1);

    thresh=0.05;
    cmplx_diff_sel=data_sel_fft1(:,2:2:end,:).*conj(data_sel_fft1(:,1:2:end,:));
    cmplx_diff_sel_no_channels=squeeze(sum(cmplx_diff_sel,3));
    mpa=max(abs(cmplx_diff_sel_no_channels));
    for i=1:size(cmplx_diff_sel_no_channels,2)
        cmplx_diff_sel_no_channels(abs(cmplx_diff_sel_no_channels(:,i))<thresh*mpa(i),i)=0;
    end
    
    %figure; plot(angle(cmplx_diff_sel_no_channels));

    % just get the mean slope of the phase
    msop1=angle(sum(cmplx_diff_sel_no_channels(2:end,:).*conj(cmplx_diff_sel_no_channels(1:end-1,:))));
    delays=msop1/2/pi/2*dt*adc_len.*sign(dk(1:2:end)); % we need thins "sign" to account for the direction (positive or negative) of the first projection in each pair
end
