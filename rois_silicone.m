
%%This script evaluates the phantom experiments. SNR values for a ROI are
%%calculated. (Actually only a signal value is determined, since the noise
%%variance is the same for all images.)

% load the images of the silicone phantom
load sil_images; %5 psis, 18 flip angles, matrix 100 x 100

%% make ROI
ROI=zeros(100);
radius=9;
circ_x=46;
circ_y=52;
for m=1:128
    for n=1:128
        x=n-circ_x;y=m-circ_y;
        if floor( sqrt(y*y+x*x))<radius+1
            ROI(m,n)=1;
        end
    end
end

%psi=169 was measured on another day, so the ROI is re-positioned
ROI169=zeros(100);
radius=9;
circ_x=47;
circ_y=57;
for m=1:128
    for n=1:128
        x=n-circ_x;y=m-circ_y;
        if floor( sqrt(y*y+x*x))<radius+1
            ROI169(m,n)=1;
        end
    end
end

signal=zeros(5,18);

for p=1:5 % 5 psis
    for k=1:18 % 18 flip angles
        if p<5
            image_cut=ROI.*squeeze(sil_images(p,k,:,:));
        else
            image_cut=ROI169.*squeeze(sil_images(p,k,:,:));
        end
        signal(p,k)=sum(image_cut(:)); %sum over ROI
    end
    signal(p,:)=signal(p,:)./max(signal(p,:)); %normalize
end


%show the ROI
figure;imagesc((1-ROI).*squeeze(sil_images(1,2,:,:)));

%print the normalized SNR values
fprintf('\nnormalized SNR for psi=50° and alpha=[5:5:90]');
signal(1,:)./max(signal(1,:))
fprintf('\nnormalized SNR for psi=115.4° and alpha=[5:5:90]');
signal(2,:)./max(signal(2,:))
fprintf('\nnormalized SNR for psi=117° and alpha=[5:5:90]');
signal(3,:)./max(signal(3,:))
fprintf('\nnormalized SNR for psi=150° and alpha=[5:5:90]');
signal(4,:)./max(signal(3,:))
fprintf('\nnormalized SNR for psi=169° and alpha=[5:5:90]');
signal(5,:)./max(signal(3,:))
fprintf('\n');
