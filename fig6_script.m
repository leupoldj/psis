%Simulation of Fig 6.

TR=50; %change this to TR=50 for bottom row in Fig. 6

%now choose H2O, GM or silicone oil:
substance = 1;
% H2O:              substance = 1
% grey matter (GM): substance = 2
% silicone oil:     substance = 3

%H2O + CuSO4
if substance == 1
T1=540;
T2=340;
D=1.93e-3; %in mm^2/s
end

%GM
if substance == 2
T1=1500;
T2=100;
D=0.8e-3; %in mm^2/s
end

%Silicone Oil
if substance == 3
T1=1290;
T2=399;
D=0.0055e-3; %in mm^2/s
end

psi_vec=[50 115.4 117 150 169];

% for H2O or GM
if (substance == 1 || substance == 2)
    res_vec=[0:0.02:5]; %divisor for voxel size, e.g. res_vec = 5 corresponds to voxel size of 300um/5 = 60um.
end

% for silicone oil: 
if (substance == 3)
    res_vec=[0:0.1:70];
end

%flip angles
flip_angles=[5:5:90];

%number of pulses
np=1000;

%allocate memory 
signal_plus=zeros(length(psi_vec),length(flip_angles),length(res_vec)); %signal S+
ernst_ampl=zeros(1,length(flip_angles)); %ernst amplitudes

%signal computation
for p=1:length(psi_vec);
    for a=1:length(flip_angles)
        tic
        [p a]
        parfor r=1:length(res_vec)
            [temp temp2]=epg_rfsp(flip_angles(a),np,T1,T2,TR,D,psi_vec(p),res_vec(r));
            signal_plus(p,a,r) = abs(temp(np));
        end %eof flip_angles
        toc
    end %eof flip_angles
end %eof psi_vec

%% compute the epsilon-values

for a=1:length(flip_angles)    
    ernst_ampl(a)=sind(flip_angles(a))*(1-exp(-TR/T1))/(1-exp(-TR/T1)*cosd(flip_angles(a)));
end

flip_indices=[1:length(flip_angles)];
epsilon=zeros(length(psi_vec),length(res_vec));

for p=1:length(psi_vec);
    for r=1:length(res_vec)     
        for a=1:length(flip_indices)
            sig=signal_plus(p,flip_indices(a),r);
            diff=sig-ernst_ampl(flip_indices(a));
            epsilon(p,r)=epsilon(p,r)+abs(diff)*(1/length(flip_indices)*1./ernst_ampl(flip_indices(a)));
        end %eof flip_angles
    end %eof flip_angles
end %eof psi_vec

figure;
for p=1:length(psi_vec)
    hold on;
    plot(res_vec,100*epsilon(p,:),'LineWidth',2);
end
legend('\psi = 50°','\psi = 115.4°','\psi = 117°','\psi = 150°','\psi = 169°','Location','northeast');
axis([0 5 0 25]);
xticks([0.3 1 2 3 4 5])
xticklabels({'1000','300','150','100','75','60'})
ylabel('epsilon [%]','FontSize',14);
xlabel('voxel size [\mum]','FontSize',14);
set(gca,'FontSize',14);
if substance == 1 
    figname = ['H2O+CuSO4, TR = ' num2str(TR) ' ms'];
end
if substance == 2 
    figname = ['GM, TR= ' num2str(TR) ' ms'];
end
if substance == 3 
    figname = ['silicone oil, TR= ' num2str(TR) ' ms'];
end
title(figname)

% only for silicone: plot also larger range
if substance == 3 
    figure;
    for p=1:length(psi_vec)
        hold on;
        plot(res_vec,100*epsilon(p,:),'LineWidth',2);
    end
    legend('\psi = 50°','\psi = 115.4°','\psi = 117°','\psi = 150°','\psi = 169°','Location','northeast');
    axis([0 70 0 20]);
    xticks([1 10 20 30 40 50 60 70])
    xticklabels({'300','30','15','10','7.5','6','5','4.3'})
    ylabel('epsilon [%]','FontSize',14);
    xlabel('voxel size [\mum]','FontSize',14);
    set(gca,'FontSize',14);
    figname = ['silicone oil, TR= ' num2str(TR) ' ms'];
    title(figname)
end