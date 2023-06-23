%parameters for M100 silicone oil
T1=1290;
T2=399;
D=0.0055e-3; %in mm^2/s

TR=20;
no_pulses=1000;
flip_angles=[5:5:90];

%allocate memory
sig117=zeros(length(flip_angles),no_pulses);
sig115=zeros(length(flip_angles),no_pulses);
sig50=zeros(length(flip_angles),no_pulses);
sig150=zeros(length(flip_angles),no_pulses);
sig169=zeros(length(flip_angles),no_pulses);
ernst=zeros(1,length(flip_angles));

%compute the signals
parfor a=1:length(flip_angles)
    [S_plus S_minus]=epg_rfsp(flip_angles(a),no_pulses,T1,T2,TR,D,50,1);
    sig50(a,:)=abs(S_plus);
    [S_plus S_minus]=epg_rfsp(flip_angles(a),no_pulses,T1,T2,TR,D,115.4,1);
    sig115(a,:)=abs(S_plus);
    [S_plus S_minus]=epg_rfsp(flip_angles(a),no_pulses,T1,T2,TR,D,117,1);
    sig117(a,:)=abs(S_plus);
    [S_plus S_minus]=epg_rfsp(flip_angles(a),no_pulses,T1,T2,TR,D,150,1);
    sig150(a,:)=abs(S_plus);
    [S_plus S_minus]=epg_rfsp(flip_angles(a),no_pulses,T1,T2,TR,D,169,1);
    sig169(a,:)=abs(S_plus);
    ernst(a)=sind(flip_angles(a))*(1-exp(-TR/T1))/(1-exp(-TR/T1)*cosd(flip_angles(a)));
end

%snr data per flip angle for silicone oil
snr50=[7.63 9.55 8.66 7.01 5.79 4.91 4.02 3.40 2.90 2.58 2.12 1.98 1.86 1.74 1.62 1.46 1.46 1.50];
snr115=[7.70 9.54 8.56 7.15 6.05 5.00 4.25 3.89 3.57 3.23 2.95 2.68 2.39 2.22 1.94 1.71 1.44 1.26];
snr117=[7.74 9.77 8.70 7.12 6.32 5.70 5.44 5.13 4.76 4.51 4.24 3.85 3.46 2.97 2.72 2.45 2.03 1.77];
snr150=[7.63 9.48 8.68 7.55 6.49 5.57 4.47 3.97 3.33 2.69 2.27 1.96 1.65 1.55 1.31 1.22 1.09 0.99];
snr169=[8.18 10.2 9.14 7.66 6.30 5.19 4.44 3.83 3.33 3.08 2.86 2.70 2.52 2.34 2.16 2.02 1.80 1.69];

%plot simulation + experimental data
figure;plot(flip_angles,sig50(:,no_pulses)/squeeze(max(sig50(:,no_pulses))),'b','LineWidth',1)
hold on;plot(flip_angles,sig115(:,no_pulses)/squeeze(max(sig115(:,no_pulses))),'r','LineWidth',1)
hold on;plot(flip_angles,sig117(:,no_pulses)/squeeze(max(sig117(:,no_pulses))),'g','LineWidth',1)
hold on;plot(flip_angles,sig150(:,no_pulses)/squeeze(max(sig150(:,no_pulses))),'cy','LineWidth',1)
hold on;plot(flip_angles,sig169(:,no_pulses)/squeeze(max(sig169(:,no_pulses))),'m','LineWidth',1)

hold on;plot(flip_angles,ernst/max(ernst),'k','LineWidth',2)
hold on;plot(flip_angles,snr50/max(snr50),'b+','MarkerSize',10);
hold on;plot(flip_angles,snr115/max(snr115),'ro','MarkerSize',10);
hold on;plot(flip_angles,snr117/max(snr117),'g*','MarkerSize',10);
hold on;plot(flip_angles,snr150/max(snr150),'cyd','MarkerSize',10);
hold on;plot(flip_angles,snr169/max(snr169),'mo','MarkerSize',6);

title('M100 silicone oil');
ylabel('Signal [a.u.]','FontSize',14);
xlabel('flip angle [deg]','FontSize',14);
set(gca,'XTick',[10:10:90],'TickDir','out')
set(gca,'FontSize',14);
legend('\psi = 50°','\psi = 115.4°','\psi = 117°','\psi = 150°','\psi = 169°','Ernst Ampl.','Location','northeast');
lgd = legend;
lgd.FontSize = 14;
%%
% compute the epsilon-values for silicone oil (see manuscript)
flip_ind=[1:18];
psi_vec=[50 115.4 117 150 169];
epsilon=zeros(1,length(psi_vec));
a=[sig50(:,no_pulses) sig115(:,no_pulses) sig117(:,no_pulses) sig150(:,no_pulses) sig169(:,no_pulses)];
b=repmat(ernst',1,length(psi_vec));
ea=acosd(exp(-TR/T1)); %Ernst angle - only for convenience, actually it is not used further
sig_ea=sind(ea)*(1-exp(-TR/T1))/(1-exp(-TR/T1)*cosd(ea)); %signal at the Ernst angle  - only for convenience, actually it is not used further
diffs=abs(a(flip_ind,:)-b(flip_ind,:))./ernst';
epsilon=100*sum(diffs,1)/18;

fprintf('\nP-values for silicone oil in percent:');
fprintf('\n50°: %f  115.4°: %f  117°: %f  150°: %f  169°: %f\n',epsilon(find(psi_vec==50)),epsilon(find(10*psi_vec==1154)),epsilon(find(psi_vec==117)),epsilon(find(psi_vec==150)),epsilon(find(psi_vec==169)));


%% now simulate H20+CuSO4. Note the variables from above (i.e. those for the silicone) are overwritten.

% parameters for H20+CuSO4
T1=540;
T2=340;
D=1.93e-3; %in mm^2/s

%allocate memory
sig117=zeros(length(flip_angles),no_pulses);
sig115=zeros(length(flip_angles),no_pulses);
sig50=zeros(length(flip_angles),no_pulses);
sig150=zeros(length(flip_angles),no_pulses);
sig169=zeros(length(flip_angles),no_pulses);
ernst=zeros(1,length(flip_angles));

% allocate memory for the D=0 simulation.
sig117D0=zeros(length(flip_angles),no_pulses);
sig115D0=zeros(length(flip_angles),no_pulses);
sig50D0=zeros(length(flip_angles),no_pulses);
sig150D0=zeros(length(flip_angles),no_pulses);
sig169D0=zeros(length(flip_angles),no_pulses);

parfor a=1:length(flip_angles)
    [S_plus S_minus]=epg_rfsp(flip_angles(a),no_pulses,T1,T2,TR,D,50,1);
    sig50(a,:)=abs(S_plus);
    [S_plus S_minus]=epg_rfsp(flip_angles(a),no_pulses,T1,T2,TR,D,115.4,1);
    sig115(a,:)=abs(S_plus);
    [S_plus S_minus]=epg_rfsp(flip_angles(a),no_pulses,T1,T2,TR,D,117,1);
    sig117(a,:)=abs(S_plus);
    [S_plus S_minus]=epg_rfsp(flip_angles(a),no_pulses,T1,T2,TR,D,150,1);
    sig150(a,:)=abs(S_plus);
    [S_plus S_minus]=epg_rfsp(flip_angles(a),no_pulses,T1,T2,TR,D,169,1);
    sig169(a,:)=abs(S_plus);
    
    [S_plus S_minus]=epg_rfsp(flip_angles(a),no_pulses,T1,T2,TR,0,50,1);
    sig50D0(a,:)=abs(S_plus);
    [S_plus S_minus]=epg_rfsp(flip_angles(a),no_pulses,T1,T2,TR,0,115.4,1);
    sig115D0(a,:)=abs(S_plus);
    [S_plus S_minus]=epg_rfsp(flip_angles(a),no_pulses,T1,T2,TR,0,117,1);
    sig117D0(a,:)=abs(S_plus);
    [S_plus S_minus]=epg_rfsp(flip_angles(a),no_pulses,T1,T2,TR,0,150,1);
    sig150D0(a,:)=abs(S_plus);
    [S_plus S_minus]=epg_rfsp(flip_angles(a),no_pulses,T1,T2,TR,0,169,1);
    sig169D0(a,:)=abs(S_plus);
    
    ernst(a)=sind(flip_angles(a))*(1-exp(-TR/T1))/(1-exp(-TR/T1)*cosd(flip_angles(a)));
end

%snr data for H20 + CuSO4
snr50=[5.14 7.92 8.58 8.28 7.48 6.70 5.98 5.36 4.75 4.27 3.85 3.53 3.26 3.03 2.78 2.61 2.49 2.32];
snr115=[5.04 7.90 8.60 8.26 7.53 6.82 6.06 5.45 4.80 4.28 3.85 3.46 3.08 2.77 2.49 2.19 1.99 1.81];
snr117=[5.05 7.89 8.65 8.33 7.66 6.86 6.13 5.45 4.83 4.28 3.86 3.40 3.12 2.77 2.47 2.22 2.02 1.79];
snr150=[5.08 7.85 8.60 8.25 7.52 6.76 6.09 5.39 4.82 4.25 3.87 3.45 3.10 2.80 2.52 2.27 2.04 1.83];
snr169=[4.99 7.67 8.41 8.13 7.51 6.78 6.05 5.49 4.95 4.45 4.08 3.69 3.37 3.07 2.77 2.52 2.31 2.08];

%plot simulation + experimental data
figure;plot(flip_angles,sig50(:,no_pulses)/squeeze(max(sig50(:,no_pulses))),'b','LineWidth',1)
hold on;plot(flip_angles,sig115(:,no_pulses)/squeeze(max(sig115(:,no_pulses))),'r','LineWidth',1)
hold on;plot(flip_angles,sig117(:,no_pulses)/squeeze(max(sig117(:,no_pulses))),'g','LineWidth',1)
hold on;plot(flip_angles,sig150(:,no_pulses)/squeeze(max(sig150(:,no_pulses))),'cy','LineWidth',1)
hold on;plot(flip_angles,sig169(:,no_pulses)/squeeze(max(sig169(:,no_pulses))),'m','LineWidth',1)
hold on;plot(flip_angles,ernst/max(ernst),'k','LineWidth',2)

hold on;plot(flip_angles,sig50D0(:,no_pulses)/squeeze(max(sig50D0(:,no_pulses))),'--b','LineWidth',1)
hold on;plot(flip_angles,sig115D0(:,no_pulses)/squeeze(max(sig115D0(:,no_pulses))),'--r','LineWidth',1)
hold on;plot(flip_angles,sig117D0(:,no_pulses)/squeeze(max(sig117D0(:,no_pulses))),'--g','LineWidth',1)
hold on;plot(flip_angles,sig150D0(:,no_pulses)/squeeze(max(sig150D0(:,no_pulses))),'--cy','LineWidth',1)
hold on;plot(flip_angles,sig169D0(:,no_pulses)/squeeze(max(sig169D0(:,no_pulses))),'--m','LineWidth',1)

hold on;plot(flip_angles,snr50/max(snr50),'b+','MarkerSize',10);
hold on;plot(flip_angles,snr115/max(snr115),'ro','MarkerSize',10);
hold on;plot(flip_angles,snr117/max(snr117),'g*','MarkerSize',10);
hold on;plot(flip_angles,snr150/max(snr150),'cyd','MarkerSize',10);
hold on;plot(flip_angles,snr169/max(snr169),'mo','MarkerSize',6);

title('H20+CuSO4');
ylabel('Signal [a.u.]','FontSize',14);
xlabel('flip angle [deg]','FontSize',14);
set(gca,'XTick',[10:10:90],'TickDir','out')
set(gca,'FontSize',14);
legend('\psi = 50°','\psi = 115.4°','\psi = 117°','\psi = 150°','\psi = 169°','Ernst Ampl.','Location','southwest');
lgd = legend;
lgd.FontSize = 14;

%compute the P-values for H20+CuSO4
flip_ind=[1:18];
psi_vec=[50 115.4 117 150 169];
epsilon=zeros(1,length(psi_vec));
a=[sig50(:,no_pulses) sig115(:,no_pulses) sig117(:,no_pulses) sig150(:,no_pulses) sig169(:,no_pulses)];
b=repmat(ernst',1,length(psi_vec));
ea=acosd(exp(-TR/T1)); %Ernst angle - only for convenience, actually it is not used further
sig_ea=sind(ea)*(1-exp(-TR/T1))/(1-exp(-TR/T1)*cosd(ea)); %signal at the Ernst angle - only for convenience, actually it is not used further
diffs=abs(a(flip_ind,:)-b(flip_ind,:))./ernst'; 
epsilon=100*sum(diffs,1)/18;

fprintf('\nP-values for H2O+CuSO4 in percent:');
fprintf('\n50°: %f  115.4°: %f  117°: %f  150°: %f  169°: %f\n',epsilon(find(psi_vec==50)),epsilon(find(10*psi_vec==1154)),epsilon(find(psi_vec==117)),epsilon(find(psi_vec==150)),epsilon(find(psi_vec==169)));

%compute the P-values for H20+CuSO4 without diffusion
flip_ind=[1:18];
psi_vec=[50 115.4 117 150 169];
epsilon=zeros(1,length(psi_vec));
a=[sig50D0(:,no_pulses) sig115D0(:,no_pulses) sig117D0(:,no_pulses) sig150D0(:,no_pulses) sig169D0(:,no_pulses)];
b=repmat(ernst',1,length(psi_vec));
ea=acosd(exp(-TR/T1)); %Ernst angle - only for convenience, actually it is not used further
sig_ea=sind(ea)*(1-exp(-TR/T1))/(1-exp(-TR/T1)*cosd(ea)); %signal at the Ernst angle - only for convenience, actually it is not used further
diffs=abs(a(flip_ind,:)-b(flip_ind,:))./ernst'; 
epsilon=100*sum(diffs,1)/18;

fprintf('\nP-values for H2O+CuSO4 without diffusion in percent:');
fprintf('\n50°: %f  115.4°: %f  117°: %f  150°: %f  169°: %f\n',epsilon(find(psi_vec==50)),epsilon(find(10*psi_vec==1154)),epsilon(find(psi_vec==117)),epsilon(find(psi_vec==150)),epsilon(find(psi_vec==169)));
