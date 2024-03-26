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

%snr data (normalized) per flip angle for silicone oil
snr50= [0.7995    1.0000    0.9065    0.7319    0.6056    0.5134    0.4214    0.3559    0.3040    0.2688    0.2219    0.2090    0.1965 0.1851    0.1726    0.1538    0.1551    0.1589];
snr115=[0.8151    1.0000    0.8960    0.7491    0.6342    0.5235    0.4442    0.4078    0.3739    0.3381    0.3093    0.2794    0.2519 0.2338    0.2011    0.1810    0.1521    0.1322];
snr117=[0.7941    1.0000    0.8916    0.7293    0.6471    0.5837    0.5575    0.5244    0.4874    0.4612    0.4339    0.3953    0.3565 0.3057    0.2769    0.2493    0.2078    0.1842];
snr150=[0.8049    1.0000    0.9146    0.7951    0.6839    0.5858    0.4703    0.4198    0.3507    0.2854    0.2399    0.2073    0.1754 0.1630    0.1401    0.1289    0.1154    0.1049];
snr169=[0.8016    1.0000    0.8965    0.7520    0.6182    0.5103    0.4354    0.3764    0.3283    0.3023    0.2798    0.2655    0.2462 0.2313    0.2121    0.1981    0.1770    0.1669];

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
set(gca,'XTick',[10:10:90],'TickDir','out');
set(gca,'FontSize',14);
axis([0 90 0 1]);
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

fprintf('\nepsilon-values for silicone oil in percent:');
fprintf('\n50°: %f  115.4°: %f  117°: %f  150°: %f  169°: %f\n',epsilon(find(psi_vec==50)),epsilon(find(10*psi_vec==1154)),epsilon(find(psi_vec==117)),epsilon(find(psi_vec==150)),epsilon(find(psi_vec==169)));


%% now simulate H20+CuSO4. Note the variables from above (i.e. those used for the silicone) are overwritten.

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

%snr data (normalized) for H20 + CuSO4
snr50= [0.5981    0.9215    1.0000    0.9656    0.8736    0.7827    0.6979    0.6271    0.5549    0.4991    0.4500    0.4123    0.3812    0.3545    0.3251    0.3057    0.2918    0.2716];
snr115=[0.5838    0.9162    1.0000    0.9592    0.8754    0.7932    0.7043    0.6349    0.5593    0.4976    0.4484    0.4035    0.3594    0.3225    0.2907    0.2562    0.2316    0.2113];
snr117=[0.5828    0.9117    1.0000    0.9643    0.8862    0.7939    0.7099    0.6315    0.5599    0.4968    0.4476    0.3948    0.3628    0.3215    0.2865    0.2582    0.2351    0.2088];
snr150=[0.5873    0.9090    1.0000    0.9567    0.8725    0.7851    0.7077    0.6258    0.5602    0.4943    0.4495    0.4016    0.3608    0.3263    0.2944    0.2647    0.2378    0.2128];
snr169=[0.5930    0.9106    1.0000    0.9673    0.8933    0.8074    0.7199    0.6532    0.5897    0.5298    0.4858    0.4396    0.4019    0.3654    0.3302    0.3001    0.2764    0.2479];

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
set(gca,'XTick',[10:10:90],'TickDir','out');
axis([0 90 0 1]);
set(gca,'FontSize',14);
legend('\psi = 50°','\psi = 115.4°','\psi = 117°','\psi = 150°','\psi = 169°','Ernst Ampl.','Location','southwest');
lgd = legend;
lgd.FontSize = 14;

%compute the epsilon-values for H20+CuSO4
flip_ind=[1:18];
psi_vec=[50 115.4 117 150 169];
epsilon=zeros(1,length(psi_vec));
a=[sig50(:,no_pulses) sig115(:,no_pulses) sig117(:,no_pulses) sig150(:,no_pulses) sig169(:,no_pulses)];
b=repmat(ernst',1,length(psi_vec));
ea=acosd(exp(-TR/T1)); %Ernst angle - only for convenience, actually it is not used further
sig_ea=sind(ea)*(1-exp(-TR/T1))/(1-exp(-TR/T1)*cosd(ea)); %signal at the Ernst angle - only for convenience, actually it is not used further
diffs=abs(a(flip_ind,:)-b(flip_ind,:))./ernst'; 
epsilon=100*sum(diffs,1)/18;

fprintf('\nepsilon-values for H2O+CuSO4 in percent:');
fprintf('\n50°: %f  115.4°: %f  117°: %f  150°: %f  169°: %f\n',epsilon(find(psi_vec==50)),epsilon(find(10*psi_vec==1154)),epsilon(find(psi_vec==117)),epsilon(find(psi_vec==150)),epsilon(find(psi_vec==169)));

%compute the epsilon-values for H20+CuSO4 without diffusion
flip_ind=[1:18];
psi_vec=[50 115.4 117 150 169];
epsilon=zeros(1,length(psi_vec));
a=[sig50D0(:,no_pulses) sig115D0(:,no_pulses) sig117D0(:,no_pulses) sig150D0(:,no_pulses) sig169D0(:,no_pulses)];
b=repmat(ernst',1,length(psi_vec));
ea=acosd(exp(-TR/T1)); %Ernst angle - only for convenience, actually it is not used further
sig_ea=sind(ea)*(1-exp(-TR/T1))/(1-exp(-TR/T1)*cosd(ea)); %signal at the Ernst angle - only for convenience, actually it is not used further
diffs=abs(a(flip_ind,:)-b(flip_ind,:))./ernst'; 
epsilon=100*sum(diffs,1)/18;

fprintf('\nepsilon-values for H2O+CuSO4 without diffusion in percent:');
fprintf('\n50°: %f  115.4°: %f  117°: %f  150°: %f  169°: %f\n',epsilon(find(psi_vec==50)),epsilon(find(10*psi_vec==1154)),epsilon(find(psi_vec==117)),epsilon(find(psi_vec==150)),epsilon(find(psi_vec==169)));
