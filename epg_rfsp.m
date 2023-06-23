function [S_plus, S_minus] = epg_rfps(alpha,np,T1,T2,TR,D,psi,res)

%This function calculates the steady state magnitude of an RF-spoiled gradient echo sequence under consideration
%of free and isotropic difusion. The used EPG-Algorithm is described in Weigel et al., JMR 205:276–285(2010)

%Input parameters:
% alpha - flip angle in deg
% np - number of pulses (i.e. sequence intervals)
% T1 - T1 in ms
% T2 - T2 in ms
% TR - TR in ms
% D - diffusion coefficient in mm^2/s. E.g., for free water at room
%     temperature use 2.0e-3
% psi - phase increment for RF-spoiling in deg
% resolution factor (see below)

% The following gradient duration d* and Amplitude A* parameters per TR are used for simulation of the EPG
% for and RF-spoiled GRE sequence with consideration of free and isotropic diffusion
% - the delay d3 might be somewhat uncommon but is a peculiarity of the Bruker implementaion
% - in slice direction, only the spoiler gradient was considered
% - gradients in phase encoding direction were neglected
%
% readout gradient:
%                                             A3 |------|
%                          A2 |------------------|      |
% |                           |                  |      |            | 
% |------d1-----|  d1  |--d3--|        d4        |  d5  |-----d6-----|
%               |      |
%               |      |
%               |------| A1
%               
% slice gradient: 
%                                               A4 |---d8---|           
% |                                                |        |
% |----------------------d7------------------------|        |---d9---|
% 

% constants
PI    = pi;
GAMMA = 2.675153362E8;   % in 1/s*1/T 

% initial calculations - in the following, SI units are used
TR = TR*1E-3;               % SI: s
T1 = T1*1E-3;               % SI: s
T2 = T2*1E-3;               % SI: s
alpha = alpha*PI/180.0;     % rad
psi = psi*PI/180.0;         % rad
e1 = exp(-TR/T1);
e2 = exp(-TR/T2);
D = D*1E-6;                 % SI: m^2/s

ampli=0.672;
% 0.672 mT/m as maximum amplitude of the used Bruker Biospec 70/20 scanner.
% The following Amplitudes A* are in percent of this maximum amplitude.

ampli=ampli*res;
% for simulation of change in in resolution. A factor res is multiplied to
% both readout and slice gradient amplitudes. E.g., when doubling
% resolution form 300um to 150um, we have res=2.
 
% gradient durations for readout gradients                
d1 = 0.00062;
d2 = 0.0021;
d3 = 0.0002;
d4 = 0.002;
d5 = 0.00437;
d6 = TR-d1-d2-d3-d4-d5;

%readout gradient amplitudes
A1 = -0.0297*ampli;
A2 = 0.0582*ampli;
A3 = 0.08*ampli;

% slice gradient durations, spoiler only
d8 = 0.0028;
d9 = 0.00062;
d7 = TR-d8-d9;

% slice gradient spoiler amplitude
A4 = 0.06*ampli;

% check if TR is long enough for total gradient duration
if (d6<0 || d7<0)
    error('invalid duration settings')
end

% dephasing parameter per readout gradient
dk1 = d2*A1*GAMMA;
dk2 = d4*A2*GAMMA;
dk3 = d5*A3*GAMMA;
dk=dk1+dk2+dk3; %delta_k between different EPG states for RO gradients
dkdk  = dk^2;

%... and per slice gradient
dks = d8*A4*GAMMA; %delta_k between different EPG states for slice gradient
dksdks = dks^2;

%allocate memory
fx=zeros(1,2*np+1);
fy=zeros(1,2*np+1);
pfx=zeros(1,2*np);
pfy=zeros(1,2*np);
pzx=zeros(1,np);
pzy=zeros(1,np);
S_plus=zeros(1,np);
S_minus=zeros(1,np);
zx=zeros(1,np);
zy=zeros(1,np);

% define M0
zx(0+1) = 1.0; 

%set pre pulse signal to zero before first pulse
S_minus(1)=0;

for l = 0:(np-1)   % RF pulse loop
    
    if (l>0)
    S_minus(l+1) = (fx(np+1)+1i*fy(np+1))*exp(-1i*phi);  %function output: the pre pulse signal
                                                    %the exp factor corrects for pulse phase as the ADC does it    
    end
    
    phi = l*(l+1)/2*psi;    % RF pulse phase for RF spoiling
    a = (cos(alpha/2.0))^2;
    b = (sin(alpha/2.0))^2;
    c = sin(alpha);
    d = cos(alpha);
    e = sin(phi);
    f = cos(phi);
    g = sin(2.0*phi);
    h = cos(2.0*phi);
    hb = h*b; gb = g*b; ec = e*c; fc = f*c;
    
    for k = 0:l     % EPG states: application of RF rotations (T)
        n = np+k; m = np-k;
        % post pulse states 
        pfx(n+1) =    a*fx(n+1) + hb*fx(m+1) + gb*fy(m+1) + ec*zx(k+1) + fc*zy(k+1);    % dephasing F real
        pfy(n+1) =    a*fy(n+1) - hb*fy(m+1) + gb*fx(m+1) - fc*zx(k+1) + ec*zy(k+1);    % dephasing F imag
        pfx(m+1) =   hb*fx(n+1) + gb*fy(n+1) +  a*fx(m+1) + ec*zx(k+1) - fc*zy(k+1);    % refoc F real
        pfy(m+1) =   gb*fx(n+1) - hb*fy(n+1) +  a*fy(m+1) - fc*zx(k+1) - ec*zy(k+1);    % refoc F imag
        pzx(k+1) = (-ec*fx(n+1) + fc*fy(n+1) - ec*fx(m+1) + fc*fy(m+1) + 2.0*d*zx(k+1)) / 2.0;  % long Z real
        pzy(k+1) = (-fc*fx(n+1) - ec*fy(n+1) + fc*fx(m+1) + ec*fy(m+1) + 2.0*d*zy(k+1)) / 2.0;  % long Z imag
    end
    
    S_plus(l+1) = (pfx(np+1)+1i*pfy(np+1))*exp(-1i*phi); %function output: the post pulse signal
                                                         %the exp factor corrects for pulse phase as the ADC does it    
    
    if (D>0) % diffusion weighting (D)
        for k = (-l):l  
            n = np+k;   
            if (k>=0)               
               %readout (RO) gradient
               %delay d1 between excitation and dephasing RO gradient
               pfx(n+1) = pfx(n+1)*exp(-k^2*dkdk*d1*D);
               pfy(n+1) = pfy(n+1)*exp(-k^2*dkdk*d1*D);
               %dephasing RO gradient with duration d2 and Amplitude A1
               pfx(n+1) = pfx(n+1)*exp(-(1/3)*(k*k*dkdk + k*dk*(k*dk+dk1) + (k*dk+dk1)*(k*dk+dk1) )*d2*D);
               pfy(n+1) = pfy(n+1)*exp(-(1/3)*(k*k*dkdk + k*dk*(k*dk+dk1) + (k*dk+dk1)*(k*dk+dk1) )*d2*D);
               %Delay d3 after dephasing gradient and readout gradient (this is a peculiarity of the Bruker Paravision implemantation) 
               pfx(n+1) = pfx(n+1)*exp(-(k*dk+dk1)*(k*dk+dk1)*d3*D);
               pfy(n+1) = pfy(n+1)*exp(-(k*dk+dk1)*(k*dk+dk1)*d3*D);
               %RO gradient with duration d4 and amplitude A2
               pfx(n+1) = pfx(n+1)*exp(-(1/3)*((k*dk+dk1)*(k*dk+dk1) + (k*dk+dk1)*(k*dk+dk1+dk2) + (k*dk+dk1+dk2)*(k*dk+dk1+dk2) )*d4*D);
               pfy(n+1) = pfy(n+1)*exp(-(1/3)*((k*dk+dk1)*(k*dk+dk1) + (k*dk+dk1)*(k*dk+dk1+dk2) + (k*dk+dk1+dk2)*(k*dk+dk1+dk2) )*d4*D);
               %spoiler gradient with duration d5 and amplitude A3
               pfx(n+1) = pfx(n+1)*exp(-(1/3)*( (k*dk+dk1+dk2)*(k*dk+dk1+dk2) +  (k*dk+dk1+dk2)*(k*dk+dk) +  (k*dk+dk)*(k*dk+dk) )*d5*D);
               pfy(n+1) = pfy(n+1)*exp(-(1/3)*( (k*dk+dk1+dk2)*(k*dk+dk1+dk2) +  (k*dk+dk1+dk2)*(k*dk+dk) +  (k*dk+dk)*(k*dk+dk) )*d5*D);
               %Delay d6 for TR filling
               pfx(n+1) = pfx(n+1)*exp(-(k*dk+dk)*(k*dk+dk)*d6*D);
               pfy(n+1) = pfy(n+1)*exp(-(k*dk+dk)*(k*dk+dk)*d6*D);
               %slice gradient
               %delay d7 before slice spoiler gradient
               pfx(n+1) = pfx(n+1)*exp(-k^2*dksdks*d7*D);
               pfy(n+1) = pfy(n+1)*exp(-k^2*dksdks*d7*D);
               %spoiler gradient in slice direction with duration d8 and amplitude A4
               pfx(n+1) = pfx(n+1)*exp(-(1/3)*(k*k*dksdks + k*dks*(k*dks+dks) + (k*dks+dks)*(k*dks+dks) )*d8*D);
               pfy(n+1) = pfy(n+1)*exp(-(1/3)*(k*k*dksdks + k*dks*(k*dks+dks) + (k*dks+dks)*(k*dks+dks) )*d8*D);
               %delay d9 for TR filling
               pfx(n+1) = pfx(n+1)*exp(-(k*dks+dks)*(k*dks+dks)*d9*D);
               pfy(n+1) = pfy(n+1)*exp(-(k*dks+dks)*(k*dks+dks)*d9*D);           
               %diffusion weighting of z-states
               pzx(k+1) = pzx(k+1)*exp(-k^2*dkdk*TR*D);
               pzy(k+1) = pzy(k+1)*exp(-k^2*dkdk*TR*D);
               pzx(k+1) = pzx(k+1)*exp(-k^2*dksdks*TR*D);
               pzy(k+1) = pzy(k+1)*exp(-k^2*dksdks*TR*D);
            end
            
            %diffusion weighting of k<0 states
            if (k<0)
               %delay d1 between excitation and dephasing RO gradient
               pfx(n+1) = pfx(n+1)*exp(-k^2*dkdk*d1*D);
               pfy(n+1) = pfy(n+1)*exp(-k^2*dkdk*d1*D);
               %dephasing RO gradient with duration d2 and Amplitude A1
               pfx(n+1) = pfx(n+1)*exp(-(1/3)*(k*k*dkdk + k*dk*(k*dk+dk1) + (k*dk+dk1)*(k*dk+dk1) )*d2*D);
               pfy(n+1) = pfy(n+1)*exp(-(1/3)*(k*k*dkdk + k*dk*(k*dk+dk1) + (k*dk+dk1)*(k*dk+dk1) )*d2*D);
               %Delay d3 after dephasing gradient and readout gradient (this is a peculiarity of the Bruker Paravision implemantation) 
               pfx(n+1) = pfx(n+1)*exp(-(k*dk+dk1)*(k*dk+dk1)*d3*D);
               pfy(n+1) = pfy(n+1)*exp(-(k*dk+dk1)*(k*dk+dk1)*d3*D); 
               %RO gradient with duration d4 and amplitude A2
               pfx(n+1) = pfx(n+1)*exp(-(1/3)*((k*dk+dk1)*(k*dk+dk1) + (k*dk+dk1)*(k*dk+dk1+dk2) + (k*dk+dk1+dk2)*(k*dk+dk1+dk2) )*d4*D);
               pfy(n+1) = pfy(n+1)*exp(-(1/3)*((k*dk+dk1)*(k*dk+dk1) + (k*dk+dk1)*(k*dk+dk1+dk2) + (k*dk+dk1+dk2)*(k*dk+dk1+dk2) )*d4*D);
               %spoiler gradient with duration d5 and amplitude A3
               pfx(n+1) = pfx(n+1)*exp(-(1/3)*( (k*dk+dk1+dk2)*(k*dk+dk1+dk2) +  (k*dk+dk1+dk2)*(k*dk+dk) +  (k*dk+dk)*(k*dk+dk) )*d5*D);
               pfy(n+1) = pfy(n+1)*exp(-(1/3)*( (k*dk+dk1+dk2)*(k*dk+dk1+dk2) +  (k*dk+dk1+dk2)*(k*dk+dk) +  (k*dk+dk)*(k*dk+dk) )*d5*D);
               %Delay d6 for TR filling
               pfx(n+1) = pfx(n+1)*exp(-(k*dk+dk)*(k*dk+dk)*d6*D);
               pfy(n+1) = pfy(n+1)*exp(-(k*dk+dk)*(k*dk+dk)*d6*D);
               %slice grad
               %delay d7 before slice spoiler gradien
               pfx(n+1) = pfx(n+1)*exp(-k^2*dksdks*d7*D);
               pfy(n+1) = pfy(n+1)*exp(-k^2*dksdks*d7*D);
               %spoiler gradient in slice direction with duration d8 and amplitude A4
               pfx(n+1) = pfx(n+1)*exp(-(1/3)*(k*k*dksdks + k*dks*(k*dks+dks) + (k*dks+dks)*(k*dks+dks) )*d8*D);
               pfy(n+1) = pfy(n+1)*exp(-(1/3)*(k*k*dksdks + k*dks*(k*dks+dks) + (k*dks+dks)*(k*dks+dks) )*d8*D);
               %delay d9 for TR filling
               pfx(n+1) = pfx(n+1)*exp(-(k*dks+dks)*(k*dks+dks)*d9*D);
               pfy(n+1) = pfy(n+1)*exp(-(k*dks+dks)*(k*dks+dks)*d9*D);
            end
        end
    end %eof (D>0)

    for k = (-l):l  % gradient evolution (S) and variable swap and relaxation (E) -
        n = np+k;   % separated from diffusion to avoid conflicts with phase history: see Weigel et al. JMR 2010
        fx(n+1+1) = pfx(n+1)*e2;   
        fy(n+1+1) = pfy(n+1)*e2;
        
        if (k>0)
            zx(k+1) = pzx(k+1)*e1;
            zy(k+1) = pzy(k+1)*e1;
        end

        if (k==0)
            zx(k+1) = pzx(k+1)*e1 + 1.0 - e1;   % T1 recovery
            zy(k+1) = pzy(k+1)*e1;
        end
    end
end %eof RF pulse loop

end %eof function





