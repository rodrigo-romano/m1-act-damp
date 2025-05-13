%
% M1 (F&M) outer loop analysis
%

%clear all;

%% General settings

% Number of frequency response points
Nom = 1000;
% Frequency points (rad/s)
w = logspace(log10(0.1),log10(75),Nom)*2*pi;
% Figure number to avoid overwrite
figNumber = 2000;
% M1 segment
seg = 1;
% Flag to load static model data
load_static_model = true; %false;%
% - - - 
% Frequency response plot options
% - - - 
% Bode plot
hbode = bodeoptions;
hbode.TickLabel.FontSize = 12;
hbode.YLabel.FontSize = 12;
hbode.XLabel.FontSize = 12;
hbode.TickLabel.FontSize = 10;
hbode.FreqUnits = 'Hz';
hbode.Grid = 'on';
hbode.PhaseMatching = 'on';
hbode.PhaseMatchingValue = -180;
% Nichols plot
hnichols = nicholsoptions;
hnichols.TickLabel.FontSize = 12;
hnichols.YLabel.FontSize = 12;
hnichols.XLabel.FontSize = 12;
hnichols.FreqUnits = 'Hz';
hnichols.Grid = 'on';
hnichols.PhaseMatching = 'on';
hnichols.PhaseMatchingValue = -180;
% - - - 


%% Load M1 controller and interface model parameters
%%

% load('controls_5pt1e_z30_dev','m1sys','fem')
% load('controls_5pt1f_z30_iTT_oad','m1sys','fem')
% load('../versions/controls_5pt1f_z30_llTT_oad','m1sys','fem')
load('../../gmt-ims/controls_5pt1g_z30_llTT_oad','m1sys','fem')
% load('../controls_5pt1g1K_z30_llTTsoftHP_oad.mat','m1sys','fem')


% m1sys{seg}.HPstiff: Hardpoint stiffness
% m1sys{seg}.LC2CG: HP: F&M convertion matrix

% Overwrite controllers
if true
    warning('Overwriting M1 controller model parameters!') %#ok<*UNRCH>
    m1sys = ovr_ofl_crtl(m1sys,seg);
end

m1_act_damp = 1800; %[Ns/m] Actuator damping
m1_dt_folder = '/Users/rromano/Documents/MATLAB/M1';
if(seg < 7)
    load(fullfile(m1_dt_folder,'OA_SupportActuatorArrayConfig'),'OA_Upsilon');
    Upsilon = OA_Upsilon;
else
    load(fullfile(m1_dt_folder,'CS_SupportActuatorArrayConfig'),'CS_Upsilon');
    Upsilon = CS_Upsilon;
end



%% Load FE plant SS model

% ===>>> User shall check file path and name
% FileFolder = fullfile(im.lfFolder,"20240403_1556_zen_30_M1_202110_FSM_202305_Mount_202305_IDOM_concreteAndFoundation_finalSI");
FileFolder = fullfile(im.lfFolder,"20250430_1659_zen_30_M1_202110_FSM_202305_Mount_202305_pier_202411_M1_actDamping");

if(~exist('inputTable','var') || 0)
    % Start loading modal model parameters
    try
        FileName = "modal_state_space_model_2ndOrder.mat";
        load(fullfile(FileFolder,FileName),'inputs2ModalF','modalDisp2Outputs',...
            'eigenfrequencies','proportionalDampingVec','inputTable','outputTable');
        % Handle modal parameters
        om2 = (2*pi*eigenfrequencies(:)).^2;
        twice_zom = 2*proportionalDampingVec(:).*(2*pi*eigenfrequencies(:));
        
        % State-space model matrices
        B = [zeros(size(inputs2ModalF));inputs2ModalF];
        C = [modalDisp2Outputs,zeros(size(modalDisp2Outputs))];
        n_m = size(inputs2ModalF,1);
        A = [zeros(n_m),eye(n_m);...
            -diag(om2), -diag(twice_zom)];
        
    % In case modal model is not available, try to load state-space model matrices    
    catch
        FileName = "modal_state_space_model.mat";     
        load(fullfile(FileFolder,FileName),'A','B','C','D','inputTable','outputTable');
        % Retrieve modal-form parameters
        [aux,temp] = spdiags(full(A));      % Indexes to get
        om2 = -aux(1:temp(end),1);          % Vector of squared eigenfrequencies
        twice_zom = -aux(temp(end)+1:end,2);    % Vector of 2*damp*om
    end
    
    if true
        % Remove the first 3 modes and adjusts SS model matrices
        [aux,temp] = spdiags(full(A));
        om2 = om2(1+3:end);
        twice_zom = twice_zom(1+3:end);
        nm = length(om2);
        A = [zeros(nm),eye(nm);-diag(om2), -diag(twice_zom)];
        B = [zeros(nm,size(B,2));B(temp(end)+1+3:end,:)];
        C = [C(:,1+3:nm+3),zeros(size(C,1),nm)];
    end
    
    if(load_static_model)
        try
            staticSolFile = fullfile(FileFolder,"static_reduction_model.mat");
            try
                load(staticSolFile,'gainMatrixMountControlled');
                gainMatrix = gainMatrixMountControlled;
            catch
                load(staticSolFile,'gainMatrix');
            end
        catch
            load_static_model = 0;
            warning('Unable to load static gain matrix\n')
        end
    end
end

om0 = sqrt(om2);
damp = 0.5 * twice_zom./sqrt(om2);
i_p0 = find(isnan(damp));
if(i_p0), damp(i_p0) = damp(i_p0(end)+1)*ones(size(i_p0)); end
        
fprintf('Model from %s\n loaded.\n', fullfile(FileFolder,FileName));



%% Extract M1 subsystem
%%
% OUTPUTS
out0 = outputTable{sprintf('M1_actuators_segment_%d_D', seg),"indices"}{1}(:);
out1a = outputTable{'OSS_Hardpoint_D',"indices"}{1}((seg-1)*12+(1:6)+6);  %face
out1b = outputTable{'OSS_Hardpoint_D',"indices"}{1}((seg-1)*12+(1:6));    %cell
out2 = outputTable{'OSS_Hardpoint_force',"indices"}{1}((seg-1)*6+(1:6));

c = [Upsilon*C(out0,:);...
    C(out1a,:)-C(out1b,:);...
    C(out2,:)];


% M1 SA force inputs
in = inputTable{sprintf('M1_actuators_segment_%d',seg),"indices"}{1};
b = B(:,in)*m1sys{seg}.Kbal;
fprintf('M1 OFL analysis using SA force inputs.\n');

if(load_static_model && exist('gainMatrix','var') && 1)
    d = [Upsilon*gainMatrix(out0,in);...
        gainMatrix(out1a,in)-gainMatrix(out1b,in);...
        gainMatrix(out2,in)]*m1sys{seg}.Kbal ...
        + c*(A\b);
else
    d = zeros(ny,nu);
    fprintf('FEM static solution not available. Residual D matrix not incorporated!\n')
end

ny = size(c,1);
nu = size(b,2);

% Pick selected submatrices related to M1 outer loop and wrap subsystem
% model into CMM object
comp_name = struct('name','M1ofl','shortname',[]);
sys_prime = struct('A',full(A),'B',b,'C',c,'D',d,'method',[]);
inputs{1} = struct('name','M1inputs','indexes',1:nu,'ctrace',[]);
outputs{1} = struct('name','M1outputs','indexes',1:ny);

objM1 = component_mode_model(comp_name,sys_prime,inputs,outputs);


%% Assess HP stiffness
%%
inHP_F_Ind = inputTable{'OSS_Harpoint_delta_F','indices'}{1}((seg-1)*6+(1:6));
HPstiffvec = 1./diag(...
        gainMatrix(out1a,inHP_F_Ind)-gainMatrix(out1b,inHP_F_Ind));
fprintf('HP stiffness:%.3g+/-%.2g [N/um]\n',...
    1e-6*mean(HPstiffvec),3e-6*std(HPstiffvec));
Hpk = mean(HPstiffvec);


%% Plant Frequency Responses
%%
% Structural plant submodel frequenct response
if(~exist('res','var') || 0)
    res = bode_second_order(objM1,w,om0,damp,1:nu,1:ny);
    res = frd(res,w,fem.Ts);
end

% Pneumatic actuator frequency response
H1_sa = freqresp(tf(10*2*pi,[1 10*2*pi])*tf(50*2*pi,[1 50*2*pi]), w);
SAdyn_Fr = frd(H1_sa, w, fem.Ts);

% Loop frequency response
feedin = 1:nu;
Z_ = zeros(length(out0));
Hd = freqresp(tf([1 0],1),w);
fb_sys = eye(6)*frd(Hd(:), w, fem.Ts)*m1_act_damp;
sys = feedback(res, fb_sys, feedin, 1:6, -1);
sys4 = feedback(res, fb_sys*4, feedin, 1:6, -1);

% lag term to represent pneumatic circuit nonlinearity
H_ = freqresp(0.2/1 *tf([1 1*2*pi],[1 0.2*2*pi]), w);

in_hp_lcD = 6+(1:6);
in_hp_lcF = 12+(1:6);
G0 = -m1sys{seg}.LC2CG* Hpk *res(in_hp_lcD,:)*SAdyn_Fr;
G = -m1sys{seg}.LC2CG* Hpk *sys(in_hp_lcD,:)*SAdyn_Fr;
% G_ = -m1sys{seg}.LC2CG* sys4(in_hp_lcF,:)*SAdyn_Fr;
G__ = -m1sys{seg}.LC2CG *Hpk *sys(in_hp_lcD,:) *frd(H_, w, fem.Ts) *SAdyn_Fr;


%% Open-loop frequency response magnitude
%%
% % Relevant for notch filter tuning

figure(figNumber*seg)
% Resize to report format
% set(gcf,'Position',[900   267   304*3/2   420*0.65]);
set(gcf,'Position',[500   267   1.5*304*3/2   420]);
label_ = {'F_x','F_y','F_z','M_x','M_y','M_z'};
for i1 = 1:6
    subplot(2,3,i1)
    [MagG,~] = bode(G(i1,i1),w);
    semilogx((w/2/pi)',20*log10(reshape(MagG,Nom,1,1)));
    hold on;
    [MagG,~] = bode(G0(i1,i1),w);
    semilogx((w/2/pi)',20*log10(reshape(MagG,Nom,1,1)));
%     [MagG,~] = bode(G_(i1,i1),w);
%     semilogx((w/2/pi)',20*log10(reshape(MagG,Nom,1,1)),'-.');
    [MagG,~] = bode(G__(i1,i1),w);
    semilogx((w/2/pi)',20*log10(reshape(MagG,Nom,1,1)),'--');
    hold off;
    if(i1 > 3), xlabel('Frequency (Hz)','fontsize',hbode.YLabel.FontSize); end
    if(mod(i1,3) == 1), ylabel('Magnitude (dB)','fontsize',hbode.YLabel.FontSize); end
    xlim([0.1,50]); grid on;
    title(label_{i1},'fontsize',hbode.YLabel.FontSize)
    
    if(i1==1)
        legend('M1 damping (1800Ns/m)',...
            'No M1 damping',...
            'M1 damping + lag filter',...
            'Location','southwest','fontsize',hbode.YLabel.FontSize-2);
        legend('boxoff');
    end
end


return
%% Loop caracterization: margins, sensitivity and bandwidth
%%
% Create M1 MIMO controller object
m1oflmimo = tf;	% Create null tf object
for ich = 1:6, m1oflmimo(ich,ich) = m1sys{seg}.ofl.SSdtC{ich}; end
% Controller frequency response
upsample_rate = m1sys{seg}.ofl.SSdtC{1}.Ts/fem.Ts;
m1MIMO_C_Fr = frd(upsample(m1oflmimo, upsample_rate),w);

GK = G__ * m1MIMO_C_Fr;

% Initialize variables
GM = zeros(6,1);
PM = zeros(6,1);
VM = zeros(6,1);
wpi = zeros(6,1);   % rad/s
w0db = zeros(6,1);  % rad/s
BW = zeros(6,1);
S = frd(zeros(6,1,Nom),w);
T = frd(zeros(6,1,Nom),w);

% Compute (SISO) Margins
for i1=1:6
    [GM(i1),PM(i1),wpi(i1),w0db(i1)] = margin(GK(i1,i1));   % Stability margins
    S(i1) = 1/(1+GK(i1,i1));        % SISO channel sensitivity function
    VM(i1) = 1/max(abs(S(i1).response(:)));
    T(i1) = GK(i1,i1)/(1+GK(i1,i1));    % SISO channel complementary sensibility function
    BW(i1) = bandwidth(T(i1));   % Channel bandwidth
end

fprintf("\nM1-S%d outer control loop (force):",seg);
header = " GM \t PM \t f_-pi \t f_0dB \t BW \t VM";
info = num2str([GM,PM,wpi/2/pi,w0db/2/pi,BW/2/pi,VM],'%.3g\t');
fprintf("\n" + header + "\n");
disp(info);     % Display margins


%% Control Loop Analysis plots
%%

%%
hTF = figure(figNumber*seg+1);
% Resize to report format
set(gcf,'Position',[910   350   304*3/2   420*1.2]);
bode(GK(1,1),GK(2,2),GK(3,3),GK(4,4),GK(5,5),GK(6,6),hbode);
hold on;
plot_w0_6db_bode(hTF,min(w0db)/2/pi)
title('M1 outer loop frequency response','fontsize',hbode.YLabel.FontSize)
legend('F_x','F_y','F_z','M_x','M_y','M_z',...
    'Location','southwest','Orientation','horizontal',...
    'NumColumns',3,'fontsize',hbode.YLabel.FontSize); legend('boxoff')
xlim([0.1, 50]);%xlim([min(w)/2/pi,max(w)/2/pi])
grid on; hold off;


%% %
figure(figNumber*seg+2)
set(gcf,'Position',[910   350   304*3/2   420]);
nicholsplot(GK(1,1),GK(2,2),GK(3,3),GK(4,4),GK(5,5),GK(6,6),hnichols);
title('M1 outer loop frequency response','fontsize',hbode.YLabel.FontSize)
hold on;
ylim([-60,25])
plot_nichols_ellipse(0.5)
legend('F_x','F_y','F_z','M_x','M_y','M_z',...
    'Location','northwest','Orientation','horizontal',...
    'NumColumns',2,'fontsize',hbode.YLabel.FontSize); legend('boxoff')

grid on; hold off;


% %
figure(figNumber*seg+3)

subplot(2,1,1)
for i1 = 1:6
    [MagS,~] = bode(S(i1),w);
    semilogx((w/2/pi)',20*log10(reshape(MagS,Nom,1,1)));
    hold on;
end
plot_S_envelope();
legend('F_x','F_y','F_z','M_x','M_y','M_z',...
    'Orientation','horizontal','Location','southeast',...
    'NumColumns',3,'fontsize',hbode.YLabel.FontSize); legend('boxoff')
grid on; hold off;
ylabel('Magnitude (dB)','fontsize',hbode.YLabel.FontSize)
title('Sensitivity functions','fontsize',hbode.YLabel.FontSize)
xlim([min(w)/2/pi,max(w)/2/pi])

subplot(2,1,2)
for i1 = 1:6
    [MagT,~] = bode(T(i1),w);
    semilogx((w/2/pi)',20*log10(reshape(MagT,Nom,1,1)));
    hold on;
end
xlim([min(w)/2/pi,max(w)/2/pi])
ylim([-60,6])

plot_cl_OAD_req(1,4); %max(BW/2/pi)
ylabel('Magnitude (dB)','fontsize',hbode.YLabel.FontSize)
title('Complementary sensitivity functions','fontsize',hbode.YLabel.FontSize)
xlabel('Frequency (Hz)','fontsize',hbode.XLabel.FontSize)
grid on; hold off;

% Resize to report format
set(gcf,'Position',[900   267   304*3/2   420*1.2]);



% %%
% SK = (eye(6)+GK) \ -eye(6);%-(m1MIMO_C_Fr);
% figure(figNumber*seg+4)
% sigmaSK = sigma(SK);
% 
% hsk = semilogx((w/2/pi)',20*log10(squeeze(sigmaSK(1,:))));
% hold on;
% semilogx((w/2/pi)',20*log10(squeeze(sigmaSK(6,:))),'color', hsk.Color);
% 
% % Resize to report format
% set(gcf,'Position',[900   267   304*3/2   420*1.2]);



%% Auxiliar functions
%%
function m1sys = ovr_ofl_crtl(m1sys,seg)  

% S-notch filter prototype. Tuning parameters: zn, zd, fd, fw
snotch = @(zn,zd,fd,fn) tf((fd/fn)^2*[1 2*zn*(2*pi*fn) (2*pi*fn)^2],...
    [1 2*zd*(2*pi*fd) (2*pi*fd)^2]); %#ok<*NASGU> 
notchF = @(fc,F,delta) tf([1, 4*pi*fc/(F*delta), 4*(pi*fc)^2],...
    [1, 4*pi*fc/F, 4*(pi*fc)^2]);
sndUDampF = @(fc,damp) zpk([],[-fc*2*pi*(damp + 1i*sqrt(1-damp^2)),...
    conj(-fc*2*pi*(damp + 1i*sqrt(1-damp^2)))],(fc*2*pi)^2);


% M1 outer force loop (OFL) controller matrix
ofl.Ts = 1/100;  % OFL sampling time
% % I-controller
% kc = 8;    % HP Stiffness change: 144N/um -> 69.5N/um
% fbH = kc * tf(1,[1 0]);

% PI+Lead-lag (SPIE2024-)
fPI = tf([0.1 5],[1 0]);
flead = tf([1 0.1*2*pi],[1 0.3*2*pi]);
flag = tf([1 15*2*pi],[1 5*2*pi]);
fbH = fPI * flead * flag;

% % Notch and low-pass filters of each channel
% % Channels := {1}:Fx; {2}:Fy; {3}:Fz; {4}:Mx; {5}:My; {6}:Mz
% % Outer axis (OA) segment
% ofl_filterOA{1} = 0.9 * notchF(11.5,1,20);
% ofl_filterOA{2} = ofl_filterOA{1};
% ofl_filterOA{3} = 1.2*notchF(23.5,1,20);
% ofl_filterOA{4} = notchF(20,1,6); %snotch(0.04,0.4,9.5,10.5)
% ofl_filterOA{5} = ofl_filterOA{4};
% ofl_filterOA{6} = notchF(18.5,1,20);
% % Center segment (CS)
% ofl_filterCS{1} = 0.9 * notchF(11.5,1,20);
% ofl_filterCS{2} = ofl_filterCS{1};
% ofl_filterCS{3} = 1.2*notchF(21.8,1,20);
% ofl_filterCS{4} = notchF(21,1,6); %snotch(0.04,0.4,9.5,10.5);
% ofl_filterCS{5} = ofl_filterCS{4};
% ofl_filterCS{6} = notchF(19,1,20);
% if(seg ~= 7), ofl_filter = ofl_filterOA;
% else, ofl_filter = ofl_filterCS;
% end

% OFL controller

% Controller channel loop
oflC_ss = cell(6,1);
plot_labels = {'F_x','F_y','F_z','M_x','M_y','M_z'};
for ich = 1:numel(oflC_ss)
    oflC_ss{seg} = balreal(ss(fbH));
    m1sys{seg}.ofl.SSdtC{ich} = c2d(oflC_ss{seg},ofl.Ts,'foh');
    
    hbode = bodeoptions;
    hbode.FreqUnits = 'Hz';
    hbode.XLabel.FontSize = 11;
    hbode.YLabel.FontSize = 11;
    hbode.Ylabel.String{1} = 'Mag';
    hbode.Title.FontSize = 12;
    figure(7000+seg)
    subplot(2,3,ich)
    bode(oflC_ss{seg},hbode); hold on;
    bode(m1sys{seg}.ofl.SSdtC{ich},'r--',hbode); hold off;
    grid on; title(sprintf('OFL %s Controller',plot_labels{ich}));
    xlim([0.05, 0.5/ofl.Ts])
end
set(gcf,'Position',[600   567   1.6*304*3/2   420*0.95]);
legend('CT','DT','Location','northwest'); legend box off;
end



%%
% Function to plot -6dB margin after crossover frequency
% (RTF).
function plot_w0_6db_bode(hTF,w0db)
children = get(hTF, 'Children');	% use this handle to obtain list of figure's children
magChild = children(end);           % Pick a handle to axes of magnitude in bode diagram
axes(magChild);                     % Make those axes current
ax = findobj(gcf,'type','axes');
xlim_ = get(ax(1),'XLim');
S_env_w = [w0db, w0db, xlim_(2)]';
S_env_mag = [0, -6, -6]';
semilogx(S_env_w,S_env_mag,'r--','Linewidth',2.5);

end

%%
function plot_S_envelope()

ax = findobj(gcf,'type','axes');
xlim_ = get(ax(1),'XLim');
S_env_w = [0.1, 1, 2, xlim_(2)]';
S_env_mag = [-20, 0, 6, 6]';
semilogx(S_env_w,S_env_mag,'r--','Linewidth',2.5);

end


%%
function plot_cl_OAD_req(CL_Bw0,CL_Bw)
ax = findobj(gcf,'type','axes');
xlim_ = get(ax(2),'XLim');
ylim_ = get(ax(2),'YLim');

semilogx([xlim_(1) CL_Bw0],[-3 -3],'r--','Linewidth',2.5);
semilogx([CL_Bw0 CL_Bw0],[-100 -3],'r--','Linewidth',2.5);
semilogx([CL_Bw xlim_(2)],[-6 -6],'r--','Linewidth',2.5);
semilogx([CL_Bw CL_Bw],[-6 2*ylim_(2)],'r--','Linewidth',2.5);
end


%% Function to plot nichols robustness boundaries
function plot_nichols_ellipse(VM)
bcolor = [0.9 0.0 0.0];

ax = findobj(gcf,'type','axes');
xlim_ = get(ax(1),'XLim');
ylim_ = get(ax(1),'YLim');

thetaS = linspace(-0.8*pi/2,0,501)';
thetaT = linspace(0,0.9*pi/2,501)';

y = (1/VM)*exp(1i*thetaS);

S = y./(1-y);
ph=angle(S)*(180/pi); ph=ph-(360*ceil(ph/360)); mag=20*log10(abs(S));
line(ph,mag,'color',bcolor,'LineStyle','--','Linewidth',2);

y = (1/VM)*exp(1i*thetaT);
T = (1-y)./y;
ph=angle(T)*(180/pi); ph=ph-(360*ceil(ph/360)); mag=20*log10(abs(T));
line(ph,mag,'color',bcolor,'LineStyle','--','Linewidth',2.5);

if 1
    GM = 6; %PM = 30;
    line([xlim_(1) -180],GM*[-1 -1],'color',bcolor,'LineStyle','--','Linewidth',2.5)
    line([-180 -180],[GM ylim_(2)],'color',bcolor,'LineStyle','--','Linewidth',2.5)
end

text(0.9*xlim_(1),-GM,sprintf('robustness\nboundary'),...
    'HorizontalAlignment','left',...
    'VerticalAlignment','bottom',...
    'color',bcolor,'Fontsize',12);

end