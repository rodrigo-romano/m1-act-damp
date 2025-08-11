function [M_oa_xyz2abc_vel, M_cs_xyz2abc_vel, M_oa_abc2xyz, M_cs_abc2xyz] =...
    calc_ort2cylTs(m1_dt_folder)

%% Function to calculate M1 actuator transformations
% function [M_oa_xyz2abc_vel, M_cs_xyz2abc_vel, M_oa_abc2xyz, M_cs_abc2xyz] = calc_ort2cylTs(m1_dt_folder)

%%
% GMT-REF-05541 / Equation 4-2: ()
tmat = @(th, phia, phib, phic) [cos(phia)*sin(th), sin(phia)*sin(th), cos(th);...
                                cos(phib)*sin(th), sin(phib)*sin(th), cos(th);...
                                cos(phic)*sin(th), sin(phic)*sin(th), cos(th)];
th = 40*pi/180;
phia = 0*pi/180;
phib = +120*pi/180;
phic = -120*pi/180;
xyz2abc_vel40 = tmat(th, phia, phib, phic);

phia = +180*pi/180;
phib = -60*pi/180;
phic = +60*pi/180;
xyz2abc_vel41 = tmat(th, phia, phib, phic);
% OA
load(fullfile(m1_dt_folder,'OA_SupportActuatorArrayConfig'), 'OA_ActData');

nr = numel(find(OA_ActData(:,5) >= 40))*3 +...
    numel(find(OA_ActData(:,5) > 5 & (OA_ActData(:,5) < 40))) +...
    numel(find(OA_ActData(:,5) == 5))*6 +...
    numel(find(OA_ActData(:,5) < 5));
nc = nr - numel(find(OA_ActData(:,5) == 5))*3;
fprintf("Dimension of Txyz2abc for OA segment: %dx%d\n",nr,nc);
M_oa_xyz2abc_vel = zeros(nr,nc);
M_oa_abc2xyzT = zeros(nr,nc);
i_cols = 1;
i_rows = 1;

for i_act = 1:size(OA_ActData,1)

    if(OA_ActData(i_act,5) == 40)
        M_oa_xyz2abc_vel((0:2)+i_rows, (0:2)+i_cols) = xyz2abc_vel40;
        M_oa_abc2xyzT((0:2)+i_rows, (0:2)+i_cols) = xyz2abc_vel40;
        i_cols = i_cols+3;
        i_rows = i_rows+3;
    elseif(OA_ActData(i_act,5) == 41)
        M_oa_xyz2abc_vel((0:2)+i_rows, (0:2)+i_cols) = xyz2abc_vel41;
        M_oa_abc2xyzT((0:2)+i_rows, (0:2)+i_cols) = xyz2abc_vel41;
        i_cols = i_cols+3;
        i_rows = i_rows+3;
    elseif(OA_ActData(i_act,5) == 5)
        M_oa_xyz2abc_vel((0:5)+i_rows, (0:2)+i_cols) = [xyz2abc_vel40;xyz2abc_vel41];
        M_oa_abc2xyzT((0:5)+i_rows, (0:2)+i_cols) = [xyz2abc_vel40;xyz2abc_vel41];
        i_cols = i_cols+3;
        i_rows = i_rows+6;
    else
        M_oa_xyz2abc_vel(i_rows,i_cols) = 1;
        M_oa_abc2xyzT(i_rows,i_cols) = 1;
        i_cols = i_cols+1;
        i_rows = i_rows+1;
    end
end

M_oa_abc2xyz = M_oa_abc2xyzT';

% CS
load(fullfile(m1_dt_folder,'CS_SupportActuatorArrayConfig'),'CS_ActData');

nr = numel(find(CS_ActData(:,5) >= 40))*3 +...
    numel(find(CS_ActData(:,5) > 5 & (CS_ActData(:,5) < 40))) +...
    numel(find(CS_ActData(:,5) == 5))*6 +...
    numel(find(CS_ActData(:,5) < 5));
nc = nr - numel(find(CS_ActData(:,5) == 5))*3;
fprintf("Dimension of Txyz2abc for CS segment: %dx%d\n",nr,nc);
M_cs_xyz2abc_vel = zeros(nr,nc);
M_cs_abc2xyzT = zeros(nr,nc);
i_cols = 1;
i_rows = 1;
for i_act = 1:size(CS_ActData,1)

    if(CS_ActData(i_act,5) == 40)
        M_cs_xyz2abc_vel((0:2)+i_rows, (0:2)+i_cols) = xyz2abc_vel40;
        M_cs_abc2xyzT((0:2)+i_rows, (0:2)+i_cols) = xyz2abc_vel40;
        i_cols = i_cols+3;
        i_rows = i_rows+3;
    elseif(CS_ActData(i_act,5) == 41)
        M_cs_xyz2abc_vel((0:2)+i_rows, (0:2)+i_cols) = xyz2abc_vel41;
        M_cs_abc2xyzT((0:2)+i_rows, (0:2)+i_cols) = xyz2abc_vel41;
        i_cols = i_cols+3;
        i_rows = i_rows+3;
    elseif(CS_ActData(i_act,5) == 5)
        M_cs_xyz2abc_vel((0:5)+i_rows, (0:2)+i_cols) = [xyz2abc_vel40;xyz2abc_vel41];
        M_cs_abc2xyzT((0:5)+i_rows, (0:2)+i_cols) = [xyz2abc_vel40;xyz2abc_vel41];
        i_cols = i_cols+3;
        i_rows = i_rows+6;
    else
        M_cs_xyz2abc_vel(i_rows,i_cols) = 1;
        M_cs_abc2xyzT(i_rows,i_cols) = 1;
        i_cols = i_cols+1;
        i_rows = i_rows+1;
    end
end

M_cs_abc2xyz = M_cs_abc2xyzT';

end

