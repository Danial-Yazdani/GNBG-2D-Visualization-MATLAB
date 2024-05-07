%***********************************GNBG***********************************
%Author: Danial Yazdani
%Last Edited: August 26, 2023
%Title: Generalized Numerical Benchmark Generator ==> Visualization in 2-D
% --------
% License:
% This program is to be used under the terms of the GNU General Public License
% (http://www.gnu.org/copyleft/gpl.html).
% Author: Danial Yazdani
% e-mail: danial DOT yazdani AT gmail DOT com
% Copyright notice: (c) 2023 Danial Yazdani
%************************************************************************** 
close all;clear all;clc;
GNBG                     = [];
GNBG.o                   = 5; % Number of components
GNBG.Dimension           = 2; % Must be set to 2 for visualization
GNBG.MinCoordinate       = -100;
GNBG.MaxCoordinate       = 100;
GNBG.MinRandOptimaPos    = -80;
GNBG.MaxRandOptimaPos    = 80;
GNBG.MinSigma            = 0;
GNBG.MaxSigma            = 20;
GNBG.MinH                = 1;
GNBG.MaxH                = 5;
GNBG.MinAngle            = -pi;
GNBG.MaxAngle            = pi;
GNBG.MinAlpha            = 0;
GNBG.MaxAlpha            = 0.5;
GNBG.MinBeta             = 5; 
GNBG.MaxBeta             = 50;
GNBG.MaxLambda           = 0.5;
GNBG.MinLambda           = 0.5;
GNBG.MLO_Sigma           = GNBG.MinSigma + (GNBG.MaxSigma-GNBG.MinSigma)*rand(GNBG.o,1);
GNBG.MLO_BottomPosition  = (GNBG.MinCoordinate + (GNBG.MaxCoordinate-GNBG.MinCoordinate)*rand(GNBG.o,GNBG.Dimension));
% GNBG.MLO_BottomPosition  = (GNBG.MinRandOptimaPos + (GNBG.MaxRandOptimaPos-GNBG.MinRandOptimaPos)*rand(GNBG.o,GNBG.Dimension));
% GNBG.MLO_H               = GNBG.MinH + (GNBG.MaxH-GNBG.MinH)*betarnd(0.2, 0.2, [GNBG.o,GNBG.Dimension]);
GNBG.MLO_H               = GNBG.MinH + (GNBG.MaxH-GNBG.MinH)*repmat(rand(GNBG.o,1),1,GNBG.Dimension);
GNBG.Alpha               = GNBG.MinAlpha + (GNBG.MaxAlpha-GNBG.MinAlpha)*rand(GNBG.o,2);
GNBG.Beta                = GNBG.MinBeta + (GNBG.MaxBeta-GNBG.MinBeta)*rand(GNBG.o,4);
GNBG.lambda              = GNBG.MinLambda + (GNBG.MaxLambda-GNBG.MinLambda)*rand(GNBG.o,1);
GNBG.Angle               = GNBG.MinAngle + (GNBG.MaxAngle-GNBG.MinAngle)*rand(GNBG.o,1);
GNBG.Rotation            = 2;%(0) Without rotation
                             %(1) Random Rotation for all main optima
                             %(2) Rotation based on the specified Angle for each main local optimum
switch GNBG.Rotation
    case 0
        GNBG.RotationMatrix=NaN(GNBG.Dimension,GNBG.Dimension,GNBG.o);
        for ii=1 : GNBG.o
            GNBG.RotationMatrix(:,:,ii) = eye(GNBG.Dimension);
        end
    case 1
        GNBG.RotationMatrix=NaN(GNBG.Dimension,GNBG.Dimension,GNBG.o);
        for ii=1 : GNBG.o
            GNBG.P=eye(GNBG.Dimension);
            upperTriangle = triu(true(GNBG.Dimension), 1); % Logical matrix for upper triangle
            GNBG.P(upperTriangle) = GNBG.MinAngle + (GNBG.MaxAngle - GNBG.MinAngle) * rand(sum(upperTriangle(:)), 1);
            [GNBG.RotationMatrix(:,:,ii)] = Rotation(GNBG.P,GNBG.Dimension);
        end
    case 2
        GNBG.RotationMatrix=NaN(GNBG.Dimension,GNBG.Dimension,GNBG.o);
        for ii=1 : GNBG.o
            GNBG.P=eye(GNBG.Dimension);
            upperTriangle = triu(true(GNBG.Dimension), 1); % Logical matrix for upper triangle
            GNBG.P(upperTriangle) = GNBG.Angle(ii);
            [GNBG.RotationMatrix(:,:,ii)] = Rotation(GNBG.P,GNBG.Dimension);
        end
    otherwise
        warning('Wrong number is chosen for GNBG.Rotation.')
end

[GNBG.OptimumValue,GlobalOptimumID] = min(GNBG.MLO_Sigma);
GNBG.OptimumPosition     = GNBG.MLO_BottomPosition(GlobalOptimumID,:);
%% Visualization data generation
T = GNBG.MinCoordinate : ( GNBG.MaxCoordinate-GNBG.MinCoordinate)/500 :  GNBG.MaxCoordinate;
L=length(T);
F=zeros(L);
for i=1:L
    for j=1:L
        F(i,j) = DummyFitness([T(i), T(j)],GNBG);
    end
end
%% Surface+contour plot
figure;
surf(T,T,F,'FaceAlpha',0.8, 'EdgeColor', 'none');
hold on;
contour(T,T,F,25,'ZLocation','zmin');
hold off;
x1lh = ylabel('x_1');
x2lh = xlabel('x_2');
zlabel('f(x_1,x_2)')
colormap jet
grid on
set(gcf,'OuterPosition',[150 150 600 550]);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
x2lh.Position(2) = x2lh.Position(2);



%% Rotation matrix generator function
function R = Rotation(teta,Dimension)
R = eye(Dimension);
for p=1 : (Dimension-1)
    for q=(p+1) : (Dimension)
        if teta(p,q)~=0
            G = eye(Dimension);
            G(p,p) = cos(teta(p,q));
            G(q,q) = cos(teta(p,q));
            G(p,q) = -sin(teta(p,q));
            G(q,p) = sin(teta(p,q));
            R = R*G;
        end
    end
end
end