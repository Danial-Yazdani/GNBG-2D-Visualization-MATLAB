%***********************************GNBG***********************************
%Author: Danial Yazdani
%Last Edited: August 26, 2023
%
% --------
% License:
% --------
% This program is to be used under the terms of the GNU General Public License
% (http://www.gnu.org/copyleft/gpl.html).
% Author: Danial Yazdani
% e-mail: danial DOT yazdani AT gmail DOT com
% Copyright notice: (c) 2023 Danial Yazdani
%**************************************************************************
function [result,GNBG] = DummyFitness(X,GNBG)
    x = X';
    f=NaN(1,GNBG.o);
    for k=1 : GNBG.o
        a = Transform((x - GNBG.MLO_BottomPosition(k,:)')'*GNBG.RotationMatrix(:,:,k)',GNBG.Alpha(k,:),GNBG.Beta(k,:));
        b = Transform(GNBG.RotationMatrix(:,:,k) * (x - GNBG.MLO_BottomPosition(k,:)'),GNBG.Alpha(k,:),GNBG.Beta(k,:));
        f(k) = GNBG.MLO_Sigma(k) + ( a * diag(GNBG.MLO_H(k,:)) * b)^GNBG.lambda(k);
    end
    result = min(f);
end

function Y = Transform(X,Alpha,Beta)
Y = X;
tmp = (X > 0);
Y(tmp) = log(X(tmp));
Y(tmp) = exp(Y(tmp) + Alpha(1)*(sin(Beta(1).*Y(tmp)) + sin(Beta(2).*Y(tmp))));
tmp = (X < 0);
Y(tmp) = log(-X(tmp));
Y(tmp) = -exp(Y(tmp) + Alpha(2)*(sin(Beta(3).*Y(tmp)) + sin(Beta(4).*Y(tmp))));
end