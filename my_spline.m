function [ yy ] = my_spline( X,Y,xx )
%MY_SPLINE does spline of Y to interolate @ xx. X is
%original data position
%   Detailed explanation goes here
% X : Ky
% xx : Ky_int
% Y : sarDataFFT

xx_row = size(xx,1);

if xx_row == 1, %xx is transposed vector
    xx=xx.';X=X.';Y=Y.';
end

Y_col=size(Y,2); %number of columns
yy=zeros(length(xx),Y_col);

for i=1:Y_col
    yy(:,i)=spline( X(:,i),Y(:,i),xx );

    yy(min(X(:,i)) > xx, i) = 1e-30;
    yy(max(X(:,i)) < xx, i) = 1e-30;
end

if xx_row == 1, %xx is transposed vector
    yy=yy.';
end

end

