%Function writen by Oisin Maguire [Nov 2020]
%
%                           Sample Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all; clear all; clc;
% Ion=importdata('CeXXII.spec');
% data_Ion=Ion.data; text_Ion=Ion.textdata;
% %% convolve all configurations
% spectoplot_pro(data_Ion,text_Ion);
% %% convolve just the lower and upper selected transitions
% lower_selection=[1 3];% lower config
% upper_selection=[1 2];% upper config
% selected_transitions=spectoplot_pro(data_Ion,text_Ion,lower_selection,upper_selection);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output=spectoplot_pro(cowandata,cowantext,lower,upper,stddev)
try stddev; catch stddev=0.02; end
if nargin==4
    config=length(cowantext)-length(cowandata);
    
    newarray=cell2num(cowantext(config+1:end,[4,9]));
    newarray= [newarray cowandata];
    
    idx = ismember(newarray(:,1),lower);
    newarray= newarray(idx,:);
    
    idx = ismember(newarray(:,2),upper);
    newarray= newarray(idx,:);
    
        gA=newarray(:,6+2);
    lamdba=newarray(:,2+2);
    output=newarray;
    
elseif nargin==2
    gA=cowandata(:,6);
    lamdba=cowandata(:,2);
end

%%Convolution...
numLines=length(gA(:));
numPoints=10000; flgArea=1; minX=min(lamdba(:)./10)*.9;
maxX=max(lamdba(:)./10)*1.1; b(:,1)=lamdba(:)./10; b(:,2)=gA(:);
incX = (maxX-minX) / numPoints; GaussAreaZero = 1 / sqrt(2*pi)/ stddev;
X = [minX:incX:maxX].'; GaussExpSum = zeros(size(X,1),1);
for i = 1:numLines
    GaussExp = exp (- ( (b(i,1) - X).^2 ) /(2*stddev^2) );
    GaussExpSum = GaussExpSum + b(i,2) * GaussExp;
end
if flgArea == 1
    GaussExpSum = GaussAreaZero * GaussExpSum;
end
plot(X,GaussExpSum)
set(gca,'FontSize',26)
ylabel('Convolved gA value')
xlabel('Wavelength [nm]')
maximize
end
