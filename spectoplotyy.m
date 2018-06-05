%Function writen by Oisin Maguire [June 2015], (last update March 2017, changes to example and description)
%Function takes 10 input arguments
%Function requires maximize: email oisin.maguire@ucdconnect.ie
%   if required: (It should still be on file exchange)
% [10] (.spec.data/.spec, cowantextdata column of config1, config1 #,
%                       cowantextdata column of config2, config2 #,
%                       wavelength upper limit, color, wavelength
%                       lower limit, experimental wavelength, spectrum)
%and returns a convolved cowan
%spectrum with a 0.02 std Gaussian, as well as an array with
%wavelength in angstroms in the first column
%gA values in the second
%gf values in the third
%Energy in the fourth
%
%                           Sample Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%close all; clear all; clc;
%Ion=importdata('Ge8_90_99.spec');
%data_Ion=Ion.data; text_Ion=Ion.textdata;
%config=length(text_Ion)-length(data_Ion);%Number of configurations within the .spec file
%for i=1:length(data_Ion(:,1))
%     I=i+config; 
%     textdata1(i)=str2num(cell2mat(text_Ion(I,4)));
%     textdata2(i)=str2num(cell2mat(text_Ion(I,9)));
%end
%experimentalData=importdata('Experimental.txt');
%X=experimentalData(:,1);%Wavelength in nm.
%Y=experimentalData(:,2);%Intensity in arbitrary units
% %spectoplotyy(data_Ion,textdata1,0,textdata2,0,200,0,10,X,Y);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% config1 / config2 ==0 means all of these configurations
% wavelength limit is in angstroms, its to help with convolution, time and
% accuracy for the normal 1:50 nm region
%
% Any bugs found please let me know, or if you fix them please send it on
% Cheers, happy coding
%
%acknowledgments to Dom, Regava, Emma and Paddy (in no particular order)
function output=spectoplotyy(cowandata,cowantext1,J1,cowantext2,J2,uplam,col,downlam,X1,Y1)
number_of_transitions=length(cowandata(:,1));
FWHM=0.02;
length_of_file=length(cowandata(1,:)); I=0; gA(1)=0;
if length_of_file==7 %ALL nargins inside............
    if nargin==1%%%%%%ALL configuration...............Working %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i=1:number_of_transitions%All
            if I==0
                disp('All transitions')
            end
            I=I+1;
            gA(I)=cowandata(i,6);
            lamdba(I)=cowandata(i,2);
            output(I,1)=cowandata(i,2);%wavelength
            output(I,2)=cowandata(i,6);%gA
            output(I,3)=cowandata(i,4);%gf
            output(I,4)=cowandata(i,1);%Energy
        end
    elseif nargin==3%%%%%%%%%%%%%%%single configurations..........Working %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if J1==0
            for i=1:number_of_transitions%All
                if I==0
                    disp('All transitions')
                end
                I=I+1;
                gA(I)=cowandata(i,6);
                lamdba(I)=cowandata(i,2);
                output(I,1)=cowandata(i,2);%wavelength
                output(I,2)=cowandata(i,6);%gA
                output(I,3)=cowandata(i,4);%gf
                output(I,4)=cowandata(i,1);%Energy
            end
        else
            for i=1:number_of_transitions%Cherry
                if I==0
                    disp('Cherry transitions')
                end
                if cowantext1(i)==J1
                    I=I+1;
                    gA(I)=cowandata(i,6);
                    lamdba(I)=cowandata(i,2);
                    output(I,1)=cowandata(i,2);%wavelength
                    output(I,2)=cowandata(i,6);%gA
                    output(I,3)=cowandata(i,4);%gf
                    output(I,4)=cowandata(i,1);%Energy
                end
            end
        end
    elseif nargin==5 %%%%%%%%%%%%%%%%%%%% mixed configurations working...:-) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if J1==0
            if J2==0
                for i=1:number_of_transitions%all all%Fine
                    if I==0
                        disp('All transitions')
                    end
                    I=I+1;
                    gA(I)=cowandata(i,6);
                    lamdba(I)=cowandata(i,2);
                    output(I,1)=cowandata(i,2);%wavelength
                    output(I,2)=cowandata(i,6);%gA
                    output(I,3)=cowandata(i,4);%gf
                    output(I,4)=cowandata(i,1);%Energy
                end
            else
                for  i=1:number_of_transitions%All cherry%Fine
                    if cowantext2(i)==J2
                        if I==0
                            disp('All Cherry')
                        end
                        I=I+1;
                        gA(I)=cowandata(i,6);
                        lamdba(I)=cowandata(i,2);
                        output(I,1)=cowandata(i,2);%wavelength
                        output(I,2)=cowandata(i,6);%gA
                        output(I,3)=cowandata(i,4);%gf
                        output(I,4)=cowandata(i,1);%Energy
                    end
                end
            end
        end
        if J2==0
            for i=1:number_of_transitions%Cherry all%Fine
                if cowantext1(i)==J1
                    if I==0
                        disp('Cherry all')
                    end
                    I=I+1;
                    gA(I)=cowandata(i,6);
                    lamdba(I)=cowandata(i,2);
                    output(I,1)=cowandata(i,2);%wavelength
                    output(I,2)=cowandata(i,6);%gA
                    output(I,3)=cowandata(i,4);%gf
                    output(I,4)=cowandata(i,1);%Energy
                end
            end
        end
        if J1>=1
            if J2>=1
                for i=1:number_of_transitions%Cherry Cherry%Fine
                    if cowantext1(i)==J1
                        if cowantext2(i)==J2
                            if I==0
                                disp('Cherry cherry')
                            end
                            I=I+1;
                            gA(I)=cowandata(i,6);
                            lamdba(I)=cowandata(i,2);
                            output(I,1)=cowandata(i,2);%wavelength
                            output(I,2)=cowandata(i,6);%gA
                            output(I,3)=cowandata(i,4);%gf
                            output(I,4)=cowandata(i,1);%Energy
                        end
                    end
                end
            end
        end
    elseif nargin==6%%%%%%%%%%%%%%%%%%%%%% mixed configurations working, needs testing, but seems to work fine%%%%%
        if J1==0
            if J2==0
                for i=1:number_of_transitions%all all%Fine
                    if cowandata(i,2)<=uplam
                        if I==0
                            disp('All transitions')
                        end
                        I=I+1;
                        gA(I)=cowandata(i,6);
                        lamdba(I)=cowandata(i,2);
                        output(I,1)=cowandata(i,2);%wavelength
                        output(I,2)=cowandata(i,6);%gA
                        output(I,3)=cowandata(i,4);%gf
                        output(I,4)=cowandata(i,1);%Energy
                    end
                end
            else
                for  i=1:number_of_transitions%All cherry%Fine
                    if cowantext2(i)==J2
                        if cowandata(i,2)<=uplam
                            if I==0
                                disp('All Cherry')
                            end
                            I=I+1;
                            gA(I)=cowandata(i,6);
                            lamdba(I)=cowandata(i,2);
                            output(I,1)=cowandata(i,2);%wavelength
                            output(I,2)=cowandata(i,6);%gA
                            output(I,3)=cowandata(i,4);%gf
                            output(I,4)=cowandata(i,1);%Energy
                        end
                    end
                end
            end
        end
        if J2==0
            for i=1:number_of_transitions%Cherry all%Fine
                if cowantext1(i)==J1
                    if cowandata(i,2)<=uplam
                        if I==0
                            disp('Cherry all')
                        end
                        I=I+1;
                        gA(I)=cowandata(i,6);
                        lamdba(I)=cowandata(i,2);
                        output(I,1)=cowandata(i,2);%wavelength
                        output(I,2)=cowandata(i,6);%gA
                        output(I,3)=cowandata(i,4);%gf
                        output(I,4)=cowandata(i,1);%Energy
                    end
                end
            end
        end
        if J1>=1
            if J2>=1
                for i=1:number_of_transitions%Cherry Cherry%Fine
                    if cowantext1(i)==J1
                        if cowantext2(i)==J2
                            if cowandata(i,2)<=uplam
                                if I==0
                                    disp('Cherry cherry')
                                end
                                I=I+1;
                                gA(I)=cowandata(i,6);
                                lamdba(I)=cowandata(i,2);
                                output(I,1)=cowandata(i,2);%wavelength
                                output(I,2)=cowandata(i,6);%gA
                                output(I,3)=cowandata(i,4);%gf
                                output(I,4)=cowandata(i,1);%Energy
                            end
                        end
                    end
                end
            end
        end
        
    elseif nargin==7 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if J1==0
            if J2==0
                for i=1:number_of_transitions%all all%Fine
                    if cowandata(i,2)<=uplam
                        if I==0
                            disp('All transitions')
                        end
                        I=I+1;
                        gA(I)=cowandata(i,6);
                        lamdba(I)=cowandata(i,2);
                        output(I,1)=cowandata(i,2);%wavelength
                        output(I,2)=cowandata(i,6);%gA
                        output(I,3)=cowandata(i,4);%gf
                        output(I,4)=cowandata(i,1);%Energy
                    end
                end
            else
                for  i=1:number_of_transitions%All cherry%Fine
                    if cowantext2(i)==J2
                        if cowandata(i,2)<=uplam
                            if I==0
                                disp('All Cherry')
                            end
                            I=I+1;
                            gA(I)=cowandata(i,6);
                            lamdba(I)=cowandata(i,2);
                            output(I,1)=cowandata(i,2);%wavelength
                            output(I,2)=cowandata(i,6);%gA
                            output(I,3)=cowandata(i,4);%gf
                            output(I,4)=cowandata(i,1);%Energy
                        end
                    end
                end
            end
        end
        if J2==0
            for i=1:number_of_transitions%Cherry all%Fine
                if cowantext1(i)==J1
                    if cowandata(i,2)<=uplam
                        if I==0
                            disp('Cherry all')
                        end
                        I=I+1;
                        gA(I)=cowandata(i,6);
                        lamdba(I)=cowandata(i,2);
                        output(I,1)=cowandata(i,2);%wavelength
                        output(I,2)=cowandata(i,6);%gA
                        output(I,3)=cowandata(i,4);%gf
                        output(I,4)=cowandata(i,1);%Energy
                    end
                end
            end
        end
        if J1>=1
            if J2>=1
                for i=1:number_of_transitions%Cherry Cherry%Fine
                    if cowantext1(i)==J1
                        if cowantext2(i)==J2
                            if cowandata(i,2)<=uplam
                                if I==0
                                    disp('Cherry cherry')
                                end
                                I=I+1;
                                gA(I)=cowandata(i,6);
                                lamdba(I)=cowandata(i,2);
                                output(I,1)=cowandata(i,2);%wavelength
                                output(I,2)=cowandata(i,6);%gA
                                output(I,3)=cowandata(i,4);%gf
                                output(I,4)=cowandata(i,1);%Energy
                            end
                        end
                    end
                end
            end
        end
    elseif nargin>=8%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Lower lamdba
        
        if J1==0
            if J2==0
                for i=1:number_of_transitions%all all%Fine
                    if cowandata(i,2)<=uplam
                        if cowandata(i,2)>=downlam
                            if I==0
                                disp('All transitions')
                            end
                            I=I+1;
                            gA(I)=cowandata(i,6);
                            lamdba(I)=cowandata(i,2);
                            output(I,1)=cowandata(i,2);%wavelength
                            output(I,2)=cowandata(i,6);%gA
                            output(I,3)=cowandata(i,4);%gf
                            output(I,4)=cowandata(i,1);%Energy
                        end
                    end
                end
            else
                for  i=1:number_of_transitions%All cherry%Fine
                    if cowantext2(i)==J2
                        if cowandata(i,2)<=uplam
                            if cowandata(i,2)>=downlam
                                if I==0
                                    disp('All Cherry')
                                end
                                I=I+1;
                                gA(I)=cowandata(i,6);
                                lamdba(I)=cowandata(i,2);
                                output(I,1)=cowandata(i,2);%wavelength
                                output(I,2)=cowandata(i,6);%gA
                                output(I,3)=cowandata(i,4);%gf
                                output(I,4)=cowandata(i,1);%Energy
                            end
                        end
                    end
                end
            end
        end
        if J2==0
            for i=1:number_of_transitions%Cherry all%Fine
                if cowantext1(i)==J1
                    if cowandata(i,2)<=uplam
                        if cowandata(i,2)>=downlam
                            if I==0
                                disp('Cherry all')
                            end
                            I=I+1;
                            gA(I)=cowandata(i,6);
                            lamdba(I)=cowandata(i,2);
                            output(I,1)=cowandata(i,2);%wavelength
                            output(I,2)=cowandata(i,6);%gA
                            output(I,3)=cowandata(i,4);%gf
                            output(I,4)=cowandata(i,1);%Energy
                        end
                    end
                end
            end
        end
        if J1>=1
            if J2>=1
                for i=1:number_of_transitions%Cherry Cherry%Fine
                    if cowantext1(i)==J1
                        if cowantext2(i)==J2
                            if cowandata(i,2)<=uplam
                                if cowandata(i,2)>=downlam
                                    if I==0
                                        disp('Cherry cherry')
                                    end
                                    I=I+1;
                                    gA(I)=cowandata(i,6);
                                    lamdba(I)=cowandata(i,2);
                                    output(I,1)=cowandata(i,2);%wavelength
                                    output(I,2)=cowandata(i,6);%gA
                                    output(I,3)=cowandata(i,4);%gf
                                    output(I,4)=cowandata(i,1);%Energy
                                end
                            end
                        end
                    end
                end
            end
        end
        
    else
        disp('I have no idea what you are doing.')
        return
    end
else
    disp('Not sure this is a spec file.');
end


if nargin==7
    if gA(1)==0
        disp('No transitions found.')
        return
    else
        %%Convolution...
        numLines=length(gA(:));
        numPoints=10000; stddev=FWHM; flgArea=1; minX=min(lamdba(:)./10)*.9;
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
        if col==0
            plot(X,GaussExpSum,'r')
        elseif col==1
            plot(X,GaussExpSum,'g')
        elseif col==2
            plot(X,GaussExpSum,'y')
        elseif col==3
            plot(X,GaussExpSum,'m')
        elseif col==4
            plot(X,GaussExpSum,'c')
        elseif col==5
            plot(X,GaussExpSum,'k')
        else
            plot(X,GaussExpSum,'b')
        end
        set(gca,'FontSize',26)
        ylabel('Convolved gA value')
        xlabel('Wavelength [nm]')
        maximize
    end
elseif nargin==8
    if gA(1)==0
        disp('No transitions found.')
        return
    else
        %%Convolution...
        numLines=length(gA(:));
        numPoints=10000; stddev=FWHM; flgArea=1; minX=min(lamdba(:)./10)*.9;
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
        if col==0
            plot(X,GaussExpSum,'r')
        elseif col==1
            plot(X,GaussExpSum,'g')
        elseif col==2
            plot(X,GaussExpSum,'y')
        elseif col==3
            plot(X,GaussExpSum,'m')
        elseif col==4
            plot(X,GaussExpSum,'c')
        elseif col==5
            plot(X,GaussExpSum,'k')
        else
            plot(X,GaussExpSum,'b')
        end
        set(gca,'FontSize',26)
        ylabel('Convolved gA value')
        xlabel('Wavelength [nm]')
        maximize
    end
elseif nargin>=9
    
    if gA(1)==0
        disp('No transitions found.')
        return
    else
        %%Convolution...
        numLines=length(gA(:));
        numPoints=10000; stddev=FWHM; flgArea=1; minX=min(lamdba(:)./10)*.9;
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
        if col==0
            [ax h1 h2]=plotyy(X,GaussExpSum,X1,Y1);
        elseif col==1
            [ax h1 h2]=plotyy(X,GaussExpSum,X1,Y1)
        elseif col==2
            [ax h1 h2]=plotyy(X,GaussExpSum,X1,Y1)
        elseif col==3
            [ax h1 h2]=plotyy(X,GaussExpSum,X1,Y1)
        elseif col==4
            [ax h1 h2]=plotyy(X,GaussExpSum,X1,Y1)
        elseif col==5
            [ax h1 h2]=plotyy(X,GaussExpSum,X1,Y1)
        else
            [ax h1 h2]=plotyy(X,GaussExpSum,X1,Y1)
        end
        set(gca,'FontSize',20)
        %ylabel('Convolved gA value') % left y-axis
        %ylabel('Experimental Spectrum')
        %set(ax(2),'xtick',get(ax(1),'xtick'),'xticklab',get(ax(1),'xticklab'))
        linkaxes(ax,'x');
        set(ax(2),'XTickLabel',[])
        axes(ax(1)); ylabel('Convolved gA value','FontSize',20);
        axes(ax(2)); ylabel('Experimental Spectrum','FontSize',20);
        xlabel('Wavelength [nm]','FontSize',20)
        maximize
    end
    
else
    %Sanity check for transitions
    if gA(1)==0
        disp('No transitions found.')
        return
    else
        %%Convolution...
        numLines=length(gA(:));
        numPoints=10000; stddev=FWHM; flgArea=1; minX=min(lamdba(:)./10)*.9;
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
end
end