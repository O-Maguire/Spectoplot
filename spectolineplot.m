%Function writen by Oisin Maguire [July 2015], (last update March 2017, changes to example and description)
%Function takes 1, 3, 5, 6, 7 or 8 input arguments
%Function requires maximize: email oisin.maguire@ucdconnect.ie
%   if required: (It should still be on file exchange)
% [1]: (.spec.data/.spec)
% [3]: (.spec.data/.spec, cowantextdata column of config, config #)
% [5]: (.spec.data/.spec, cowantextdata column of config1, config1 #,
%                       cowantextdata column of config2, config2 #))
% [6]: (.spec.data/.spec, cowantextdata column of config1, config1 #,
%                       cowantextdata column of config2, config2 #,
%                       wavelength upper limit)
% [7] (.spec.data/.spec, cowantextdata column of config1, config1 #,
%                       cowantextdata column of config2, config2 #,
%                       wavelength upper limit, color)
% [8] (.spec.data/.spec, cowantextdata column of config1, config1 #,
%                       cowantextdata column of config2, config2 #,
%                       wavelength upper limit, color, wavelength
%                       lower limit)
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
% %spectolineplot(data_Ion);% line plot version
% %spectolineplot(data_Ion,textdata1,0);% cherry picked line plot 
% %spectolineplot(data_Ion,textdata1,0,textdata2,3);% Cherry cherry line plot
% %spectolineplot(data_Ion,textdata1,0,textdata2,4,200);% with upper wavelength limit in angstroms
% %spectolineplot(data_Ion,textdata1,0,textdata2,4,200,-1);% with colour
% %spectolineplot(data_Ion,textdata1,0,textdata2,4,200,0,10);% lower wavelength limit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% config1 / config2 ==0 means all of these configurations
% wavelength limit is in angstroms
%
% Colour:-0=red,1=green,2=yellow,3=magenta,4=cyan,5=black, else=blue
%
% Any bugs found please let me know, or if you fix them please send it on
% Cheers, happy coding
%
%acknowledgments to Dom, Regava, Emma and Paddy (in no particular order)
function output=spectolineplot(cowandata,cowantext1,J1,cowantext2,J2,uplam,col,downlam)
number_of_transitions=length(cowandata(:,1));
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
    elseif nargin==8%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Lower lamdba
        
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


if gA(1)==0
    disp('No transitions found.')
    return
else
    if nargin>=7%color
        gA(:);lamdba(:)./10;
        for i=1:length(gA(:))
            x=[lamdba(i)./10 lamdba(i)./10];
            y=[0 gA(i)];
            if col==0
                plot(x,y,'r')
            elseif col==1
                plot(x,y,'g')
            elseif col==2
                plot(x,y,'y')
            elseif col==3
                plot(x,y,'m')
            elseif col==4
                plot(x,y,'c')
            elseif col==5
                plot(x,y,'k')
            else
                plot(x,y,'b')
            end
            hold on
        end
        set(gca,'FontSize',26)
        ylabel('gA value')
        xlabel('Wavelength [nm]')
        maximize  
    else%No color i.e. blue
        gA(:);lamdba(:)./10;
        for i=1:length(gA(:))
            x=[lamdba(i)./10 lamdba(i)./10];
            y=[0 gA(i)];
            plot(x,y)
            hold on
        end
        set(gca,'FontSize',26)
        ylabel('gA value')
        xlabel('Wavelength [nm]')
        maximize
    end
end