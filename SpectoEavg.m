close all; clear all; clc;
%% configuables
Ion=importdata('DyXX_5s15p1.spec');
lower_state='3';
J_lower='4.5';
upper_state='9';
J_upper='5.5';

%% constants
text_Ion=Ion.textdata; data_Ion=Ion.data;
config=length(text_Ion)-length(data_Ion);

%% removing header in textdata
text_Ion(1:config,:) = [];

%% selecting lower config
idx = ismember(text_Ion(:,4),lower_state); 
text_Ion= text_Ion(idx,:);

%% selecting upper config
idx = ismember(text_Ion(:,9),upper_state); 
text_Ion= text_Ion(idx,:);

%% selecting lower J
idx = ismember(text_Ion(:,3),J_lower); 
text_Ion= text_Ion(idx,:);

%% selecting upper J
idx = ismember(text_Ion(:,8),J_upper); 
text_Ion= text_Ion(idx,:);

%% should have clean data now...
%% creating delta E line
DeltaE=str2double(text_Ion(:,7))-str2double(text_Ion(:,2));
Eavg=mean(DeltaE)



