%{
Title: Assignment 1
Course Code: AESM1511
Name: Karla Echeverry Caro
Studentnumber: 4537270
Date Created: 22 September 2022
Date modified: 26 September 2022
Mail: K.V.EcheverryCaro@student.tudelft.nl
Assignment 1: Linear Time-Invariant systems and the Principle of
Superposition
%}
%% Initial Settings 
close all;clc;clear;
%% Initial Inputs
%Original Signal
ini=(-pi/2); % Start point of the signal
fin=((3*pi)/2); % End point of the signal
steps=1257; % Sampling interval
x=linspace(ini,fin,steps); % Interval 
MasterSignalValues=sin(x); % Type of Signal
%%
%Question 2
len=(abs(fin)+abs(ini))/steps;  %Calculates the sampling length 
fprintf('The start point is %.3f.\n\nThe end point is %.3f.\n',ini,fin);
fprintf('\nThe sampling length is %.3f.\n',len);
fprintf('\nThe length of the signal is %d.\n',length(x));
%%
%Plot Original signal
figure();
plot(x,MasterSignalValues)
title('Original Signal')
ylabel('Amplitude')
xlabel('Angle')
legend('Original Signal','Location','NorthEastOutside')
set(gca,'XTick',ini:pi/2:fin) 
set(gca,'XTickLabel',{'-π/2','0','π/2','π','3π/2'})
%%
%Question 3
%Scaled Signals
scale=1:2:11; % Odd numbers
resh=reshape(scale,6,1); 
%{ 
Here we use the fact that the multiplication os matrices form a new matrix 
with the following dimensions [mxn]*[nxp]=[m*p] 
original signal = [1 x 1257]
Scale Vector = [1 x 6]
original signal  * Scale Vector = not possible  but reshaping the scale vector to
orginal signal = [6 x 1] will create a matrix with size [6 x 1257]
Every signal is in a new row.
Linear Algebra rocks :D
%}
newsignal=sin((x.*resh))./resh;
%%
%Question 4
sum_of_first_2=sum(newsignal(1:1:2, :),1); %The sum of the first and the second column/signals
%The same follows below, since for loops are not allowed :P
sum_of_first_3=sum(newsignal(1:1:3,:),1);
sum_of_first_4=sum(newsignal(1:1:4,:),1);
sum_of_first_5=sum(newsignal(1:1:5,:),1);
sum_of_first_6=sum(newsignal(1:1:6,:),1);
sum_of_six=sum(newsignal,1);
%% Plot the scaled signals

figure();
hold on 
% Orginal Signal
plot(x,newsignal(1,:))
% Cummulative sum
plot(x,sum_of_first_2)
plot(x,sum_of_first_3)
plot(x,sum_of_first_4)
plot(x,sum_of_first_5)
plot(x,sum_of_first_6)
title('Scaled Signals')
ylabel('Amplitude')
xlabel('Angle')
legend('Original Signal','Sum of scale 1, 3','Sum of scale 1, 3, 5',... 
    'Sum of scale 1, 3, 5, 7','Sum of scale 1, 3, 5, 7, 9',...
    'Sum of scale 1, 3, 5, 7, 9, 11','Location','NorthEastOutside')
set(gca,'XTick',ini:pi/2:fin) 
set(gca,'XTickLabel',{'-π/2','0','π/2','π','3π/2'})

hold off
%% Question 4 Theoretical Questions
pic1=figure();
hold on 
% Orginal Signal
subplot(2,2,1);
plot(x,newsignal(1,:))
title('Original Signal')
ylabel('Amplitude')
xlabel('Angle')
set(gca,'XTick',ini:pi/2:fin) 
set(gca,'XTickLabel',{'-π/2','0','π/2','π','3π/2'})

% Scaled Signal *11
subplot(2,2,2);
plot(x,newsignal(6,:))
title('Scaled Signal with scaling factor 11')
ylabel('Amplitude')
xlabel('Angle')
set(gca,'XTick',ini:pi/2:fin) 
set(gca,'XTickLabel',{'-π/2','0','π/2','π','3π/2'})

%Summation of Scaled Signal*11 + Original Signal
subplot(2,2,3);
plot(x,sum(newsignal(1:5:6,:)))
title('Scaled Signal with scaling factor 11 + Original Signal')
ylabel('Amplitude')
xlabel('Angle')
set(gca,'XTick',ini:pi/2:fin) 
set(gca,'XTickLabel',{'-π/2','0','π/2','π','3π/2'})

hold off

fprintf(['\nWhat happens to the original master signal? \n' ...
    'The original Signal is being scaled, which makes the angular frequency.\n'...
    'As can be seen in figure %g.\n'],pic1);

fprintf(['\nHow does it change after adding additional signals?\n' ...
    'It becomes the summation of the both signal.\n'...
    'As can be seen in figure %g.\n'],pic1);

fprintf(['\nWhat do you think the summed signal would look like if you '...
           'continue adding more and\nmore signals reaching, for example,' ...
           '1000 summed signals?\n'...
           'Adding more signals will result in a periodic signal in the ',...
           'from of a unit step function,\nlooking similarly to an even ',...
           'function of the following form:-u(x+π)+u(x)-u(x-π) repeating\n',...
           'it self constantly with a fundamental period of 2π. In a dirac ', ...
           'pulse train function it would\nlook like this Σ(δ(x-mπ+π)-δ(x-mπ))*u(x)' ...
           ' ∀ m=(-∞,∞)\n']);
%%
% Question 5  
sum_of_odd_rows=sum(newsignal(1:2:end,:),1); % Summation of odd rows in the matrix
sum_of_even_rows=sum(newsignal(2:2:end,:),1);% Summation of even rows in the matrix
sum_of_odd_even=sum_of_odd_rows+sum_of_even_rows;% Summation of odd and even rows in the matrix
%% 
%Plot Rows

figure();
hold on 
plot(x,sum_of_odd_rows)
plot(x,sum_of_even_rows)
plot(x,sum_of_odd_even)

title('Rows Summation')
ylabel('Amplitude')
xlabel('Angle')
legend('Sum of Odd rows','Sum of Even rows','Sum of all rows','Location','NorthEastOutside')
set(gca,'XTick',ini:pi/2:fin) 
set(gca,'XTickLabel',{'-π/2','0','π/2','π','3π/2'})

hold off
%%
% Question 6

figure();
hold on 
plot(x,sum_of_odd_even)
plot(x,sum_of_first_6)

title('Comparissons')
ylabel('Amplitude')
xlabel('Angle')
legend('Sum of odd and even signals','Sum of all the signals','Location','NorthEastOutside')
set(gca,'XTick',ini:pi/2:fin) 
set(gca,'XTickLabel',{'-π/2','0','π/2','π','3π/2'})

hold off

fprintf(['\nDoes the order of summation matter in obtaining the final result?\n',...
           'For Linear functions the order in which you add, multiply or subtract '...
           'does not change the \n final result, so it must satisfy the additivity '...
           'and homogeneity meaning f(x + y) = f(x) + f(y)\n and f(αx) = α f(x) ∀ α.\n'])


fprintf(['\nDoes the order of summation matter in obtaining intermediate results?\n',...
           'The order of summation has a result in the intermediate results but not'...
           'the final results.\n'])

% ------------- END OF CODE --------------
