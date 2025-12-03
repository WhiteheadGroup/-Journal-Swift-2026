%Inverter Model
%SDS

clear
clc

%Set Initial Conditions
kd1 = 115.6; %nM
kd2 = 2.1; %nM
kdoff = 200; %nM
p = 200; %nM
x = 20; %nM
h = 2; %nM

V = 1;

%Create Variables
syms c s pl

%Set Ligand Range
ligand = round(logspace(-2,4,10),1,'decimals');

%High
Solutions  = zeros(length(ligand),5);
%Create For Loop
for L = 1:length(ligand)
    %Solve System of EQs
    S = solve(c^2 == c*(kdoff+x+h-s)+x*(s-h), s == pl*(h-c)/(kd2+pl),pl^2 == pl*(p+ligand(L)-2*s+kd1)+(s-p)*(ligand(L)-s));    

    %Store Solutions
    Solutions(L,:) = [ligand(L),vpa(S.c(V, 1)),vpa(S.pl(V, 1)),vpa(S.s(V, 1)),vpa(S.c(V, 1)) + vpa(S.s(V, 1))];
end

%Graph
figure
x = Solutions(:,1);
y = Solutions(:,2);
semilogx(x,y)

title('Inverter')
xlabel('Ligand [nM]') 
ylabel('Complex [nM]')
