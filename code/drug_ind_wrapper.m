function [tcrit1] = drug_ind_wrapper(ptest,Cdoxvec)
Cdoxmax = 1000;
tvec1=[0:3:450];
tscRNAseq = 1656;
tgen1 =[0:1:tvec1(end)];
tvec2 =[0:3:tscRNAseq];
tgen2 = [0:1:tscRNAseq];
tdrug = 1;
k = 0.5;
kdrug = 0.0175;

for i = 1:length(Cdoxvec)
U1(:,i)=k*Cdoxvec(i)*exp(-kdrug*(tgen1))/(0.1*Cdoxmax);
dt = 1;
tdrug = 1;
N0 = 2e3;
p1 = [ptest(1:2), ptest(4:8)];
% Use this to find the critical times as a function of dose in 96 well
[Nsrdat1(:,:,i), tcrit1(i), Ncrit1(i)]=fwd_Greene_model2(p1, tvec1, N0, U1(:,i), dt,tdrug);
% Use this to find the phi(@tscRNAseq) as a function of dose in expansion
tgen2 = [0:1:tscRNAseq];
U2(:,i)=k*Cdoxvec(i)*exp(-kdrug*(tgen2))/(0.1*Cdoxmax);
p2 =p1;
% replace carcaph
p2(2) = ptest(3);
[Nsrdat2(:,:,i), tcrit2(i), Ncrit2(i)]=fwd_Greene_model2(p2, tvec2, N0, U2(:,i), dt,tdrug);
% Use this to find phi@tSCRNAseq
% This could be found at any time, we just pick last bc last is the time of
% scNAseq
phi_sc(i) = Nsrdat2(end,2,i)/Nsrdat2(end,1,i);
end
end

