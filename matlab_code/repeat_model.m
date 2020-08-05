function [p1repeat, exp_Nreps]= repeat_model(N0,q)

for i = 1:N0
    probi(i) = i/q;
end

p1repeat = sum(probi);

for j =1:N0
    expnrepeat(j)=j.*((p1repeat).^j);

end

exp_Nreps = sum(expnrepeat);







end