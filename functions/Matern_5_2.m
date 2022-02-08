function y = Matern_5_2(d,beta)
    x = sqrt(5)*beta.*d;
    y = (1+x+x.^2/3).*exp(-x);
end