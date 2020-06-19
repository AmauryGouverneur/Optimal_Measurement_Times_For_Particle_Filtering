function res = clip(x,a,b)

res = x.*(x>a & x<b) + a.*(x<=a) + b.*(x>=b);
end