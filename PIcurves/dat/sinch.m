function y = sinch(x)

if x == 0
   y = 1;
   else
   y = (1 - exp(-x))./x;
end
