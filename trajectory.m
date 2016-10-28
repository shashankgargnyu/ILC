t=(0:0.01:10);
l=length(t);
for i=0:l-1
h(i+1)=t(i+1)^3*(4-0.3*t(i+1))*0.01;
end
plot(t,h);