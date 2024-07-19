

tic; start_time=cputime;
for i=1:1000000,
   
end;
e_for=toc; e_for_cpu=cputime-start_time;



tic; start_time=cputime;
while i<1000001,
   
   i=i+1;
   
end;
e_while=toc; e_while_cpu=cputime-start_time;