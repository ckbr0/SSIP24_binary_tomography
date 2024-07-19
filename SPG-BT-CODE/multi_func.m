function value=multi_func(vec)


[n m]=size(vec);
value=0;
for i=1:n,
    value=value+vec(i);
    
end;

