a = ecoli

[D,L] = oc_set(a,1);

d = D.data(L==1,:);

size(d)

output = [D.data,L];

size(output)

csvwrite('ecoli.csv',output)