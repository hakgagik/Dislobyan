x=centers(:,1);
y=centers(:,2);
z=centers(:,3);

testPts=(0:0.01:3);
response=zeros(length(testPts),1);
for k=1:length(response);

test_vector=c_star*testPts(k);

response(k)=abs(sum(exp(-2i*pi*(test_vector(1)*x+test_vector(2)*y+test_vector(3)*z))));
end

plot(testPts,response)