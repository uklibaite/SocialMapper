function [y] = findAngleDiff(vec1, vec2)

n = length(vec1);
y = zeros(n,1);

for i = 1:n
    y(i) = acos(dot(vec1(i,:),vec2(i,:))./(norm(vec1(i,:)).*norm(vec2(i,:))));
end