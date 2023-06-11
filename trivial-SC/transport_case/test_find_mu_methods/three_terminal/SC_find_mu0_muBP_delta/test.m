clearvars
tic
aa=rand(4000);
bb = aa\eye(4000) -inv(aa);
% bb = inv(aa);
% cc=bb* rand(4000);
toc