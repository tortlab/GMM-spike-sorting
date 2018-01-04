

%%

% 
% xpoly = vx;
% ypoly = vy;
% xline = xvalues;
% yline = yvalues;

function isIn = linexpoly(xpoly,ypoly,xline,yline)

xmin= min(xpoly);
xmax= max(xpoly);
%%
valid = [find(xline(1,:)<xmin,1,'last') find(xline(1,:)>xmax,1)];
xline = xline(1,valid(1):valid(2));
yline = yline(:,valid(1):valid(2));
%%
isIn=zeros(1,size(yline,1));
for ispike=1:size(yline,1)
      isIn(ispike)=intersect_poly(xpoly,ypoly,xline,yline(ispike,:));
end

end
%%

function [P] = intersect_line(x1,y1,x2,y2,x3,y3,x4,y4)

P=[];

% Ax = b
A = [(y2-y1),-(x2-x1);
    (y4-y3), -(x4-x3)];

b = [x1*(y2-y1)-y1*(x2-x1);
    x3*(y4-y3)-y3*(x4-x3)];

if det(A)~=0
    P = inv(A)*b;
end
end


function cross = intersect_poly(xpoly,ypoly,xline,yline)
cross = 0;
for ipoly = 1:length(xpoly)-1
    x1 = xpoly(ipoly);
    x2 = xpoly(ipoly+1);
    y1 = ypoly(ipoly);
    y2 = ypoly(ipoly+1);
    
    for i = 1:length(xline)-1
        x3=xline(i);
        x4=xline(i+1);
        y3=yline(i);
        y4=yline(i+1);
        
        p = intersect_line(x1,y1,x2,y2,x3,y3,x4,y4);
        
        if ~isempty(p) & p(1)>=min(x1,x2) & p(1)<=max(x1,x2) & p(1)>=min(x3,x4) & p(1)<=max(x3,x4)
            cross = 1;
            return
        end
    end
end
end