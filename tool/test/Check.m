function tag = Check(a, b, c, p)
% global pointA pointB pointC;
ab = b - a;
ac = c - a;
ap = p - a;
L(1) = dot(ab,ab);
L(2) = dot(ab,ac);
L(3) = dot(ac,ac);
L(4) = L(2)*L(2) - L(1)*L(3);

RV1 = dot(ap,ab);
RV2 = dot(ap,ac);
RV3 = RV1 - L(1);
RV4 = RV2 - L(2);
RV5 = RV1 - L(2);
RV6 = RV2 - L(3);
RVu = (RV2*L(2)-RV1*L(3))/L(4);
RVv = (RV1*L(2)-RV2*L(1))/L(4);

Dx = 1.0;
if ~and(and(RVu>=-0.8*Dx/sqrt(L(1)),RVv>=-0.8*Dx/sqrt(L(3))),(RVu+RVv)<=(1+0.8*Dx/sqrt(L(1)+0.8*Dx/sqrt(L(3)))))
    tag = 8;
else
    % Corner A
    if (and(RV1<=0,RV2<=0))
        tag = 1;
        return;
    end
    
    % Corner B
    if (and(RV3>=0,RV4<=RV3))
        tag = 2;
        return;
    end
    
    % Corner C
    if (and(RV6>=0,RV5<=RV6))
        tag = 3;
        return;
    end
    
    % Edge AB
    if (and(and(RV1>0,RV3<0),RVv<0))
        tag = 4;
        return;
    end
    
    % Edge BC
    if (and(RV4>RV3,and(RV5>RV6,RVu+RVv>1)))
        tag = 5;
        return;
    end
    
    % Edge AC
    if (and(and(RV2>0,RV6<0),RVu<0))
        tag = 6;
        return;
    end
    
    % Inner
    if (and(RVu+RVv<=1,and(RVu>=0,RVv>=0)))
        tag = 7;
        return;
    end
end
% if u >= 0
%     if v <= 0
%         if u+v < 1
%             disp('region 1');
%             tag = 1;
%         else
%             disp('region 2');
%             tag = 2;
%         end
%     else
%         if u+v < 1
%             disp('region 3');
%             tag = 3;
%         else
%             disp('region 4');
%             tag = 4;
%         end
%     end
% else
%     if v <= 0
%         disp('region 5');
%         tag = 5;
%     else
%         if u+v < 1
%             disp('region 6');
%             tag = 6;
%         else
%             disp('region 7');
%             tag = 7;
%         end
%     end 
end
% tag = 2;
% Dx = 0.5;
% if u>(-0.8*Dx/norm(pointC-pointA))
%     if v>(-0.8*Dx/norm(pointB-pointA))
%         %         if u*norm(pointC-pointA)/(norm(pointC-pointA)+1.6*Dx)+v*norm(pointB-pointA)/(norm(pointB-pointA)+1.6*Dx) < 1
%         if u+v<1+0.8*Dx/norm(pointB-pointA)+0.8*Dx/norm(pointC-pointA)
%             tag = 1;
%         end
%     end
% end
% 
% end
% .AND. ((RVu+RVv).LE.(1.0+0.8*Dx/sqrt(OMP_trimeshWallLength(2))+0.8*Dx/sqrt(OMP_trimeshWallLength(4)))
