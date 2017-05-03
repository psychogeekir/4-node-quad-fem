function [NodeNumber,ElementNumber] = GetCurrentElement(x,y,NodeAssociateWithElement,NodeCoordinate)

x_unique = unique(NodeCoordinate(:,1));
y_unique = unique(NodeCoordinate(:,2));

xind = find(x_unique >= x,1);
yind = find(y_unique >= y,1);

if (xind > 1 && yind > 1)
    NodeNumber = find((NodeCoordinate(:,1)==x_unique(xind) | NodeCoordinate(:,1)==x_unique(xind-1))...
        & (NodeCoordinate(:,2)==y_unique(yind) | NodeCoordinate(:,2)==y_unique(yind-1)));
elseif (xind == 1 && yind > 1)
    NodeNumber = find((NodeCoordinate(:,1)==x_unique(xind) | NodeCoordinate(:,1)==x_unique(xind+1))...
        & (NodeCoordinate(:,2)==y_unique(yind) | NodeCoordinate(:,2)==y_unique(yind-1)));
elseif (xind > 1 && yind == 1)
    NodeNumber = find((NodeCoordinate(:,1)==x_unique(xind) | NodeCoordinate(:,1)==x_unique(xind-1))...
        & (NodeCoordinate(:,2)==y_unique(yind) | NodeCoordinate(:,2)==y_unique(yind+1)));
elseif (xind == 1 && yind == 1)
    NodeNumber = find((NodeCoordinate(:,1)==x_unique(xind) | NodeCoordinate(:,1)==x_unique(xind+1))...
        & (NodeCoordinate(:,2)==y_unique(yind) | NodeCoordinate(:,2)==y_unique(yind+1)));
end

ElementNumber = find(ismember(sort(NodeAssociateWithElement)',NodeNumber','rows'));

end