function [a0_perm] = perm_a_K(a_K, a0_1)


I = length(a_K);
a0_perm = zeros(1,I);

if I == 2
    [~, select1] = min( abs(  a_K(1) - a0_1 )  );
    a0_perm(1) = a0_1(select1);
    a0_1(select1) = [];
    [~, select2] = min (abs( a_K(2) - a0_1) );
    a0_perm(2) = a0_1(select2);
elseif I == 3
    [~, select1] = min( abs(  a_K(1) - a0_1 )  );
    a0_perm(1) = a0_1(select1);
    a0_1(select1) = [];
    [~, select2] = min( abs(  a_K(2) - a0_1 )  );
    a0_perm(2) = a0_1(select2);
    a0_1(select2) = [];
    a0_perm(3) = a0_1;

elseif I == 4
    [~, select1] = min( abs(  a_K(1) - a0_1 )  );
    a0_perm(1) = a0_1(select1);
    [~, select2] = min( abs(  a_K(2) - a0_1 )  );
    a0_perm(2) = a0_1(select2);
    a0_1(select2) = [];
    [~, select3] = min (abs( a_K(3) - a0_1) );
    a0_perm(3) = a0_1(select3);
    a0_1(select3) = [];
    a0_perm(4) = a0_1;
elseif I == 5
    [~, select1] = min( abs(  a_K(1) - a0_1 )  );
    a0_perm(1) = a0_1(select1);
    [~, select2] = min( abs(  a_K(2) - a0_1 )  );
    a0_perm(2) = a0_1(select2);
    [~, select3] = min (abs( a_K(3) - a0_1) );
    a0_perm(3) = a0_1(select3);
    a0_1(select3) = [];
    [~, select4] = min( abs( a_K(4) - a0_1) );
    a0_perm(4) = a0_1(select4);
    a0_1(select4) = [];
    a0_perm(5) = a0_1;
end
end
