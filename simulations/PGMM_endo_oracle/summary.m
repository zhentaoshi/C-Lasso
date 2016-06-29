K = 3;
p = 2;
n_lambda = 5;
Rep = size(myRep,1);

A = zeros( [ K p n_lambda Rep ] );
A_corr = A;
A_post = A;
A_post_corr = A;

for i = 1:Rep
    for j = 1:n_lambda
        A(:, :, j, i) =  myRep(i,j).H.a;
        A_post(:,:,j,i) = myRep(i,j).H_post.post_a;
        A_corr(:, :, j, i) =  myRep(i,j).H_post.a_corr;
        A_post_corr(:, :, j, i) = myRep(i,j).H_post.post_a_corr;
    end
end

mean([A A_corr], 4)

mean([A_post A_post_corr], 4)

mean(myRepClass, 3)

