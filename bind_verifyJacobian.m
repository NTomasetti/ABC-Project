clear
clc

rng(12345)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% David T. Frazier and Gael. M. Martin 
% ABC for MA(2) example:
% y_t=\theta_1\epsilon_{t-1}+\theta_2\epsilon_{t-2}+\epsilon_{t}
% This file verifies injectivity of the binding function by simulating a
% very long sample and evaluating the summary statistics at various points
% in the region of high concentration. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
T = 10000;                                                                 % Sample size: needs to be very large
d_theta = 2;                                                               % Number of parameters
d_eta = 4;                                                                 % Number of summary statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose appropriate parameter values to span the space and construct all
% possible combinations therein. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta1 = -.99:.4:.99; %grid for first parameter, will need to be relatively small unfortunately. 
theta2 = -.99:.4:.99; %grid for second parameter, will need to be relatively small unfortunately. 

all_comb=combvec(theta1,theta2); %all possible combinations of grid points
theta1p=all_comb(1,:)'; %Pulling first coordinate in combination
theta2p=all_comb(2,:)'; %Pulling second coordinate in combination

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulating a very long chain of errors that will be fixed throughtout the
% procedure: controls Monte Carlo variation. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nu = randn(T,1); %Simulating once the errors
nu = nu*ones(1,size(all_comb,2)); 
%The above creates a T-by-(size all combinations) vector of errors. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical evaluation of the Jacobian across different points in the
% parameter space. If Jacobian is Positive definite across all these points
% then the function must be injective. This is a sufficient condition only
% and is not necessary. 
%
% Central finite differences is used for evaluation the Jacobian. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
h = .0000001; %Step size for the numerical derivative procedure

for i=1:length(all_comb) %Parallel loop for speed purposes. 
score1 = zeros(1,4); %Using the parallel loop requires re-allocation of memory.
score2 = zeros(1,4); %Using the parallel loop requires re-allocation of memory.
  
    theta = all_comb(:,i); % Pulling one combination. 
    grid1 = [theta(1)-.5*h,theta(1)+.5*h]; % Adding and substracting first dimension for the numerical derivative: 
    score1(1,:) =( bind_MA2_simple(T,grid1(2),theta(2),nu) - bind_MA2_simple(T,grid1(1),theta(2),nu))/h; % Central numerical difference. 

    grid2 = [theta(2)-.5*h,theta(2)+.5*h]; % Adding and substracting first dimension for the numerical derivative: 
    score2(1,:) =( bind_MA2_simple(T,theta(1),grid2(2),nu) - bind_MA2_simple(T,theta(1),grid2(1),nu))/h; % Central numerical difference.
    
    Score(:,:,i) = [score1;score2]; %Combining the scores across the different combinations of grid points. The result here is a three dimensional matrix.
    
end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All Deriviatives have been obtained. The resul is a
% Dim(theta)xDim(Summaries) matrix. We need to form all possible
% Dim(theta)xDim(theta) combinations of this matrix. Then, we must check if
% each of these combinations is positive definite. Restrict ourself to
% square matrices b/c that is the only way to apply the sufficient
% condition. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 b = nchoosek(1:1:size(Score,2),d_theta); % All possible combinations
 
% Need all possible combinations of a Dim(theta)xDim(theta) matrix over the
% grid of parameters and across the different summary statistic
% combinations. I've chosen to do this by constructing a four dimensional
% matrix, but other ways are possibly easier and do not require such
% abstract constructions...

Mats2x2 = zeros(length(b),d_theta,d_theta,length(all_comb)); 
% Generating the four-dimensional matrix by pulling the correct elements
% from the score function
    for i=1:size(b,1)
         bb = b(i,:);
         Mats2x2(i,:,:,:) = [Score(:,bb(1),:) Score(:,bb(2),:)]; 
     
    end




for i=1:size(b,1)
   
temp1 = squeeze(Mats2x2(i,:,:,:)); % temporary matrix to house the requesed Dim(theta)xDim(theta) matrix. 

    for j=1:length(all_comb)
        pos_def(j) = all(eig(temp1(:,:,j))>0); % Determining if the temporary matrix is positiive definite over all different parameter combinations.
    end
% Recording whether or not the Jacobian is pos. def. over the different
% parameter combinations. This process is repeated over the different
% possible summary combinations. 

pp(i,:) = pos_def; 
end


% The sufficient condition is that over any value of the parameters is a
% convex set, if the Jacobian is positive definite, then the function is
% injective on the convex set. This means that, over all the points we've
% tried, the function will be injective if over all the combinations, one
% of the submatrices is always positive definite. 

% To find this, we simply % sum the entries of pp across the number of
% parameter combinations dimension. 
PD = sum(pp,2);

% Finds if any of the matrices was pos. def. over all parameter values.
[index] = find(PD==length(all_comb));

% yields the indices of b that satisfy the sufficient condition. 
Satisfaction = b(index,:)
% this entry could be empty!

toc
