lambda  =  -0.789792;
lambda2 = -0.115003;
for n = 553:553
     %n1 = strcat('./full_undistorted_images/', num2str(n),'/', 'B_undistorted_k0_is_', num2str(lambda,6), '.jpg');
     %n2 = strcat('./full_undistorted_images/', num2str(n + 1),'/', 'B_undistorted_k0_is_', num2str(lambda,6), '.jpg');
     n1 = strcat('./full_undistorted_images/', num2str(n),'/', 'C_', 'und', '.JPG');
     n2 = strcat('./full_undistorted_images/', num2str(n+1),'/', 'C_', 'und', '.JPG');
     %Irgb1 = imread(n1);
     %Irgb2 = imread(n2);
     %I1 = rgb2gray(Irgb1);
     %I2 = rgb2gray(Irgb2);
     %points1 = detectSURFFeatures(I1);
     %points2 = detectSURFFeatures(I2);

      
        
     %[features1,valid_points1] = extractFeatures(I1,points1);
     %[features2,valid_points2] = extractFeatures(I2,points2);
        
     %indexPairs = matchFeatures(features1,features2);

     %matchedPoints1 = valid_points1(indexPairs(:,1),:);
     %matchedPoints2 = valid_points2(indexPairs(:,2),:);
     %[fLMedS, inliers] = estimateFundamentalMatrix(matchedPoints1,matchedPoints2,'Method', 'LTS', 'NumTrials',2000, 'InlierPercentage', 20);
     %matchedPoints1 = matchedPoints1(inliers,:);
     %matchedPoints2 = matchedPoints2(inliers,:);
     %fLMedS = fLMedS';
     %u2 = (matchedPoints2.Location)';
     %u1 = (matchedPoints1.Location)';
     %sz1 = size(u1);
     %sz2 = size(u2);
     %u1 = [u1; ones(1, sz1(2))];
     %u2 = [u2; ones(1, sz2(2))];
     %ll1 = (fLMedS*u2);
     %ll2 = (u1'*fLMedS)';
     %errors = (ones(1, sz1(2))./(ll1(1,:).^2 + ll1(2,:).^2) + ones(1, sz2(2))./(ll2(1,:).^2 + ll2(2,:).^2))'.*(diag(u1'*fLMedS*u2)).^2;
     %figure;
     %showMatchedFeatures(I1, I2, matchedPoints1,matchedPoints2,'montage','PlotOptions',{'ro','go','y--'});
     %title('Matlab Point matches after outliers were removed');
    
     myF =  dlmread(strcat('/home/danielbord/CLionProjects/AutomaticSolver/pipeline/', 'fundamental_matrices_C_2'));
     myF = myF(4:6, 1:3);
    
     mpd1 = dlmread(strcat('/home/danielbord/CLionProjects/AutomaticSolver/pipeline/automatic_solver_results_C/',num2str(n),'_left'));
     mpd2 = dlmread(strcat('/home/danielbord/CLionProjects/AutomaticSolver/pipeline/automatic_solver_results_C/',num2str(n),'_right'));
     w = 7360.0; h = 4912.0;
     r =sqrt(w^2/4 + h^2/4);
     szz1 = size(mpd1);
     mpd1(:,1) = mpd1(:,1) - ones(szz1(1), 1)*w/2;
     mpd1(:,2) = mpd1(:,2) -  ones(szz1(1), 1)*h/2;
     mpd2(:,1) = mpd2(:,1) -  ones(szz1(1), 1)*w/2;
     mpd2(:,2) = mpd2(:,2) -  ones(szz1(1), 1)*h/2;
     mpd1 = mpd1/r;
     mpd2 = mpd2/r;
     d = 7360;
     al = (16 + 4*lambda * (d / r) * (d / r) + lambda2*(d/r)^4)/16;
     rd1 = (sum(mpd1.^2,2));
     rd2 = (sum(mpd2.^2,2));
     mp1 = mpd1;
     mp2 = mpd2;
     mp1(:, 1) = al*r*mpd1(:, 1)./( ones(szz1(1), 1) + lambda*rd1) +  ones(szz1(1), 1)*w/2;
     mp1(:, 2) = al*r*mpd1(:, 2)./( ones(szz1(1), 1) + lambda*rd1)+  ones(szz1(1), 1)*h/2;
     
     mp2(:, 1) = al*r*mpd2(:, 1)./( ones(szz1(1), 1) + lambda*rd2) +  ones(szz1(1), 1)*w/2;
     mp2(:, 2) = al*r*mpd2(:, 2)./( ones(szz1(1), 1) + lambda*rd2) + ones(szz1(1), 1)* h/2;
    % figure;
    % k =25;
     %showMatchedFeatures(Irgb1, Irgb2, mp1(1:1+k,:),mp2(10:10+k,:),'montage','PlotOptions',{'ro','go','y--'});
     %title('My undistorted points');
     %epiLines = epipolarLine(myF', mp2(1:1+k,:));
     %points = lineToBorderPoints(epiLines,size(I1));
     %line(points(:,[1,3])',points(:,[2,4])');
     
   
     
     sz1 = size(mp1);
     sz2 = size(mp2);
     mp1 = [mp1'; ones(1, sz1(1))];
     mp2 = [mp2'; ones(1, sz2(1))];
     u1 = mp1;
     u2 = mp2;
     myF = myF';
     ll1 = (myF*u2);
     ll2 = (u1'*myF)';
     errors = (ones(1, sz1(1))./(ll1(1,:).^2 + ll1(2,:).^2) + ones(1, sz2(1))./(ll2(1,:).^2 + ll2(2,:).^2))'.*(diag(u1'*myF*u2)).^2;

end
