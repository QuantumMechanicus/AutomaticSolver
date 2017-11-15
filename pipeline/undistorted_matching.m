addpath '/home/danielbord/Documents/MATLAB/Italyanskaya_street/'
for n = 552:552
     n1 = strcat('pictures/', num2str(n),'/', '_undistorted_l_is_-0.854862.jpg');
     n2 = strcat('pictures/', num2str(n + 1),'/', '_undistorted_l_is_-0.854862.jpg');
     I1 = rgb2gray(imread(n1));
     I2 = rgb2gray(imread(n2));
     points1 = detectSURFFeatures(I1);
     points2 = detectSURFFeatures(I2);

      
        
     [features1,valid_points1] = extractFeatures(I1,points1);
     [features2,valid_points2] = extractFeatures(I2,points2);
        
     indexPairs = matchFeatures(features1,features2);

     matchedPoints1 = valid_points1(indexPairs(:,1),:);
     matchedPoints2 = valid_points2(indexPairs(:,2),:);
     [fLMedS, inliers] = estimateFundamentalMatrix(matchedPoints1,matchedPoints2,'Method', 'LTS', 'NumTrials',2000, 'InlierPercentage', 20);
     matchedPoints1 = matchedPoints1(inliers,:);
     matchedPoints2 = matchedPoints2(inliers,:);
     fLMedS = fLMedS';
     u2 = (matchedPoints2.Location)';
     u1 = (matchedPoints1.Location)';
     sz1 = size(u1);
     sz2 = size(u2);
     u1 = [u1; ones(1, sz1(2))];
     u2 = [u2; ones(1, sz2(2))];
     ll1 = (fLMedS*u2);
     ll2 = (u1'*fLMedS)';
     errors = (ones(1, sz1(2))./(ll1(1,:).^2 + ll1(2,:).^2) + ones(1, sz2(2))./(ll2(1,:).^2 + ll2(2,:).^2))'.*(diag(u1'*fLMedS*u2)).^2;
     figure;
     showMatchedFeatures(I1, I2, matchedPoints1,matchedPoints2,'montage','PlotOptions',{'ro','go','y--'});
     title('Point matches after outliers were removed');

end
