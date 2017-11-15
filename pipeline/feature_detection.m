names = 'ABCD';
for n = 552:552
    for k = 1:4
        n1 = strcat('images/', num2str(n),'/', names(k), '.JPG');
        n2 = strcat('images/', num2str(n + 1),'/', names(k), '.JPG');
        I1 = rgb2gray(imread(n1));
        I2 = rgb2gray(imread(n2));

        points1 = detectSURFFeatures(I1);
        points2 = detectSURFFeatures(I2);
        
        [features1,valid_points1] = extractFeatures(I1,points1);
        [features2,valid_points2] = extractFeatures(I2,points2);
        
        indexPairs = matchFeatures(features1,features2);

        matchedPoints1 = valid_points1(indexPairs(:,1),:);
        matchedPoints2 = valid_points2(indexPairs(:,2),:);
        
        figure;
        ax = axes; 
        
        showMatchedFeatures(I1,I2,matchedPoints1,matchedPoints2, 'montage' ,'Parent',ax);
        title(ax, 'Candidate point matches');
        legend(ax, 'Matched points 1','Matched points 2');
        points_correspondence_dir = './points_correspondence/';
        if (isdir(points_correspondence_dir) == 0)
            mkdir(points_correspondence_dir)
        end
        
        
        p =  strcat(points_correspondence_dir, names(k), '/');
        
        if (isdir(p) == 0)
            mkdir(p)
        end
        
        n1 = strcat('points_',num2str(n), '_left');
        n2 = strcat('points_',num2str(n), '_right');
        
        dlmwrite(strcat(p,n1), matchedPoints1.Location, 'delimiter','\t');
        dlmwrite(strcat(p,n2), matchedPoints2.Location, 'delimiter','\t');
        
       
    end
end    
    
