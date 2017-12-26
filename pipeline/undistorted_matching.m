lambda  =-0.990042;
%lambda2 = 0.00703376;
for n = 552:552
     %n1 = strcat('./full_undistorted_images/', num2str(n),'/', 'B_undistorted_k0_is_', num2str(lambda,6), '.jpg');
     %n2 = strcat('./full_undistorted_images/', num2str(n + 1),'/', 'B_undistorted_k0_is_', num2str(lambda,6), '.jpg');
     n1 = strcat('./full_undistorted_images/', num2str(n),'/', 'A_', 'und_t', '.JPG');
     n2 = strcat('./full_undistorted_images/', num2str(n+1),'/', 'A_', 'und_t', '.JPG');
     Irgb1 = imread(n1);
     Irgb2 = imread(n2);
     myF =  dlmread(strcat('/home/danielbord/CLionProjects/AutomaticSolver/pipeline/automatic_solver_results_A/',num2str(n),'_recomp_F'));
     myF = myF(1:3, 1:3);
    
    
     
    
     mpd1 = dlmread(strcat('/home/danielbord/CLionProjects/AutomaticSolver/pipeline/automatic_solver_results_A/',num2str(n),'_left'));
     mpd2 = dlmread(strcat('/home/danielbord/CLionProjects/AutomaticSolver/pipeline/automatic_solver_results_A/',num2str(n),'_right'));
     
     
     
     u1 = dlmread(strcat('/home/danielbord/CLionProjects/AutomaticSolver/pipeline/automatic_solver_results_A/',num2str(n),'_undistorted_left'));
     u2 = dlmread(strcat('/home/danielbord/CLionProjects/AutomaticSolver/pipeline/automatic_solver_results_A/',num2str(n),'_undistorted_right'));
     

     w = 7360.0; h = 4912.0;
     r =sqrt(w^2/4 + h^2/4);
    
    
     dImg = imread(strcat('/home/danielbord/CLionProjects/AutomaticSolver/pipeline/images/', num2str(n), '/A.JPG'));
      figure;
     imshow(dImg);
     hold on;
     plot(mpd1(:,1), mpd1(:,2), 'r.', 'LineWidth', 2, 'MarkerSize', 3);
     rret_mpd1 = zeros(size(mpd1));
     rret_mpd2 =  zeros(size(mpd2));
     epiLines = epipolarLine(myF', u2);
     qqq = size(u1);
     d_r = 7360;
     r_img = sqrt((7360/2)^2 + (4912/2)^2);
     dd = zeros(qqq(1));
     %closest_point2 = zeros(qqq);
    %  for kk = 1:qqq(1)
   %       dist = abs(epiLines(kk,1)*u1(kk, 1) + epiLines(kk,2)*u1(kk, 2) + epiLines(kk,3))/sqrt(epiLines(kk,1)^2 + epiLines(kk,2)^2);
   %       vec = [epiLines(kk,1), epiLines(kk,2)];
   %       vec = vec/norm(vec);
  %        closest_point = u1(kk, :) +  dist*vec;
   %       closest_point2(kk,:) = closest_point;
  %        dd(kk) = dist;
  %        alpha = 4/(4+lambda*(d_r/r_img)^2);
  %        closest_point(1) = alpha*(closest_point(1) - 7360/2)/r_img;
   %        closest_point(2) =alpha*(closest_point(2) - 4912/2)/r_img;
   %        m_r_u =sum(closest_point.^2);
  %         m_r_d = (1 - sqrt(1 - 4 * lambda * m_r_u))/ (2 * lambda * sqrt(m_r_u));
  %         closest_point(1) = r_img*closest_point(1)*(1 + lambda*m_r_d^2) +  7360/2;
  %         closest_point(2) = r_img*closest_point(2)*(1 + lambda*m_r_d^2) +  4912/2;
  %         rret_mpd1(kk,1) = closest_point(1);
  %         rret_mpd1(kk,2) = closest_point(2);
  %         
   %  end
       rret_mpd1 = dlmread(strcat('/home/danielbord/CLionProjects/AutomaticSolver/pipeline/automatic_solver_results_A/',num2str(n),'_curve_left'));
       rret_mpd2 =  dlmread(strcat('/home/danielbord/CLionProjects/AutomaticSolver/pipeline/automatic_solver_results_A/',num2str(n),'_curve_right'));
     plot(rret_mpd1(:,1), rret_mpd1(:,2), 'g.', 'LineWidth', 2, 'MarkerSize', 3);
     for kk = 1:qqq(1)
          line([mpd1(kk, 1), rret_mpd1(kk,1)],[mpd1(kk,2), rret_mpd1(kk,2)]);
     end
   
     figure;
     showMatchedFeatures(Irgb1, Irgb2, u1,u2,'montage','PlotOptions',{'ro','go','y--'});
     title('My undistorted points');
     epiLines = epipolarLine(myF', u2);
     points = lineToBorderPoints(epiLines,size(Irgb1));
     %line(points(:,[1,3])',points(:,[2,4])');
     
     epiLines2 = epipolarLine(myF, u1);
     shift = ones(size(u2));
     shift(:, 1) = w;
     shift(:, 2) = 0;
     points2 = lineToBorderPoints(epiLines2,size(Irgb2)) + [shift, shift];
     %line(points2(:,[1,3])',points2(:,[2,4])');

end
