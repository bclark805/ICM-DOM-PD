function [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out,Y_out);
% Function to calculate statistics for biogeochemical models
%input X (observed data) Y (modeled output) to get stats and axis
%dimensions for plotting
%B Clark UMCES, Jan 2019
%MEF   
    part1 = sum((X_out-mean(X_out)).^2);
    part2 = sum((Y_out-X_out).^2);
    
       MEF =(part1-part2)./part1;     
 % r   

     bottom = sqrt(sum((X_out-mean(X_out)).^2).*sum((Y_out-mean(Y_out)).^2));
    top=sum((X_out-mean(X_out)).*(Y_out-mean(Y_out)));

r=top./bottom;

%Wilmott skill
    top=sum((Y_out-X_out).^2);
    bottom= sum((abs(Y_out-mean(X_out))+abs(X_out - mean(X_out))).^2);
    WMS = 1-top./bottom;

%rmse
RMSE = sqrt((sum((Y_out-X_out).^2))./length(X_out));

   axis_min=min(min([X_out Y_out])) ;
   axis_max=max(max([X_out Y_out]));

   %mean percent error
   
   MPE=mean((X_out-Y_out)./X_out)*100;
   
   %reliability index
   
   RI = exp((mean((log(X_out./Y_out)).^2)).^0.5);
   
  
end