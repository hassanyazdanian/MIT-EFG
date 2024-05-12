function v = pgauss(k)
% This function returns a matrix with 4 gauss points and their weights
if k==4
   v(1,1) =-.861136311594052575224;
   v(1,2) =-.339981043584856264803;
   v(1,3) =-v(1,2);
   v(1,4) = -v(1,1);
   v(2,1) =.347854845137453857373;
   v(2,2) =.652145154862546142627;
   v(2,3) = v(2,2);
   v(2,4) = v(2,1);
end
