function mx = THG_find_local_max_1D_20150117(data,tolerance)

mx = [];
cnt  = 1;
% start at tolerance, to ensure that algorithm works. This also means that
% only data point falling into the tolerance at the beginning and end will
% be checked
for j = 1+tolerance : length(data)-tolerance
    
    check1 = sum(data(j) > data(j-tolerance:j-1));
    check2 = sum(data(j) > data(j+1:j+tolerance));
    
    if check1 + check2 == tolerance * 2
       
        mx(cnt) = j;
        cnt = cnt + 1;
       
    end
    
    clear check1 check2
    
end;
        

% THG 17/01/2015

% JQK: find peak with tolerance specifing the amount of adjacent data
% points that have to be exceeded