function mat2 = matResize(mat1,szMat2)
% setup dimensions and length 
szMat1  = size(mat1);
lenSz1  = length(szMat1);
lenSz2  = length(szMat2);
if(lenSz1<lenSz2)
    mat1    = repmat(mat1,[ones(1,lenSz1),szMat2(lenSz1+1:end)]);
    szMat1  = size(mat1);
    lenSz1  = length(szMat1);
end
% fill newsize with oldsize
szMat2n = szMat1;
szMat2n(1:lenSz2) = szMat2(1:lenSz2);
szMat2  = szMat2n;
% new matrix
if(islogical(mat1))
    mat2    = zeros(szMat2)>0;
else
    mat2    = zeros(szMat2,'like',mat1);
end
% setup index
szMatM  = min(szMat1,szMat2);
for id = 1:lenSz1
    idx1{id} = cenInd(szMat1(id),szMatM(id));
    idx2{id} = cenInd(szMat2(id),szMatM(id));
end
% command 
cmd1 = 'mat1(idx1{1}';
cmd2 = 'mat2(idx2{1}';
for id=2:lenSz1
    cmd1 = sprintf('%s,idx1{%d}',cmd1,id);
    cmd2 = sprintf('%s,idx2{%d}',cmd2,id);
end
cmd1 = sprintf('%s)',cmd1);
cmd2 = sprintf('%s)',cmd2);
cmd  = [cmd2,'=',cmd1,';'];
% eval operation
eval(cmd);




