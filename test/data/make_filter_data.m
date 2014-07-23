
% in Octave, using Wavelab850

addpath('/home/gfa/forr/octave/Wavelab850/Orthogonal')

name = 'filter';
n = 64;
x = randn(n,1);
L = 0;
n2 = 8;
x2 = randn(n2,n2);
L2 = 0;

wname = {'Daubechies','Coiflet','Haar','Symmlet','Battle','Vaidyanathan','Beylkin'};
wnum = {[4:2:20],[2:5],[0],[4:10],[1,3,5],[0],[0]};

for i = 1:length(wname)
    for num = 1:length(wnum{i})
        qmf = MakeONFilter(wname{i},wnum{i}(num));
        if strcmp(wname{i},'Symmlet')
            qmf=reverse(qmf);
        end
        y = FWT_PO(x,L,qmf);
        y2 = FWT2_PO(x2,L2,qmf);
        
        %norm(IWT_PO(y,L,qmf)-x)
        % 1-d
        fid=fopen([name, '1d_', wname{i}, num2str(wnum{i}(num)),'.txt'], 'w');
        fprintf(fid,"%.16e\n",y);
        fclose(fid);
        % 2-d
        fid=fopen([name, '2d_', wname{i}, num2str(wnum{i}(num)),'.txt'], 'w');
        dlmwrite(fid,y2,"delimiter", "\t", "newline", "\\n", "precision","%.16e");
        fclose(fid);
    end
end

fid=fopen([name, '1d_', 'data','.txt'], 'w');
fprintf(fid,"%.16e\n",x);
fclose(fid);

fid=fopen([name, '2d_', 'data','.txt'], 'w');
dlmwrite(fid,x2,"delimiter", "\t", "newline", "\\n", "precision","%.16e");
fclose(fid);


