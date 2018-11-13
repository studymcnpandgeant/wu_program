clear all
fidin=fopen('material.dat'); 
N1=6; 
N2=6; 
N3=180; 
DATA=zeros(N1,N2,N3);
% while ~feof(fidin)   
% tline=fgetl(fidin);  
% tline=fgetl(fidin);      
for i=1:1:N3
           tline=fgetl(fidin);      
           for k=1:1:N1
             tline=fgetl(fidin);       
             str=tline;     
             sline=sscanf(str,'%f',N2); 
             DATA(k,:,i)=sline;  
           end 
   end
% end
fclose(fidin);
N=N1*N2*N3;
vecdata = [1:1:N];
for i=1:1:N3
   for k=1:1:N1
      for j=1:1:N2
          index = (i-1)*36 +(k-1)*6+j;
         vecdata(index) = DATA(k,j,i); 
      end
   end
end

data1=load('corereal.dat',' ');
NN=16*35;
veccore = [1:1:NN];
for i=1:1:35
    for j=1:1:16
       index2 = (i-1)*16 + j;
       veccore(index2) = data1(i,j);
    end
end

data2=load('r.dat',' ');
R=18;
vecr = [1:1:R];
for i=1:1:3
    for j=1:1:6
       index3 = (i-1)*6 +j;
       vecr(index3) = data2(i,j);
    end
end

data3=load('z.dat',' ');
Z=36;
vecz = [1:1:Z];
for i=1:1:6
    for j=1:1:6
       index4 = (i-1)*6 +j;
       vecz(index4) = data3(i,j);
    end
end

dx = [1:1:18];
dx(1) = 10;
for i=2:1:17
   dx(i) = vecr(i) - vecr(i-1); 
end
dx(18) = 10;

dy = [1:1:37];
dy(1) = 10;
for i=2:1:36
   dy(i) = vecz(i)-vecz(i-1); 
end
dy(37) = 10;

PetscBinaryWrite('wu-sn2dtest6',vecdata,veccore,dx,dy,vecr,vecz, 'indices','int32','precision','float64');