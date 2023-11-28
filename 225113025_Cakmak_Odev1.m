%% Odev1 225113025 Yusuf Ziya Çakmak

A=[1 1 1;2 2 2;3 3 3];

s=3;
cosofA=CosMatrix(A,s,6);

%pi değerini hesaplar
function result = p_i(x)
    result = ((-1)^x)/(factorial(2*x)); 
end

%k değerine göre mk ve qk değerlerini seçer
function [mk,qk]=k(i)
    mk_table=[2,4,6,9,12,16,20,25,30,36,42,49];
    qk_table=[2,2,3,3,4,4,5,5,6,6,7,7];
    if(i>12 || i==0)
        error('0<k<=12 araliginda bir deger seciniz');
    else
        mk=mk_table(i);
        qk=qk_table(i);
    end
end

%mk değerine göre pi değerlerini oluşturur
function filledArray=fill_pi_values(mk)
    temp_Array=zeros(1,mk);
    for i=1:1:mk
        temp_array(i)=p_i(i);
    end
    filledArray=temp_array;
end

%B'nin qk ya göre bütün kuvvetlerini oluşturur
function B_powers=fill_B_powers(B,qk)
    [row,col]=size(B);
    temp_Array=zeros(row,col*qk);
    j=1;
    for i=1:col:qk*col
        temp_array(1:row,i:i+col-1)=B^j;
        j=j+1;
    end
    B_powers=temp_array;
end

%Pmk fonksiyonunu hesaplar
function Pmk_res=Pmk(B,mk,qk,pi_values,B_powers)
    [row,col]=size(B);
    identityMatrix=eye(row);
    tempMatrix=zeros(row);
    counter=mk;
    while(counter>=1)
        for j=qk:-1:1
            jcol=j*col;
            if(counter==mk && j==qk)
                tempMatrix=tempMatrix+B_powers(1:row,jcol-col+1:jcol)*pi_values(counter);
                counter=counter-1;
            else
                tempMatrix=tempMatrix+B_powers(1:row,jcol-col+1:jcol)*pi_values(counter);
                counter=counter-1;
            end
        end
        if(counter==0)
            tempMatrix=tempMatrix+identityMatrix*1;
        else
            tempMatrix=tempMatrix+identityMatrix*pi_values(counter);
            tempMatrix=cross(tempMatrix,B_powers(1:row,qk*col-col+1:qk*col));
        end     
    end
    Pmk_res=tempMatrix;
end

%Yarım açı formüllerinden geri gelir
function value=revertHalf(B,s)

    [row,col]=size(B);
    %C=zeros(row);
    C=B;
    identityMatrix=eye(row);
    for i=1:s
        C=2*C-identityMatrix;
    end
    value=C;
end

%Bütün tanımlı fonksiyonları kullanarak matrisin cosinusunu hesaplar
function val=CosMatrix(A,s,kval)
    [mk,qk]=k(kval);
    pi_values=fill_pi_values(mk);
    B=A^2;
    B=B/(4^s);
    B_powers=fill_B_powers(B,qk);
    B=Pmk(B,mk,qk,pi_values,B_powers);
    rv=revertHalf(B,s);
    val=rv;
end