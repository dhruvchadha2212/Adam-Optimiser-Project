function f=funcx(x)
    f=(1.5-x(1)+x(1)*x(2))^2+(2.25-x(1)+x(1)*x(2)^2)^2+(2.625-x(1)+x(1)*(x(2)^3))^2;
endfunction

function g=grad1(x)  
g=[2*((1.5-x(1)+x(1)*x(2))*(-1+x(2))+(2.25-x(1)+x(1)*x(2)^2)*(-1+x(2)^2)+(2.625-x(1)+x(1)*(x(2)^3))*(-1+x(2)^3)),2*((1.5-x(1)+x(1)*x(2))*x(1)+(2.25-x(1)+x(1)*x(2)^2)*(2*x(1)*x(2))+(2.625-x(1)+x(1)*(x(2)^3))*(3*x(1)*x(2)^2))];
endfunction

xprev=[2,3];
xnew=[3,2];

dim=2;
G=zeros(1,dim);
prev=.00000001;
step_sz=.1;
counter=0;

  while(abs(funcx(xprev)-funcx(xnew))>.000000001)
      xprev=xnew;
      grad=grad1(xprev);
      
      for i=1:dim
          G(i)=G(i)+grad(i)*grad(i)';
          xnew(i)=xprev(i)-(step_sz/(sqrt(G(i))+prev))*grad(i);
      end
      
      counter=counter+1;
  end
  printf("%d\n",counter);
   printf("%2.6f\n",funcx(xprev));
