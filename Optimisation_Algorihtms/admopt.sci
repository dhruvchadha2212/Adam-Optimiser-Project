function f=funcx(x)
f=(1.5-x(1)+x(1)*x(2))^2+(2.25-x(1)+x(1)*x(2)^2)^2+(2.625-x(1)+x(1)*(x(2)^3))^2;
endfunction

function g=grad(x)
     g=[2*((1.5-x(1)+x(1)*x(2))*(-1+x(2))+(2.25-x(1)+x(1)*x(2)^2)*(-1+x(2)^2)+(2.625-x(1)+x(1)*(x(2)^3))*(-1+x(2)^3)),2*((1.5-x(1)+x(1)*x(2))*x(1)+(2.25-x(1)+x(1)*x(2)^2)*(2*x(1)*x(2))+(2.625-x(1)+x(1)*(x(2)^3))*(3*x(1)*x(2)^2))];
     endfunction

xprev=[-10,-5];
xnew=[4,4];

beta1=.9;
beta2=.99;
alpha=.01;
t=0;
mprev=0;
vprev=0;
prev=0.00000001;
counter=0

  while(abs(funcx(xnew)-funcx(xprev))>0.000000001)
      counter=counter+1;
      xprev=xnew;
      t=t+1;
      g=grad(xprev);
      mnew=beta1*mprev+(1-beta1)*g;
     // disp(g);
      vnew=beta2*vprev+(1-beta2)*(g*g');
      mcorr=mnew/(1-(beta1^t));
      vcorr=vnew/(1-(beta2^t));
      xnew=xprev-alpha*mcorr/(sqrt(vcorr)+prev);
      vprev=vnew;
      mprev=mnew;
  end;
  
  printf("%d\n",counter);
  printf("%2.6f\n",funcx(xprev));
