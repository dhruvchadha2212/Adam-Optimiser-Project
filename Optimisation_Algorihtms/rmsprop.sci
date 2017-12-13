function f=funcx(x)
    f=(1.5-x(1)+x(1)*x(2))^2+(2.25-x(1)+x(1)*x(2)^2)^2+(2.625-x(1)+x(1)*(x(2)^3))^2;
endfunction
function g=grad1(x)
    g=[2*((1.5-x(1)+x(1)*x(2))*(-1+x(2))+(2.25-x(1)+x(1)*x(2)^2)*(-1+x(2)^2)+(2.625-x(1)+x(1)*(x(2)^3))*(-1+x(2)^3)),2*((1.5-x(1)+x(1)*x(2))*x(1)+(2.25-x(1)+x(1)*x(2)^2)*(2*x(1)*x(2))+(2.625-x(1)+x(1)*(x(2)^3))*(3*x(1)*x(2)^2))];
endfunction
xprev=[-10,-5];
xnew=[4,4];
step_sz=.001;
G=0;
prev=.00000001;
counter=0;
  while(abs(funcx(xprev)-funcx(xnew))>.000000001)
      xprev=xnew;
      g=grad1(xprev);
      G=0.9*G+0.1*g*g';
      xnew=xprev-((step_sz*g)/(sqrt(G+prev)));
      counter=counter+1;
  end;
  printf("%d\n",counter);
  printf("%2.6f\n",funcx(xnew));
