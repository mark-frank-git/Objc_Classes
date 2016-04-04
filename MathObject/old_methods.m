
int level,maxlev,levmax,feval;   
/*##################################*
  Adaptive integration routine      *
  uses simpson's rule adaptively    *
  from Handbook of C Tools for      *
  Scientists and Engineers by       *
  L. Baker                          *
 *##################################*/
- (double) integrate1From:(double) ax to:(double)bx withTol:(double)tol from:(id)sender
{
  double a,b,eps;
  double aa,epss,absarea,ans,fa,fb,fm,range,est;
  double mid;
  
  a = ax;
  b = bx;
  eps = tol;
  if(![ sender respondsTo:sel_getUid("functionToIntegrate:") ])
  {
    printf("Can't call functionToIntegrate in integrateFrom:\n");
    return 0.;
  }
  
  levmax=1;
  feval=3;
  level=1;
  maxlev=50;
  aa=a;epss=eps;
  absarea=1.;
  est=1.;
  range=b-a;
  mid = .5*(a+b);
  fa = [sender functionToIntegrate:a];
  fb = [sender functionToIntegrate:b];
  fm = 4.*[sender functionToIntegrate:mid];
  ans=[self simpsons:aa range:range fa:fa fm:fm fb:fb area:absarea estimate:est
      withTol:epss from:sender];
  return ans;
}

   
/*##################################*
  Simpson integration routine       *
  from Handbook of C Tools for      *
  Scientists and Engineers by       *
  L. Baker                          *
 *##################################*/
- (double) simpsons:(double)a range:(double)da  fa:(double)fa fm:(double)fm fb:(double)fb
           area:(double)area estimate:(double)est withTol:(double)eps from:(id)sender
{
  double absarea, arg;
  double dx,x1,x2,est1,est2,est3,f1,f2,f3,f4,sum,epss,norm=.588;
  absarea=area;

  if(![ sender respondsTo:sel_getUid("functionToIntegrate:") ])
  {
    printf("Can't call functionToIntegrate in simpsons:\n");
    return 0.;
  }
  dx=.333333333*da;
  epss=eps*norm;
  x1=a+dx;
  x2=x1+dx;
  arg = a + 0.5*dx;
  f1 = 4.*[sender functionToIntegrate:arg];
  f2=[sender functionToIntegrate:x1];
  f3=[sender functionToIntegrate:x2];
  arg = a + 2.5*dx;
  f4 = 4.*[sender functionToIntegrate:arg];
  feval += 4;
  est1=(fa+f1+f2)*dx*.166666666666;
  est2=(f2+fm+f3)*dx*.166666666666;
  est3=(f3+f4+fb)*dx*.166666666666;
  absarea=area-ABS(est)+ABS(est1)+ABS(est2)+ABS(est3);
  sum=est1+est2+est3;
  level++;
  if(level > levmax)
  {
    levmax = level;
  }
  if ( ((ABS(est-sum)> eps*absarea)||(est==1.)) && (level<maxlev))
  {
     sum= [self simpsons:a range:dx fa:fa fm:f1 fb:f2 area:absarea estimate:est1
           withTol:epss from:sender]
         +[self simpsons:x1 range:dx fa:f2 fm:fm fb:f3 area:absarea estimate:est2
           withTol:epss from:sender]
         +[self simpsons:x2 range:dx fa:f3 fm:f4 fb:fb area:absarea estimate:est3
           withTol:epss from:sender];
  }
 level--;
  return sum;
}
