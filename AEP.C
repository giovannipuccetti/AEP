# include <iostream>
# include <math.h>
# include <time.h>
# include <sys/time.h>
typedef double NT; //working class number type
using namespace std;


NT JDF(NT* x, NT h);
void Printresults(int N, NT *Pn, NT *PnStar, NT *Tn);

const int DIM = 3; //dimension

//loads the cubecoord, cubesign matrices
//# include "cubecoord.h"
int cubecoord2[8] ={1,1, 1,0, 0,1, 0,0};
NT cubesign2[4] = {1., -1.,-1., 1.};
int cubecoord3[24] = {1,1,1, 1,1,0, 1,0,1, 
  0,1,1, 1,0,0, 0,1,0, 0,0,1, 0,0,0};
NT cubesign3[8] = {1., -1.,-1.,-1., 1.,1.,1., -1.}; 
int cubecoord4[64]=
  {1,1,1,1, 1,1,1,0, 1,1,0,1, 1,0,1,1, 0,1,1,1, 
  1,1,0,0, 1,0,1,0, 1,0,0,1, 0,1,1,0, 0,1,0,1, 0,0,1,1, 
  1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1, 0,0,0,0};
NT cubesign4[16]=
  {1, -1,-1,-1,-1, 1,1,1,1,1,1, -1,-1,-1,-1, 1};
int cubecoord5[160]=
  {1,1,1,1,1,
  1,1,1,1,0, 1,1,1,0,1, 1,1,0,1,1, 1,0,1,1,1, 0,1,1,1,1,
  1,1,1,0,0, 1,1,0,1,0, 1,0,1,1,0, 0,1,1,1,0, 1,1,0,0,1, 
  1,0,1,0,1, 0,1,1,0,1, 1,0,0,1,1, 0,1,0,1,1, 0,0,1,1,1,
  1,1,0,0,0, 1,0,1,0,0, 1,0,0,1,0, 1,0,0,0,1, 0,1,1,0,0,
  0,1,0,1,0, 0,1,0,0,1, 0,0,1,1,0, 0,0,1,0,1, 0,0,0,1,1,
  1,0,0,0,0, 0,1,0,0,0, 0,0,1,0,0, 0,0,0,1,0, 0,0,0,0,1,
  0,0,0,0,0};
NT cubesign5[32]=
  {1, -1,-1,-1,-1,-1, 1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,  1,1,1,1,1, -1};

int* ccIT, *cubecoordITfix;
NT* cubesignIT, *cubesignITfix; 
NT a, qqh;  // JDF evaluations on the corner
NT fm1[2], fm2[2], fm3[2], fm4[2], fm5[2]; //marginal dist
int i4, pow2dim = int(pow(2,DIM));
int simplexfac[6] = {-1,-1, 3,4,15,21};
NT alphaDIM[6] = {0.,0., 2./3., 1./2., 2./5., 1./3.};
NT expolfac[6] = {0.,0., 9./8., 4./3., 625./384., 81./40.};

int main(){
  NT s = 1.5;
  int Ntot = 11;
  std::cout<< "\n Ntot: ";//number of iterations
  std::cin>> Ntot;
  std::cout<< "\n s: "; //threshold
  std::cin>> s;
  std::cout<< "\n STATE =  s: "<<s<<"  Ntot: "<<Ntot;

  NT dp, h, as , x[DIM], P = 0., alpha = alphaDIM[DIM];
  int iterationlength, ms, i, i2, i3, coordindicator = 1; // ms = matrixsize
  NT Pn[Ntot], PnStar[Ntot], Tn[Ntot]; 

  //initializing coordinate lists
  NT* list1 = new NT[DIM];       NT* coordIT = list1;
  for( i2 = 0; i2 < DIM; ++i2){*coordIT++ = 0;}
  NT* list1length = new NT[1]; NT* lengthIT = list1length;
  *list1length = s;
  NT* list1sign = new NT[1];   NT* signIT = list1length;
  *list1sign = 1;
 
  NT* list2       = new NT[1]; NT* list2length = new NT[1];
  NT* list2sign   = new NT[1];
  NT* NEWcoordIT, *NEWlengthIT, *NEWsignIT, *xIT, *xITfix = x;

  // assigns correct udmatrix and corresponding pointers acording to DIM
  //# include "udmatrices.h" 
  NT *udmatrixIT, *udmatrixITfix, *udsignIT, *udsignITfix;
  NT *udlengthIT, *udlengthITfix;
  
  // dimension 2  
  NT  udmatrix2[6] = {2./3.,2./3., 2./3.,0, 0,2./3.};
  NT  udsign2[3]   = {-1.,1.,1.};
  NT  udlength2[3] = {-1./3.,1./3., 1./3.};
  // dimension 3
  NT  udmatrix3[12]= {1./2,1./2,1./2,1./2,0,0,0,1./2,0,0,0,1./2};
  NT  udsign3[4]   = {-1.,1.,1.,1.};
  NT  udlength3[4] = {-1./2,1./2,1./2,1./2};
  // dimension 4
   h = 2./5.;  //abusing h
  NT udmatrix4[60]=  {h,h,h,h, h,h,h,0, h,h,0,h, h,0,h,h, 0,h,h,h, h,h,0,0, h,0,h,0, 0,h,h,0, h,0,0,h, 0,h,0,h, 0,0,h,h, h,0,0,0, 0,h,0,0, 0,0,h,0, 0,0,0,h};
  NT udsign4[15] = {-1, 1,1,1,1, -1,-1,-1,-1,-1,-1, 1,1,1,1};
  h = 1./5.;
  NT udlength4[15]={-3*h,-h,-h,-h,-h,h,h,h,h,h,h,3*h,3*h,3*h,3*h};
  // dimension 5
  h = 1./3.; //abusing h
  NT udmatrix5[21*5]=   {h,h,h,h,h, h,h,h,h,0, h,h,h,0,h, h,h,0,h,h, h,0,h,h,h, 0,h,h,h,h, 0,0,0,h,h, 0,0,h,0,h, 0,h,0,0,h, h,0,0,0,h, 0,0,h,h,0, 0,h,0,h,0, h,0,0,h,0, 0,h,h,0,0, h,0,h,0,0, h,h,0,0,0, h,0,0,0,0, 0,h,0,0,0, 0,0,h,0,0, 0,0,0,h,0, 0,0,0,0,h};
  NT udsign5[21]={-1, 1,1,1,1,1, -1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,1,1,1,1};
  NT udlength5[21]= {-2*h, -h,-h,-h,-h,-h, h,h,h,h,h,h,h,h,h,h, 2*h,2*h,2*h,2*h,2*h};

  if(DIM == 2){udmatrixITfix = udmatrix2; udsignITfix = udsign2;   udlengthITfix = udlength2;
    cubecoordITfix = cubecoord2; cubesignITfix = cubesign2;}
  if(DIM == 3){udmatrixITfix = udmatrix3; udsignITfix = udsign3;   udlengthITfix = udlength3;
    cubecoordITfix = cubecoord3; cubesignITfix = cubesign3;}
  if(DIM == 4){udmatrixITfix = udmatrix4; udsignITfix = udsign4;   udlengthITfix = udlength4;
    cubecoordITfix = cubecoord4; cubesignITfix = cubesign4;}
  if(DIM == 5){udmatrixITfix = udmatrix5; udsignITfix = udsign5;   udlengthITfix = udlength5;
    cubecoordITfix = cubecoord5; cubesignITfix = cubesign5;}

  struct timeval TimeAtStart, TimeNow;
  gettimeofday(&TimeAtStart, NULL);

  //iterations
  for(int n = 0; n <= Ntot-1; ++n){
    // initialize new matrices
    ms = int(pow(simplexfac[DIM],n+1)); //matrixsize
    if(n == Ntot-1){ms = 1;}
    if(coordindicator == 1){
      delete[] list2; delete[] list2length; delete[] list2sign;
      list2=new NT[ms*DIM]; list2length= new NT[ms]; list2sign  = new NT[ms];
      coordIT = list1;    lengthIT   = list1length; signIT     = list1sign;
      NEWcoordIT = list2; NEWlengthIT= list2length; NEWsignIT  = list2sign;
      coordindicator = 2;}
    else{
      delete[] list1; delete[] list1length; delete[] list1sign;
      list1=new NT[ms*DIM]; list1length= new NT[ms]; list1sign  = new NT[ms];
      coordIT    = list2; lengthIT   = list2length; signIT     = list2sign;
      NEWcoordIT = list1; NEWlengthIT= list1length; NEWsignIT  = list1sign;
      coordindicator = 1;}

    // create new coordinates and calculate mass of actual cube
    dp = 0.;
    iterationlength = int(pow(simplexfac[DIM],n));
    for( i = 0; i < iterationlength; ++i){
      xIT = xITfix;
      for( i2 = 0; i2 < DIM; ++i2){*xIT++ = *coordIT++;}
      h  = *lengthIT++; //actual simplexlength
      as = *signIT++; //actual sign
      
      
      if(n != Ntot-1){ //only if not in last iteration
        //creating new coordinates
        udmatrixIT = udmatrixITfix; udsignIT = udsignITfix; 
        udlengthIT = udlengthITfix;
        for( i3 = 0; i3 < simplexfac[DIM]; ++i3){
          xIT           = xITfix;
          for( i2 = 0; i2 < DIM; ++i2){*NEWcoordIT++=*xIT++ +*udmatrixIT++ *h;}
          *NEWlengthIT++= *udlengthIT++ * h;
          *NEWsignIT++  = *udsignIT++ * as;}
      }
      
      // calculate mass of actual cube
      dp += as*JDF(x,h*alpha);
    }// end for i
    Pn[n] = P + dp;
    PnStar[n] = P + expolfac[DIM]*dp;
    
    gettimeofday(&TimeNow, NULL);
    Tn[n] = (TimeNow.tv_sec-TimeAtStart.tv_sec)*1.0 + (TimeNow.tv_usec-TimeAtStart.tv_usec)/1000000.0;

    P += dp;
  }// end for n



  delete[] list2; delete[] list2length; delete[] list2sign;
  delete[] list1; delete[] list1length; delete[] list1sign;


  // time measurement
  Printresults(Ntot, Pn, PnStar, Tn);
  std::cout<< "\n" <<"time needed: "<<Tn[Ntot-1]<<" seconds\n";
  return 0;
}
 

NT JDF(NT* y, NT qh){//computes probability weight of cube Q(b,h)
  //std::cout<< " " << y[0]<< " " << y[1] << " " << y[2]<< " " << qh << "\n" ;
  ccIT = cubecoordITfix; cubesignIT = cubesignITfix;
  a = 0;  qqh = 0.;
  for( i4 = 0; i4 < 2; ++i4){
    //Marginal distributions
    fm1[i4] = fmax(y[0]+ qqh,0)-fmax(y[0]+ qqh -1,0);
    fm2[i4] = fmax(y[1]+ qqh,0)-fmax(y[1]+ qqh -1,0);
    fm3[i4] = fmax(y[2]+ qqh,0)-fmax(y[2]+ qqh -1,0);
    //fm4[i4] = fmax(y[3]+ qqh,0)-fmax(y[3]+ qqh -1,0);
    //fm5[i4] = fmax(y[4]+ qqh,0)-fmax(y[4]+ qqh -1,0);
    qqh = qh;}
  for( i4 = 0; i4 < pow2dim; ++i4){
    a  += (*cubesignIT++)*fm1[ccIT[0]]*fm2[ccIT[1]]*fm3[ccIT[2]];
    //a  += (*cubesignIT++)*fm1[ccIT[0]]*fm2[ccIT[1]]*fm3[ccIT[2]]*fm4[ccIT[3]];
    //a  += (*cubesignIT++)*fm1[ccIT[0]]*fm2[ccIT[1]]*fm3[ccIT[2]]*fm4[ccIT[3]]*fm5[ccIT[4]];
    //cw += (*cubesignIT++)*pow( pow(fm1[ccIT[0]],pt)+pow(fm2[ccIT[1]],pt)+pow(fm3[ccIT[2]],pt)-2,1/pt); // pt <0
    //cw += (*cubesignIT++)*exp(-pow(pow(-log(fm1[ccIT[0]]),pt)+pow(-log(fm2[ccIT[1]]),pt)+pow(-log(fm3[ccIT[2]]),pt),1/pt)); //theta >= 1
    ccIT += DIM;}
  return fabs(a);
  //return a; // even dim
}
 

void Printresults(int N, NT *Pn, NT *PnStar, NT *Tn){
  char bchar = ' '; //character between numbers
  char endchar = ' ';//character at the end
  NT corr = 0.5;  
  std::cout<< "\n n            Pn        PnStar      Time      Error ExtraError\n";
  for(int i = 0; i < N; ++i){
    //cout<<i+1<<" "<<Pn[i]<<" "<<Ln[i]<<" "<<Un[i]<<" "<<Dn[i]<<" "<<Tn[i]<<"\n";
    cout.width(2); cout.fill(' ');
    cout<<i+1<<bchar;
    cout.width(13); cout.precision(9); cout.fill(' ');
    cout<<Pn[i]<<bchar;
    cout.width(13); cout.precision(9); cout.fill(' ');
    cout<<PnStar[i]<<bchar;
    cout.width(9); cout.precision(3); cout.fill(' ');
    cout<<Tn[i]<<bchar;
    cout.width(10); cout.precision(4); cout.fill(' ');
    cout<<Pn[i]-corr<<bchar;
    cout.width(10); cout.precision(4); cout.fill(' ');
    cout<<PnStar[i]-corr<<endchar<<"\n";
  }
} 

