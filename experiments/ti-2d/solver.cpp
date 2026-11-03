// Periodic 2-D Euler + positive thermal conduction + KI radiative source.
// Conserved units: n [cm^-3], velocity [pc/Myr], pressure [m (pc/Myr)^2].
// No density/temperature floors, mass sinks or imposed condensations.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
using A=std::array<double,4>;
constexpr double KB=1.380649e-16,MH=1.6735575e-24,PC=3.0856775814913673e18,MYR=31557600e6;
struct Solver {
 int nx,ny,N,mode; double lx,ly,dx,dy,g=5./3.,mass=1.37*MH,b,heat,kunit,kbase,kcold,tref;
 std::vector<A> U,W,sx,sy,fx,fy,R0,R1,V;
 std::vector<double> temp,kap,tab;
 long long retries=0,fallbacks=0; double table_lo=std::log(2.),table_hi=std::log(30000.);
 Solver(int x,int y,double Lx,double Ly,int m,double kb,double kc,double tr):nx(x),ny(y),N(x*y),mode(m),lx(Lx),ly(Ly),dx(Lx/x),dy(Ly/y),kbase(kb),kcold(kc),tref(tr),U(N),W(N),sx(N),sy(N),fx(N),fy(N),R0(N),R1(N),V(N),temp(N),kap(N),tab(16385) {
  double vu=PC/MYR;b=KB/(mass*vu*vu);heat=2e-26*MYR/(mass*vu*vu);kunit=MYR/(PC*PC*mass*vu*vu);
  for(size_t i=0;i<tab.size();i++){double T=std::exp(table_lo+(table_hi-table_lo)*i/(tab.size()-1));tab[i]=std::log(1e7*std::exp(-118400/(T+1000))+.014*std::sqrt(T)*std::exp(-92/T));}
 }
 int at(int x,int y) const {return ((y+ny)%ny)*nx+(x+nx)%nx;}
 static double mc(double a,double b) {if(a*b<=0)return 0;return std::copysign(std::min({2*std::abs(a),2*std::abs(b),.5*std::abs(a+b)}),a);}
 double ratio(double T) const {if(!(T>=2 && T<=30000))throw std::runtime_error("Temperature outside cooling-table domain, without extrapolation");double z=(std::log(T)-table_lo)*(tab.size()-1)/(table_hi-table_lo);int i=std::min(int(z),int(tab.size())-2);return std::exp(tab[i]+(z-i)*(tab[i+1]-tab[i]));}
 bool prim(const std::vector<A>& a) {
  for(int i=0;i<N;i++){double n=a[i][0];if(!(n>0))return false;double u=a[i][1]/n,v=a[i][2]/n,p=(g-1)*(a[i][3]-.5*n*(u*u+v*v));if(!(p>0) || !std::isfinite(p))return false;
   W[i]={n,u,v,p};temp[i]=p/(b*n);if(!(temp[i]>=2 && temp[i]<=30000))return false;kap[i]=kunit*(kbase+kcold*(tref/temp[i])*(tref/temp[i]));}
  return true;
 }
 A flux(A l,A r,int dir) {
  if(dir==1){std::swap(l[1],l[2]);std::swap(r[1],r[2]);}
  double nl=l[0],nr=r[0],ul=l[1],ur=r[1],vl=l[2],vr=r[2],pl=l[3],pr=r[3];
  double el=pl/(g-1)+.5*nl*(ul*ul+vl*vl),er=pr/(g-1)+.5*nr*(ur*ur+vr*vr);
  A ql={nl,nl*ul,nl*vl,el},qr={nr,nr*ur,nr*vr,er},fl={nl*ul,nl*ul*ul+pl,nl*ul*vl,(el+pl)*ul},fr={nr*ur,nr*ur*ur+pr,nr*ur*vr,(er+pr)*ur};
  double cl=std::sqrt(g*pl/nl),cr=std::sqrt(g*pr/nr),sl=std::min(ul-cl,ur-cr),sr=std::max(ul+cl,ur+cr);A f;
  if(sl>=0)f=fl;else if(sr<=0)f=fr;else {
   double sm=(pr-pl+nl*ul*(sl-ul)-nr*ur*(sr-ur))/(nl*(sl-ul)-nr*(sr-ur));double ps=pl+nl*(sl-ul)*(sm-ul);
   if(!(ps>0) || !std::isfinite(sm)) {fallbacks++;for(int z=0;z<4;z++)f[z]=(sr*fl[z]-sl*fr[z]+sl*sr*(qr[z]-ql[z]))/(sr-sl);}
   else {bool left=sm>=0;double ss=left?sl:sr,uu=left?ul:ur,nn=left?nl:nr,pp=left?pl:pr,vv=left?vl:vr,ee=left?el:er;double ns=nn*(ss-uu)/(ss-sm);
    A star={ns,ns*sm,ns*vv,((ss-uu)*ee-pp*uu+ps*sm)/(ss-sm)};for(int z=0;z<4;z++)f[z]=(left?fl[z]:fr[z])+ss*(star[z]-(left?ql[z]:qr[z]));}
  }
  if(dir==1)std::swap(f[1],f[2]);return f;
 }
 void rhs(const std::vector<A>& a,std::vector<A>& out,double& qsum) {
  if(!prim(a))throw std::runtime_error("Invalid state in RHS");
  bool hydro=(mode==0||mode==1),cond=(mode==0||mode==2),cool=(mode==0||mode==3);
  if(hydro){
   for(int y=0;y<ny;y++)for(int x=0;x<nx;x++){int i=at(x,y);for(int z=0;z<4;z++){sx[i][z]=mc(W[i][z]-W[at(x-1,y)][z],W[at(x+1,y)][z]-W[i][z]);sy[i][z]=mc(W[i][z]-W[at(x,y-1)][z],W[at(x,y+1)][z]-W[i][z]);}}
   for(int y=0;y<ny;y++)for(int x=0;x<nx;x++){int i=at(x,y),ix=at(x+1,y),iy=at(x,y+1);A l,r;for(int z=0;z<4;z++){l[z]=W[i][z]+.5*sx[i][z];r[z]=W[ix][z]-.5*sx[ix][z];}fx[i]=flux(l,r,0);for(int z=0;z<4;z++){l[z]=W[i][z]+.5*sy[i][z];r[z]=W[iy][z]-.5*sy[iy][z];}fy[i]=flux(l,r,1);}
  }
  qsum=0;
  for(int y=0;y<ny;y++)for(int x=0;x<nx;x++){int i=at(x,y),il=at(x-1,y),ib=at(x,y-1),ir=at(x+1,y),it=at(x,y+1);A f={0,0,0,0};
   if(hydro)for(int z=0;z<4;z++)f[z]=-(fx[i][z]-fx[il][z])/dx-(fy[i][z]-fy[ib][z])/dy;
   if(cond){double ki=kap[i];auto face=[&](int j){return 2*ki*kap[j]/(ki+kap[j])*(temp[j]-temp[i]);};f[3]+=(face(ir)+face(il))/(dx*dx)+(face(it)+face(ib))/(dy*dy);}
   if(cool){double q=heat*W[i][0]*(W[i][0]*ratio(temp[i])-1);f[3]-=q;qsum+=q;}
   out[i]=f;
  }
 }
 double timestep() {
  if(!prim(U))throw std::runtime_error("Invalid accepted state");double rate=0;
  for(int i=0;i<N;i++) {
   double n=W[i][0],p=W[i][3],c=std::sqrt(g*p/n),h=0,d=0,s=0;
   if(mode==0||mode==1)h=((std::abs(W[i][1])+c)/dx+(std::abs(W[i][2])+c)/dy)/.35;
   if(mode==0||mode==2){double chi=kap[i]*(g-1)/(b*n);d=chi*(1/(dx*dx)+1/(dy*dy))/.18;}
   if(mode==0||mode==3){double q=std::abs(heat*n*(n*ratio(temp[i])-1));s=q/(.12*p/(g-1));}
   rate=std::max(rate,std::max({h,d,s}));
  }return rate>0?1/rate:1e-3;
 }
 double advance(double dt) {
  double q0,q1;rhs(U,R0,q0);
  for(int attempt=0;attempt<24;attempt++){
   for(int i=0;i<N;i++)for(int z=0;z<4;z++)V[i][z]=U[i][z]+dt*R0[i][z];
   if(!prim(V)){dt*=.5;retries++;continue;}
   rhs(V,R1,q1);
   for(int i=0;i<N;i++)for(int z=0;z<4;z++)V[i][z]=.5*(U[i][z]+V[i][z]+dt*R1[i][z]);
   if(!prim(V)){dt*=.5;retries++;continue;}
   U.swap(V);last_loss=.5*dt*(q0+q1)*dx*dy;return dt;
  }throw std::runtime_error("Positivity retry limit");
 }
 double last_loss=0;
 void frame(const std::string& dir,int n) {char s[100];std::snprintf(s,sizeof(s),"/density-%04d.bin",n);std::ofstream f(dir+s,std::ios::binary);for(auto& u:U){float v=float(u[0]);f.write(reinterpret_cast<char*>(&v),sizeof(v));}}
 void state(const std::string& path) {std::ofstream f(path,std::ios::binary);for(auto& u:U)f.write(reinterpret_cast<char*>(u.data()),4*sizeof(double));}
 void diagnostic(std::ofstream& f,double t,long long steps,double rad) {
  prim(U);double mn=1e99,mx=0,tmin=1e99,tmax=0,field=1e99,p=0,cold=0,mean=0,emax=0;A sum={0,0,0,0};
  for(int i=0;i<N;i++){double n=W[i][0],T=temp[i],pp=W[i][3];mn=std::min(mn,n);mx=std::max(mx,n);tmin=std::min(tmin,T);tmax=std::max(tmax,T);p+=pp/b;mean+=n;if(T<300)cold+=n;
   double k=(kbase+kcold*(tref/T)*(tref/T)),lambda=2e-26*ratio(T);field=std::min(field,std::sqrt(k*T/(n*n*lambda))/PC/std::min(dx,dy));emax=std::max(emax,std::sqrt(W[i][1]*W[i][1]+W[i][2]*W[i][2])/std::sqrt(g*pp/n));
   for(int z=0;z<4;z++)sum[z]+=U[i][z];}
  f<<std::setprecision(16)<<t<<','<<steps<<','<<sum[0]*dx*dy<<','<<sum[1]*dx*dy<<','<<sum[2]*dx*dy<<','<<sum[3]*dx*dy<<','<<rad<<','<<p/N<<','<<mn<<','<<mx<<','<<tmin<<','<<tmax<<','<<cold/mean<<','<<field<<','<<emax<<','<<retries<<','<<fallbacks<<'\n';f.flush();
 }
};
int main(int argc,char**argv){try {
 if(argc!=13)throw std::runtime_error("usage: solver nx ny Lx Ly tend mode kbase kcold Tref input_state output_directory snapshots");
 int nx=std::stoi(argv[1]),ny=std::stoi(argv[2]),mode=std::stoi(argv[6]),frames=std::stoi(argv[12]);double Lx=std::stod(argv[3]),Ly=std::stod(argv[4]),end=std::stod(argv[5]);std::string dir=argv[11];std::filesystem::create_directories(dir);
 Solver s(nx,ny,Lx,Ly,mode,std::stod(argv[7]),std::stod(argv[8]),std::stod(argv[9]));std::ifstream input(argv[10],std::ios::binary);for(auto& u:s.U)input.read(reinterpret_cast<char*>(u.data()),4*sizeof(double));if(!input || !s.prim(s.U))throw std::runtime_error("Invalid input state");
 std::ofstream hist(dir+"/history.csv");hist<<"time_Myr,steps,mass,momentum_x,momentum_y,total_energy,integrated_net_radiative_loss,mean_P_over_kB,min_n,max_n,min_T,max_T,cold_mass_fraction_T_below_300K,min_Field_proxy_cells,max_Mach,retries,HLLE_fallbacks\n";
 double t=0,rad=0;long long steps=0;s.diagnostic(hist,t,steps,rad);s.frame(dir,0);int next=1;double diagnostic_at=end/200.;
 while(t<end){double goal=end*next/frames;double dt=std::min(s.timestep(),goal-t);if(dt<=1e-14){s.frame(dir,next);next++;continue;}
  double used=s.advance(dt);t+=used;rad+=s.last_loss;steps++;
  if(t>=diagnostic_at-1e-13){s.diagnostic(hist,t,steps,rad);diagnostic_at+=end/200.;}
  if(t>=goal-1e-13){s.frame(dir,next);next++;}
 }
 s.diagnostic(hist,t,steps,rad);s.state(dir+"/final.bin");std::cout<<"Completed "<<nx<<'x'<<ny<<" at "<<t<<" Myr, "<<steps<<" steps, "<<s.retries<<" retries\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
