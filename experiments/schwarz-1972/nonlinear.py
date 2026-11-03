"""Conservative1D moving-mass-grid prototype for Schwarz et al.1972.

No fitted time shift, temperature floor, density clipping, or imposed cloud.
See config-nonlinear.json for graph-input extrapolations and solver scope.
"""
from pathlib import Path
import csv,json,time
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.linalg import solve_banded
from reproduce import ROOT,KB,MH,PC,MYR,axis_fit,recombination


class Rates:
    def __init__(self,linear,config,cold_electron_power=-.5):
        d=linear['cooling_digitization'];sx,bx,_=axis_fit(d['x_ticks']);sy,by,_=axis_fit(d['y_ticks'])
        rows=list(csv.DictReader((ROOT/'cooling-nonlinear-readings.csv').open()));self.xs=np.array([.001,.01,.1]);self.curves=[]
        for x in self.xs:
            r=sorted((z for z in rows if float(z['electron_fraction'])==x),key=lambda z:float(z['x_px']))
            ts=10**np.array([sx*float(z['x_px'])+bx for z in r]);ls=10**np.array([sy*float(z['y_px'])+by for z in r])
            self.curves.append(PchipInterpolator(np.log(ts),np.log(ls),extrapolate=False))
        self.tmin=float(min(ts));self.tmax=float(max(ts));self.cold_power=cold_electron_power;self.config=config
    def evaluate(self,t,x):
        t,x=np.broadcast_arrays(np.asarray(t,float),np.asarray(x,float));lt=np.log(np.clip(t,self.tmin,self.tmax))
        vals=np.stack([np.exp(f(lt)) for f in self.curves]);ders=np.stack([vals[i]*f.derivative()(lt)/np.clip(t,self.tmin,self.tmax) for i,f in enumerate(self.curves)])
        i=np.where(x<=.01,0,1);w=(x-self.xs[i])/(self.xs[i+1]-self.xs[i]);idx=np.expand_dims(i,0)
        l0=np.take_along_axis(vals,idx,axis=0)[0];l1=np.take_along_axis(vals,idx+1,axis=0)[0]
        d0=np.take_along_axis(ders,idx,axis=0)[0];d1=np.take_along_axis(ders,idx+1,axis=0)[0]
        slope=(l1-l0)/(self.xs[i+1]-self.xs[i]);value=l0+w*(l1-l0);dt=d0+w*(d1-d0)
        # Low-ionization continuation: neutral plus electron collider terms.
        neutral=np.maximum(vals[0]-.001*(vals[1]-vals[0])/.009,0.)
        electron=neutral+(vals[1]-vals[0])/.009
        lowx=x<.001
        value=np.where(lowx,neutral*(1-x)+electron*x,value)
        dt=np.where(lowx,ders[0]+(x-.001)*(ders[1]-ders[0])/.009,dt)
        slope=np.where(lowx,electron-neutral,slope)
        # Continuous cold asymptote, with the collider mix estimated from the
        # lowest two x curves. This is an exposed model extension, not a floor.
        cold=t<self.tmin
        if np.any(cold):
            if not self.config['cold_extension']['enabled']:raise ValueError('Cold source-data boundary reached')
            b=np.array([float(np.exp(f(np.log(self.tmin)))) for f in self.curves]);s=(b[1]-b[0])/.009;an=max(b[0]-.001*s,0.);ae=an+s
            den=an*(1-x)+ae*x;fe=np.divide(ae*x,den,out=np.zeros_like(x),where=den>0)
            fed=np.divide(ae*an,den**2,out=np.zeros_like(x),where=den>0)
            factor=np.exp(-self.config['cold_extension']['threshold_K']*(1/t-1/self.tmin));h=(t/self.tmin)**self.cold_power;mix=1-fe+fe*h
            dtcold=value*factor*(self.config['cold_extension']['threshold_K']/t**2*mix+fe*h*self.cold_power/t)
            dxcold=factor*(slope*mix+value*fed*(h-1))
            value=np.where(cold,value*factor*mix,value);dt=np.where(cold,dtcold,dt);slope=np.where(cold,dxcold,slope)
        return value,dt,slope


def initial_state(ncells,c,amplitude=None):
    length=c['wavelength_pc']*PC;mean=c['nH_mean_cm3'];amp=c['density_amplitude'] if amplitude is None else amplitude;k=2*np.pi/length
    q=np.linspace(0,length,ncells+1);edge=q.copy()
    for _ in range(10):edge-=(edge+amp/k*(1-np.cos(k*edge))-q)/(1+amp*np.sin(k*edge))
    dx=np.diff(edge);column=mean*length/ncells;n=column/dx;x=np.full(ncells,c['electron_fraction']);u=np.zeros(ncells)
    p0=mean*(1+c['helium_per_H']+c['electron_fraction'])*KB*c['T_background_K']
    energy=p0/(n*(c['gamma']-1))
    return dx,u,energy,x,0.,column


def primitive(dx,u,e,x,column,c):
    n=column/dx;t=e*(c['gamma']-1)/(KB*(1+c['helium_per_H']+x));p=(c['gamma']-1)*n*e
    return n,t,p,np.sqrt(c['gamma']*p/(1.4*MH*n))


def hydro(dx,u,e,x,column,c,dt):
    n,t,p,cs=primitive(dx,u,e,x,column,c);du=np.roll(u,-1)-u;compression=np.maximum(-du,0)
    visc=1.4*MH*n*(c['viscosity_quadratic']*compression**2+c['viscosity_linear']*cs*compression)
    stress=p+visc;acc=-(stress-np.roll(stress,1))/(1.4*MH*column)
    uh=u+.5*dt*acc;dh=dx+.5*dt*du;eh=e-.5*dt*stress*du/column
    if np.any(dh<=0) or np.any(eh<=0):raise ArithmeticError('Invalid midpoint state')
    n,t,p,cs=primitive(dh,uh,eh,x,column,c);duc=np.roll(uh,-1)-uh;compression=np.maximum(-duc,0)
    stress=p+1.4*MH*n*(c['viscosity_quadratic']*compression**2+c['viscosity_linear']*cs*compression)
    un=u-dt*(stress-np.roll(stress,1))/(1.4*MH*column);um=.5*(u+un);change=dt*(np.roll(um,-1)-um)
    dn=dx+change;en=e-stress*change/column
    if np.any(dn<=0) or np.any(en<=0):raise ArithmeticError('Invalid hydro state')
    return dn,un,en,dt*um[0]


def conduction(dx,e,x,column,c,dt):
    if dt==0 or c['conductivity_cgs_prefactor']==0:return e.copy()
    capacity=KB*(1+c['helium_per_H']+x)/(c['gamma']-1);t=e/capacity;k=c['conductivity_cgs_prefactor']*np.sqrt(t)
    kp=np.roll(k,-1);g=2*k*kp/(k+kp)/(.5*(dx+np.roll(dx,-1)))
    right=dt*g/(column*capacity);left=dt*np.roll(g,1)/(column*capacity);diag=1+right+left
    band=np.zeros((3,len(t)));band[1]=diag;band[0,1:]=-right[:-1];band[2,:-1]=-left[1:]
    u1=np.zeros_like(t);u2=np.zeros_like(t);u1[0]=-left[0];u2[-1]=-right[-1]
    sol=solve_banded((1,1),band,np.column_stack([t,u1,u2]),check_finite=False)
    y,z1,z2=sol.T;mat=np.array([[1+z1[-1],z2[-1]],[z1[0],1+z2[0]]]);weights=np.linalg.solve(mat,[y[-1],y[0]])
    tn=y-z1*weights[0]-z2*weights[1]
    if np.any(tn<=0):raise ArithmeticError('Nonpositive conductive state')
    return tn*capacity


def cooling(dx,e,x,column,c,rates,dt):
    n=column/dx;h=1+c['helium_per_H'];g=c['gamma'];old=e*(g-1)/(KB*(h+x));lo=.001*old;hi=old*(h+x)/h;t=old.copy()
    for _ in range(16):
        mid=.5*(old+t);alpha,sa=recombination(mid);a=n*x*alpha*dt;xn=x/(1+a);xp=-xn*a/(1+a)*sa/(2*mid);xm=.5*(x+xn)
        lam,lt,lx=rates.evaluate(mid,xm);loss=dt*n*lam;rec=KB*mid*(x-xn)
        residual=KB*t*(h+xn)/(g-1)+loss+rec-e
        deriv=KB*(h+xn+t*xp)/(g-1)+dt*n*(.5*lt+.5*lx*xp)+KB*(.5*(x-xn)-mid*xp)
        error=np.max(abs(residual)/e)
        if error<1e-11:return KB*t*(h+xn)/(g-1),xn,loss,rec
        lo=np.where(residual<0,t,lo);hi=np.where(residual>=0,t,hi);trial=t-residual/deriv
        t=np.where((trial>lo)&(trial<hi)&np.isfinite(trial),trial,.5*(lo+hi))
    raise ArithmeticError('Cooling solve did not converge')


def total_energy(u,e,column):return column*np.sum(e+.5*1.4*MH*u*u)


def measure(state,c,rates):
    dx,u,e,x,origin,column=state;n,t,p,cs=primitive(dx,u,e,x,column,c);edge=origin+np.r_[0,np.cumsum(dx)];pos=.5*(edge[:-1]+edge[1:])
    lam,_,_=rates.evaluate(t,x);field=np.sqrt(c['conductivity_cgs_prefactor']*np.sqrt(t)*t/(n*n*lam))
    cloud=(t<100)&(n>10*c['nH_mean_cm3'])
    return dict(n=n,t=t,p=p,u=.5*(u+np.roll(u,-1)),x=x,pos=pos,dx=dx,field_cells=field/dx,
                contrast=float(n.max()/n.min()),peak_n=float(n.max()),min_T=float(t.min()),warm_T=float(t[np.argmin(n)]),cloud_width_pc=float(np.sum(dx[cloud])/PC),density_FWHM_pc=float(np.sum(dx[n>=n.min()+.5*(n.max()-n.min())])/PC),cloud_mass_fraction=float(np.mean(cloud)),cold_extension_mass_fraction=float(np.mean(t<rates.tmin)),low_x_mass_fraction=float(np.mean(x<rates.xs[0])))


def simulate(c,rates,ncells,timestep_scale=1.,amplitude=None,end_Myr=None,conduction_on=True):
    c=dict(c)
    if not conduction_on:c['conductivity_cgs_prefactor']=0.
    state=initial_state(ncells,c,amplitude);dx,u,e,x,origin,column=state;initial_energy=total_energy(u,e,column);length=dx.sum();losses=0.;rec_losses=0.;time_s=0.;steps=0;retries=0
    targets=np.unique(np.r_[c['snapshot_times_Myr'],end_Myr if end_Myr is not None else c['time_end_Myr']]);end=c['time_end_Myr'] if end_Myr is None else end_Myr;targets=targets[targets<=end+1e-12];snapshots=[];history=[];max_energy_error=0.;max_length_error=0.;max_momentum=0.;ever_cold=np.zeros(ncells,bool);ever_lowx=np.zeros(ncells,bool);min_field=np.inf;first_cold=None;first_lowx=None;contrast100=None
    def record():
        nonlocal min_field,max_energy_error,max_length_error,max_momentum
        m=measure((dx,u,e,x,origin,column),c,rates);ever_cold[:]|=m['t']<rates.tmin;ever_lowx[:]|=x<rates.xs[0];min_field=min(min_field,float(m['field_cells'].min()))
        enerr=abs((total_energy(u,e,column)+losses+rec_losses)/initial_energy-1);max_energy_error=max(max_energy_error,enerr);max_length_error=max(max_length_error,abs(dx.sum()/length-1));max_momentum=max(max_momentum,abs(np.mean(u))/1e6)
        history.append(dict(time_Myr=time_s/MYR,contrast=m['contrast'],peak_n_cm3=m['peak_n'],min_T_K=m['min_T'],warm_T_K=m['warm_T'],cloud_width_pc=m['cloud_width_pc'],density_FWHM_pc=m['density_FWHM_pc'],cloud_mass_fraction=m['cloud_mass_fraction'],cold_extension_mass_fraction=m['cold_extension_mass_fraction'],low_x_mass_fraction=m['low_x_mass_fraction'],energy_budget_error=enerr,min_Field_length_cells=float(m['field_cells'].min())))
        return m
    record();snapshots.append((0.,measure(state,c,rates)));index=1;wall=time.perf_counter()
    while time_s<end*MYR*(1-1e-14):
        n,t,p,cs=primitive(dx,u,e,x,column,c);lam,_,_=rates.evaluate(t,x);alpha,_=recombination(t)
        sound_dt=c['cfl']*np.min(dx/(cs+abs(np.roll(u,-1)-u)))
        source_dt=c['cooling_fraction']*np.min(e/(n*lam+n*x*x*alpha*KB*t))
        target=targets[index] if index<len(targets) else end;min_field=min(min_field,float(np.min(np.sqrt(c['conductivity_cgs_prefactor']*np.sqrt(t)*t/(n*n*lam))/dx)))
        dt=min(timestep_scale*sound_dt,timestep_scale*source_dt,target*MYR-time_s)
        if dt<=0:raise ArithmeticError('Nonpositive timestep')
        for retry in range(20):
            try:
                ee=conduction(dx,e,x,column,c,.5*dt);ee,xx,l1,r1=cooling(dx,ee,x,column,c,rates,.5*dt)
                dd,uu,ee,move=hydro(dx,u,ee,xx,column,c,dt);ee,xx,l2,r2=cooling(dd,ee,xx,column,c,rates,.5*dt);ee=conduction(dd,ee,xx,column,c,.5*dt)
                tt=ee*(c['gamma']-1)/(KB*(1+c['helium_per_H']+xx))
                if np.any(tt>rates.tmax*1.00001):raise ValueError('Accepted state exceeds warm cooling input')
                if np.any(tt<1.):raise ValueError('Accepted state below1K validity guard')
                break
            except ArithmeticError:dt*=.5;retries+=1
        else:raise ArithmeticError('Step rejection budget exhausted')
        dx,u,e,x=dd,uu,ee,xx;origin+=move;losses+=column*np.sum(l1+l2);rec_losses+=column*np.sum(r1+r2);time_s+=dt;steps+=1
        ever_cold[:]|=tt<rates.tmin;ever_lowx[:]|=x<rates.xs[0]
        bracket=[float((time_s-dt)/MYR),float(time_s/MYR)]
        if first_cold is None and np.any(tt<rates.tmin):first_cold=bracket
        if first_lowx is None and np.any(x<rates.xs[0]):first_lowx=bracket
        if contrast100 is None and np.max(dx)/np.min(dx)>=100:contrast100=bracket
        max_energy_error=max(max_energy_error,abs((total_energy(u,e,column)+losses+rec_losses)/initial_energy-1));max_length_error=max(max_length_error,abs(dx.sum()/length-1));max_momentum=max(max_momentum,abs(np.mean(u))/1e6)
        if steps%c['history_interval_steps']==0:record()
        if abs(time_s/MYR-target)<1e-10:
            m=record();snapshots.append((float(target),m));index+=1
        if steps>500000:raise ArithmeticError('Timestep budget exceeded')
    stats=dict(first_contrast100_time_bracket_Myr=contrast100,first_cold_extension_time_bracket_Myr=first_cold,first_low_x_extension_time_bracket_Myr=first_lowx,cells=ncells,timestep_scale=timestep_scale,steps=steps,rejected_steps=retries,wall_seconds=time.perf_counter()-wall,max_energy_budget_relative_error=max_energy_error,max_box_length_relative_error=max_length_error,max_mean_velocity_over_10kms=max_momentum,minimum_Field_length_cells=min_field,mass_fraction_ever_below_input_T=float(np.mean(ever_cold)),mass_fraction_ever_below_input_x=float(np.mean(ever_lowx)),cooling_T_domain_K=[rates.tmin,rates.tmax],initial_energy_per_area=initial_energy,radiated_energy_per_area=float(losses),recombination_energy_per_area=float(rec_losses),end_time_Myr=float(time_s/MYR))
    return dict(snapshots=snapshots,history=history,stats=stats)


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--cells',type=int,default=64);p.add_argument('--end',type=float,default=1.05);p.add_argument('--scale',type=float,default=1.);args=p.parse_args()
    c=json.loads((ROOT/'config-nonlinear.json').read_text());r=Rates(json.loads((ROOT/'config.json').read_text()),c);result=simulate(c,r,args.cells,args.scale,end_Myr=args.end)
    print(json.dumps(result['stats'],indent=2));print(json.dumps([{k:v for k,v in h.items() if k in ['time_Myr','contrast','peak_n_cm3','min_T_K','warm_T_K','cloud_width_pc','cloud_mass_fraction']} for h in result['history'][::max(1,len(result['history'])//12)]],indent=2))
