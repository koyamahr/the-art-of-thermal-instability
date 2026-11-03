"""Run the nonlinear prototype, conservation controls, and resolution checks."""
import csv,hashlib,json,platform
import time as timer
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from nonlinear import Rates,initial_state,primitive,hydro,conduction,cooling,total_energy,simulate,ROOT,KB,MH,PC,MYR,recombination


def main():
    c=json.loads((ROOT/'config-nonlinear.json').read_text());linear=json.loads((ROOT/'config.json').read_text());rates=Rates(linear,c);out=ROOT/'output';out.mkdir(exist_ok=True);checks=[];start=timer.perf_counter()
    def check(name,error,tol):checks.append(dict(name=name,error=float(error),tolerance=tol,passed=bool(error<tol)))
    def condition(name,v):checks.append(dict(name=name,passed=bool(v)))
    def profile(run,time):return next(m for t,m in run['snapshots'] if abs(t-time)<1e-9)
    def l1(coarse,fine):
        qc=(np.arange(len(coarse))+.5)/len(coarse);qf=(np.arange(len(fine))+.5)/len(fine)
        interpolated=np.interp(qf,qc,coarse,period=1.)
        return float(np.sum(abs(interpolated-fine))/np.sum(abs(fine)))
    # Hydro energy cancellation is independent of a cooling implementation.
    dx,u,e,x,origin,col=initial_state(128,c);n,t,p,cs=primitive(dx,u,e,x,col,c)
    dd,uu,ee,move=hydro(dx,u,e,x,col,c,.05*min(dx/cs))
    check('An initially isobaric static mesh stays static without sources',np.max(abs(uu))/1e6,1e-12)
    u=1e5*np.sin(2*np.pi*np.arange(len(dx))/len(dx));e=e*(1+.05*np.cos(4*np.pi*np.arange(len(dx))/len(dx)))
    energy=total_energy(u,e,col);dd,uu,ee,move=hydro(dx,u,e,x,col,c,.05*min(dx/cs))
    check('Compatible hydro work conserves thermal plus kinetic energy',abs(total_energy(uu,ee,col)/energy-1),1e-12)
    check('Periodic hydro preserves box length',abs(dd.sum()/dx.sum()-1),1e-13)
    check('Periodic hydro preserves total nodal momentum',abs(np.mean(uu)-np.mean(u))/1e6,1e-13)
    # Linearized periodic heat equation: independent discrete Fourier damping.
    nc=128;length=3*PC;dx=np.full(nc,length/nc);col=.3*length/nc;x=np.full(nc,.05);cv=KB*(1.1+x)/(c['gamma']-1);phase=2*np.pi*(np.arange(nc)+.5)/nc;t0=1000.;amp=1e-5;tt=t0*(1+amp*np.cos(phase));e=cv*tt
    dt=.02*MYR;en=conduction(dx,e,x,col,c,dt);tn=en/cv
    projected=2*np.mean((tn/t0-1)*np.cos(phase))/amp
    eigen=4*c['conductivity_cgs_prefactor']*np.sqrt(t0)*np.sin(np.pi/nc)**2/(col*cv[0]*dx[0]);expected=1/(1+dt*eigen)
    check('Periodic conduction agrees with backward-Euler Fourier damping',abs(projected/expected-1),1e-7)
    check('Periodic conduction conserves summed thermal energy',abs(en.sum()/e.sum()-1),1e-12)
    condition('Implicit conduction respects temperature maximum principle',tn.min()>=tt.min()-1e-8 and tn.max()<=tt.max()+1e-8)
    # Local cooling/recombination versus an independent adaptive ODE integrator.
    init=initial_state(8,dict(c,density_amplitude=0.));dx,u,e,x,origin,col=init
    def source_rhs(time,y):
        temp,ion=np.exp(y);lam,_,_=rates.evaluate(temp,ion);alpha,_=recombination(temp);rr=MYR*.3*ion*alpha
        return [-MYR*.3*(c['gamma']-1)*lam/((1.1+ion)*KB*temp)+(2-c['gamma'])*ion/(1.1+ion)*rr,-rr]
    ref=solve_ivp(source_rhs,[0,.05],np.log([8200.,.05]),method='DOP853',rtol=1e-12,atol=1e-13).y[:,-1]
    for j in range(50):e,x,loss,rec=cooling(dx,e,x,col,c,rates,.001*MYR)
    final=e*(c['gamma']-1)/(KB*(1.1+x))
    check('Nonlinear source temperature agrees with adaptive ODE',float(np.max(abs(final/np.exp(ref[0])-1))),2e-6)
    check('Quadratic hydrogen recombination agrees with adaptive ODE',float(np.max(abs(x/np.exp(ref[1])-1))),2e-6)
    # Small-amplitude nonlinear run must approach the separately derived linear equations.
    from reproduce import rhs
    class Adapter:
        def evaluate(self,t,x):
            lam,lt,lx=rates.evaluate(t,x);return float(lam),float(lt*t/lam),float(lx*x/lam)
    cc=dict(linear,wavelength_pc=[3.]);reflinear=solve_ivp(lambda t,y:rhs(t,y,cc,Adapter()),[0,.4],np.r_[np.log([8200.,.05]),1.,0.,0.,0.],rtol=2e-10,atol=1e-12,method='DOP853').y[2,-1]
    cc=dict(c,snapshot_times_Myr=[0.,.4],time_end_Myr=.4);small=simulate(cc,rates,256,amplitude=1e-5,end_Myr=.4)
    m=profile(small,.4);q=(np.arange(256)+.5)/256;ampl=2*np.mean((m['n']/.3-1)*np.sin(2*np.pi*q))/1e-5
    check('Small-amplitude nonlinear density agrees with linear3pc solution',abs(ampl/reflinear-1),.005)
    runs=[]
    for nc in c['resolutions']:
        print('Running mass-grid cells',nc,flush=True);run=simulate(c,rates,nc);runs.append(run)
        check(str(nc)+': radiative plus recombination energy budget',run['stats']['max_energy_budget_relative_error'],c['convergence_tolerances']['energy_budget_relative_error'])
        check(str(nc)+': periodic box length',run['stats']['max_box_length_relative_error'],1e-11)
        check(str(nc)+': total momentum',run['stats']['max_mean_velocity_over_10kms'],1e-11)
        condition(str(nc)+': nonlinear density contrast exceeds100',max(h['contrast'] for h in run['history'])>100)
        condition(str(nc)+': finite positive primitive states',all(np.all(m['n']>0)&np.all(m['t']>0)&np.all(m['x']>0) for t,m in run['snapshots']))
        condition(str(nc)+': ionization does not exceed initial fraction',all(np.max(m['x'])<=.05*(1+1e-10) for t,m in run['snapshots']))
    print('Running auxiliary grid',c['auxiliary_resolution'],flush=True);aux=simulate(c,rates,c['auxiliary_resolution'])
    check('Auxiliary grid: global energy budget',aux['stats']['max_energy_budget_relative_error'],1e-6)
    check('Auxiliary grid: periodic box length',aux['stats']['max_box_length_relative_error'],1e-11)
    check('Auxiliary grid: total momentum',aux['stats']['max_mean_velocity_over_10kms'],1e-11)
    compare_time=c['comparison_time_Myr'];spatial=[]
    for a,b in zip(runs[:-1],runs[1:]):
        aa=profile(a,compare_time);bb=profile(b,compare_time);spatial.append(dict(coarse_cells=a['stats']['cells'],fine_cells=b['stats']['cells'],density_L1_mass_coordinate=l1(aa['n'],bb['n']),relative_peak_density=float(aa['n'].max()/bb['n'].max()-1),relative_minimum_temperature=float(aa['t'].min()/bb['t'].min()-1)))
    finest=runs[-1]
    resolution_diagnostics=dict(scope='Historical mesh comparison, not a continuum-convergence claim',finest_historical_pair_density_L1=spatial[-1]['density_L1_mass_coordinate'],five_percent_profile_criterion_met=bool(spatial[-1]['density_L1_mass_coordinate']<.05),historical_vs_auxiliary_density_L1=l1(profile(finest,compare_time)['n'],profile(aux,compare_time)['n']),four_Field_cells_criterion_met=bool(aux['stats']['minimum_Field_length_cells']>=4),auxiliary_stats=aux['stats'])
    print('Running half timestep and cold-asymptote controls',flush=True)
    timed=simulate(c,rates,c['time_refinement_cells'],timestep_scale=.5)
    base=next(r for r in runs if r['stats']['cells']==c['time_refinement_cells']);temporal=l1(profile(base,compare_time)['n'],profile(timed,compare_time)['n'])
    check('Half timestep: density-profile L1',temporal,c['convergence_tolerances']['half_step_mass_coordinate_density_L1'])
    alt=simulate(c,Rates(linear,c,cold_electron_power=0.),c['cold_sensitivity_cells'])
    am=profile(alt,compare_time);bm=profile(base,compare_time);cold_probe=dict(scope='Electron T-power0 versus-1/2 below source Tmin; a sensitivity probe, not an uncertainty bound',cells=c['cold_sensitivity_cells'],relative_peak_density=float(am['n'].max()/bm['n'].max()-1),relative_min_temperature=float(am['t'].min()/bm['t'].min()-1),density_profile_L1=l1(am['n'],bm['n']))
    primary_run=next(r for r in runs if r['stats']['cells']==c['primary_cells']);primary=profile(primary_run,compare_time)
    source_anchor=dict(time_Myr=.985,cloud_temperature_K=15,intercloud_temperature_K=1800,cloud_x=.002,cloud_width_pc=.007,cloud_mass_fraction=.2,source='SMS1972 p682 text, approximate cloud quantities, not exact profile samples or equivalently defined width criteria')
    comparison=dict(source=source_anchor,model=dict(peak_nH_cm3=primary['peak_n'],density_contrast=primary['contrast'],min_temperature_K=primary['min_T'],warm_temperature_K=primary['warm_T'],minimum_x=float(primary['x'].min()),cloud_width_pc=primary['cloud_width_pc'],density_FWHM_pc=primary['density_FWHM_pc'],cloud_mass_fraction=primary['cloud_mass_fraction']),note='Measured differences remain. Cloud selection is T<100K and n>10nmean; FWHM uses half the density excess. These need not equal the original cloud-edge definition. No coefficients or time shifts fitted.')
    condition('Condensation reaches contrast100 before the cold extrapolation starts',finest['stats']['first_contrast100_time_bracket_Myr'][1]<finest['stats']['first_cold_extension_time_bracket_Myr'][0])
    if not all(x['passed'] for x in checks):
        (out/'validation-nonlinear-failed.json').write_text(json.dumps(dict(checks=checks,spatial=spatial,temporal=temporal,stats=[r['stats'] for r in runs]),indent=2)+'\n')
        raise RuntimeError(json.dumps([x for x in checks if not x['passed']],indent=2))
    colors=['#777777','#598c82','#c08a45','#685489'];fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained');q=(np.arange(c['primary_cells'])+.5)/c['primary_cells']*100
    for (time,m),color in zip(primary_run['snapshots'][:4],colors):
        for ax,key,scale in zip(axes.flat,['n','u','t','x'],['log','linear','log','log']):
            ax.plot(q,m[key]/1e5 if key=='u' else m[key],color=color,lw=1.4,label=str(time)+' Myr');ax.set_yscale(scale)
    for ax,label in zip(axes.flat,[r'$n_H$ (cm$^{-3}$)','Velocity (km/s)','Temperature (K)','Ionized fraction x']):
        ax.set(xlabel='Equal-H-mass coordinate (%)',ylabel=label);ax.grid(alpha=.15);ax.spines[['top','right']].set_visible(False)
    axes[0,0].legend(frameon=False,fontsize=8);fig.suptitle('Schwarz, McCray & Stein (1972): nonlinear 1D prototype (96 mass cells)')
    fig.supxlabel('Periodic 3 pc, 10% initially isobaric perturbation; compatible moving-mass-grid hydro + conduction + recombination.\nGraph-reconstructed cooling with declared cold/low-ionization extensions. No fitted time shift; historical match remains unvalidated.',fontsize=8.5)
    savefig(fig,out,'nonlinear');plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(14,4.7),layout='constrained')
    for run in runs+[aux]:
        history=run['history'];times=[h['time_Myr'] for h in history]
        axes[0].semilogy(times,[h['contrast'] for h in history],label=str(run['stats']['cells'])+' cells')
        axes[1].semilogy(times,[h['min_T_K'] for h in history]);axes[2].semilogy(times,[h['min_Field_length_cells'] for h in history])
    axes[0].legend(frameon=False,fontsize=8);axes[1].axhline(rates.tmin,color='k',ls='--',lw=1);axes[2].axhline(4,color='k',ls='--',lw=1)
    for ax,label in zip(axes,['Density contrast','Minimum T (K)','Minimum Field-length proxy / cell width']):ax.set(xlabel='Time (Myr)',ylabel=label);ax.grid(alpha=.15)
    fig.suptitle('Historical mesh comparison and source-domain exposure');fig.supxlabel('Field-length proxy: sqrt[K(T)T/(nH² Lambda)]; the four-cell line is a stated diagnostic, not a universal theorem.',fontsize=8.5);savefig(fig,out,'nonlinear-convergence');plt.close(fig)
    with (out/'nonlinear-profiles.csv').open('w',newline='') as f:
        w=csv.writer(f,lineterminator='\n');w.writerow(['time_Myr','cell','mass_coordinate_percent','physical_position_pc','cell_width_pc','nH_cm-3','velocity_kms','temperature_K','x','P_over_k_K_cm-3','Field_proxy_cells'])
        for t,m in primary_run['snapshots']:
            for i in range(len(q)):w.writerow([t,i,q[i],m['pos'][i]/PC,m['dx'][i]/PC,m['n'][i],m['u'][i]/1e5,m['t'][i],m['x'][i],m['p'][i]/KB,m['field_cells'][i]])
    with (out/'nonlinear-history.csv').open('w',newline='') as f:
        fields=['cells']+list(finest['history'][0]);w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader()
        for run in runs+[aux]:
            for h in run['history']:w.writerow(dict(cells=run['stats']['cells'],**h))
    names=['nonlinear.py','reproduce-nonlinear.py','reproduce.py','config.json','cooling-readings.csv','cooling-nonlinear-readings.csv']
    report=dict(scope=c['scope'],primary_cells=c['primary_cells'],checks=checks,runs=[r['stats'] for r in runs],spatial_comparison_at_Myr=compare_time,spatial_comparison=spatial,resolution_diagnostics=resolution_diagnostics,temporal_comparison=dict(cells=c['time_refinement_cells'],density_L1=temporal,half_step_stats=timed['stats']),cold_extension_sensitivity=cold_probe,source_comparison=comparison,limitations=c['limitations'],wall_seconds=timer.perf_counter()-start,environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__),code_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names},config_sha256=hashlib.sha256((ROOT/'config-nonlinear.json').read_bytes()).hexdigest())
    (out/'validation-nonlinear.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(checks_passed=len(checks),resolution_diagnostics=resolution_diagnostics,spatial=spatial,temporal=temporal,finest=finest['stats'],comparison=comparison,wall_seconds=report['wall_seconds']),indent=2))


def savefig(fig,out,name):
    matplotlib.rcParams['svg.hashsalt']='sms72-nonlinear';fig.savefig(out/(name+'.png'),dpi=180);fig.savefig(out/(name+'.svg'),metadata={'Date':None});p=out/(name+'.svg');p.write_text('\n'.join(line.rstrip() for line in p.read_text().splitlines())+'\n')


if __name__=='__main__':main()
