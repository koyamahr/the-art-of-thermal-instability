"""Turn actual solver frames into fixed-scale media; archive data and diagnostics."""
from model import *
import hashlib,shutil,sys
from PIL import Image,ImageDraw,ImageFont
import imageio_ffmpeg
from matplotlib.colors import LinearSegmentedColormap

COLORS=['#080f1b','#173d44','#4f8078','#bdc89f','#f2dfad']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    c=json.loads((HERE/'config.json').read_text());raw=HERE/'output/raw';out=HERE/'output';nx,ny=c['grid'];v=c['video'];w,hh=v['width'],v['height']
    h=np.genfromtxt(raw/'history.csv',delimiter=',',names=True);coarse=np.genfromtxt('/tmp/ti2d-coarse/history.csv',delimiter=',',names=True)
    u=np.fromfile(raw/'final.bin',dtype='<f8').reshape(ny,nx,4);n,t,vx,vy,p=primitive(u)
    np.savez_compressed(out/'states.npz',initial_conserved=initial(c),final_conserved=u,x_pc=(np.arange(nx)+.5)*c['domain_pc'][0]/nx,y_pc=(np.arange(ny)+.5)*c['domain_pc'][1]/ny)
    shutil.copy(raw/'history.csv',out/'history.csv');shutil.copy('/tmp/ti2d-coarse/history.csv',out/'history-coarse.csv')
    cmap=LinearSegmentedColormap.from_list('thermal',COLORS)
    paths=sorted(raw.glob('density-*.bin'));assert len(paths)==c['snapshots']+1
    writer=imageio_ffmpeg.write_frames(str(out/'condensation.mp4'),(w,hh),fps=v['fps'],codec='libx264',pix_fmt_in='rgb24',pix_fmt_out='yuv420p',quality=8,macro_block_size=16,output_params=['-movflags','+faststart'])
    writer.send(None)
    try:font=ImageFont.truetype('/System/Library/Fonts/Helvetica.ttc',22)
    except OSError:font=ImageFont.load_default()
    poster_index=round(.78*(len(paths)-1))
    for i,path in enumerate(paths):
        d=np.fromfile(path,dtype='<f4').reshape(ny,nx)
        z=(np.log10(d)-np.log10(v['density_min_cm3']))/np.log10(v['density_max_cm3']/v['density_min_cm3'])
        im=Image.fromarray((cmap(np.clip(z,0,1))[...,:3]*255).astype('uint8')[::-1]).resize((w,hh),Image.Resampling.BICUBIC)
        draw=ImageDraw.Draw(im);draw.rectangle((0,hh-55,w,hh),fill=(8,15,27));draw.text((24,hh-41),f"{i*c['end_time_Myr']/c['snapshots']:.2f} Myr   |   log density   |   4.8 x 3.0 pc   |   enhanced conduction",font=font,fill='#e8ecd9')
        writer.send(np.asarray(im).tobytes())
        if i==poster_index:im.save(out/'poster.png')
    writer.close()
    scale=np.linspace(0,1,600)[None,:];Image.fromarray((cmap(np.repeat(scale,24,axis=0))[...,:3]*255).astype('uint8')).save(out/'density-scale.png')
    stats=dict(grid=c['grid'],dx_pc=c['domain_pc'][0]/nx,end_Myr=float(h[-1]['time_Myr']),frames=len(paths),poster_time_Myr=poster_index*c['end_time_Myr']/c['snapshots'],min_n_cm3=float(h['min_n'].min()),max_n_cm3=float(h['max_n'].max()),min_T_K=float(h['min_T'].min()),max_T_K=float(h['max_T'].max()),min_Field_proxy_cells_sampled=float(h['min_Field_proxy_cells'].min()),max_Mach_sampled=float(h['max_Mach'].max()),final_cold_mass_fraction=float(h[-1]['cold_mass_fraction_T_below_300K']),coarse_final_cold_mass_fraction=float(coarse[-1]['cold_mass_fraction_T_below_300K']),final_mean_P_over_kB=float(h[-1]['mean_P_over_kB']),coarse_final_mean_P_over_kB=float(coarse[-1]['mean_P_over_kB']),relative_mass_drift=float(np.max(abs(h['mass']/h[0]['mass']-1))),relative_energy_ledger_residual=float(np.max(abs(h['total_energy']-h[0]['total_energy']+h['integrated_net_radiative_loss']))/h[0]['total_energy']),retries=int(h[-1]['retries']),HLLE_fallbacks=int(h[-1]['HLLE_fallbacks']),model_landmarks=landmarks(c))
    (out/'diagnostics.json').write_text(json.dumps(stats,indent=2)+'\n')
    report=json.loads((HERE/'validation.json').read_text())
    for name,val,limit in [('movie run mass conservation',stats['relative_mass_drift'],1e-10),('movie run radiative energy ledger',stats['relative_energy_ledger_residual'],1e-9),('positive temperatures',0 if stats['min_T_K']>0 else 1,.5),('movie has all actual frames',0 if len(paths)==c['snapshots']+1 else 1,.5)]:
        report['checks'].append(dict(name=name,value=val,limit=limit,passed=val<limit))
    report['code_sha256']={p.name:sha(p) for p in HERE.iterdir() if p.suffix in ('.py','.cpp') and p.name!='publish.py'}
    report['config_sha256']=sha(HERE/'config.json');(HERE/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(stats,indent=2))
if __name__=='__main__':main()
