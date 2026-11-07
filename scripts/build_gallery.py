"""Curate archival movies; regenerate with --source-root, or check committed media."""
from pathlib import Path
import argparse,hashlib,json,re,subprocess,shutil
from html import escape
ROOT=Path(__file__).resolve().parents[1];DEST=ROOT/'docs/gallery';CONFIG=ROOT/'media/gallery/selection.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify():
 m=json.loads((DEST/'manifest.json').read_text());assert len(m['clips'])==8
 assert sha(CONFIG)==m['selection_sha256']
 assert sha(Path(__file__))==m['builder_sha256']
 assert len(set(c['source_sha256'] for c in m['clips']))==8
 for name,h in m['assets_sha256'].items():assert sha(DEST/name)==h,name
 for c in m['clips']:assert c['source_frames']==c['output_frames'] and c['output_dimensions']==[c['width'],c['height']]
 assert (DEST/'selection.json').read_bytes()==CONFIG.read_bytes()
 print('Gallery: eight byte-distinct sources, native dimensions/frame counts and committed media hashes verified. This is a technical media check, not scientific reproduction validation.')
def run(ff,args):
 r=subprocess.run([ff,'-hide_banner','-v','error','-xerror',*args],capture_output=True,text=True,errors='replace',check=True)
 return r

def encode(source):
 import imageio_ffmpeg
 ff=imageio_ffmpeg.get_ffmpeg_exe();selection=json.loads(CONFIG.read_text());clips=[]
 DEST.mkdir(parents=True,exist_ok=True)
 for c0 in selection['clips']:
  c=dict(c0);p=source/c['source_relative_path'];assert sha(p)==c['source_sha256'],str(p)
  c['video']=c['id']+'.mp4';c['poster']=c['id']+'.png'
  r=run(ff,['-i',str(p),'-map','0:v:0','-an','-progress','pipe:1','-f','null','-']);c['source_frames']=int(re.findall(r'^frame=(\d+)',r.stdout,re.M)[-1])
  run(ff,['-y','-i',str(p),'-map','0:v:0','-an','-fps_mode','passthrough','-map_metadata','-1','-c:v','libx264','-preset','slow','-crf','16','-pix_fmt','yuv420p','-movflags','+faststart',str(DEST/c['video'])])
  r=run(ff,['-i',str(DEST/c['video']),'-map','0:v:0','-an','-progress','pipe:1','-f','null','-']);c['output_frames']=int(re.findall(r'^frame=(\d+)',r.stdout,re.M)[-1]);assert c['source_frames']==c['output_frames']
  r=subprocess.run([ff,'-hide_banner','-i',str(DEST/c['video'])],capture_output=True,text=True,errors='replace');line=next(s for s in r.stderr.splitlines() if 'Video:' in s);wh=re.search(r'(\d{2,5})x(\d{2,5})',line);c['output_dimensions']=[int(wh[1]),int(wh[2])];assert c['output_dimensions']==[c['width'],c['height']]
  c['poster_time_s']=round(c['duration_s']*c['poster_fraction'],6)
  run(ff,['-y','-ss',str(c['poster_time_s']),'-i',str(p),'-frames:v','1',str(DEST/c['poster'])]);c['output_bytes']=(DEST/c['video']).stat().st_size
  clips.append(c);print(c['id'],c['source_frames'],'frames,',c['output_bytes'],'bytes',flush=True)
 reel=showreel(ff,source,clips)
 shutil.copy(CONFIG,DEST/'selection.json')
 shutil.copy(ROOT/'media/gallery/README.md',DEST/'README.md')
 assets=[DEST/c[key] for c in clips for key in ['video','poster']]+[DEST/'showreel.mp4',DEST/'showreel.png',DEST/'selection.json',DEST/'README.md']
 manifest=dict(review_date=selection['review_date'],scope=selection['scope'],editing=selection['editing'],provenance_limit=selection['provenance_limit'],selection_sha256=sha(CONFIG),builder_sha256=sha(Path(__file__)),ffmpeg_version=subprocess.run([ff,'-version'],capture_output=True,text=True).stdout.splitlines()[0],clips=clips,showreel=reel,assets_sha256={p.name:sha(p) for p in assets})
 (DEST/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n');page(manifest);verify()

def showreel(ff,source,clips):
 import imageio_ffmpeg
 from PIL import Image,ImageDraw,ImageFont
 import matplotlib,numpy as np
 fontpath=Path(matplotlib.get_data_path())/'fonts/ttf/DejaVuSans.ttf'
 def font(size):return ImageFont.truetype(str(fontpath),size)
 writer=imageio_ffmpeg.write_frames(str(DEST/'showreel.mp4'),(800,700),fps=15,codec='libx264',pix_fmt_out='yuv420p',quality=8,macro_block_size=1,output_params=['-movflags','+faststart']);writer.send(None)
 def card(title,subtitle):
  im=Image.new('RGB',(800,700),'#080f1b');d=ImageDraw.Draw(im);d.text((44,270),title,font=font(32),fill='#f0eddc');d.text((44,325),subtitle,font=font(17),fill='#bdcbbb');return im
 opening=card('An archive in motion','Four edited excerpts from historical simulation recordings')
 for _ in range(15):writer.send(np.asarray(opening).tobytes())
 segments=[]
 starts=[2.,5.,3.,12.]
 for index,(c,start) in enumerate(zip(clips[:4],starts)):
  top=136 if c['id']=='ti-mhd-3d' else (101 if c['id']=='ti-conduction-viscosity-3d' else 103)
  frames=imageio_ffmpeg.read_frames(str(source/c['source_relative_path']),pix_fmt='rgb24',input_params=['-ss',str(start)],output_params=['-t','4','-vf','fps=15']);meta=next(frames);count=0
  for raw in frames:
   frame=Image.frombytes('RGB',tuple(meta['size']),raw);cropped=frame.crop((0,top,frame.width,frame.height));cropped.thumbnail((700,568))
   im=Image.new('RGB',(800,700),'#000000');im.paste(cropped,((800-cropped.width)//2,82+(568-cropped.height)//2));d=ImageDraw.Draw(im)
   d.text((28,18),c['title'],font=font(25),fill='#f0eddc');d.text((28,52),c['category']+'  |  edited excerpt',font=font(15),fill='#b9c9c0')
   d.text((28,657),'Source frame credit: Koyama & Inutsuka 2002',font=font(14),fill='#bac6bd');d.text((736,22),f'{index+1}/4',font=font(17),fill='#bac6bd')
   writer.send(np.asarray(im).tobytes());count+=1
  assert count==60,(c['id'],count)
  segments.append(dict(id=c['id'],source_start_s=start,source_end_s=start+4,source_crop_xywh=[0,top,c['width'],c['height']-top],output_frames=count))
 ending=card('Explore the original recordings','Full-frame movies, native resolution, documented provenance')
 for _ in range(15):writer.send(np.asarray(ending).tobytes())
 writer.close();run(ff,['-y','-ss','6.5','-i',str(DEST/'showreel.mp4'),'-frames:v','1',str(DEST/'showreel.png')]);run(ff,['-i',str(DEST/'showreel.mp4'),'-f','null','-'])
 return dict(video='showreel.mp4',poster='showreel.png',width=800,height=700,fps=15,duration_s=18,segments=segments,editing='Hard cuts between four distinct clips; four seconds from each source, with original top title frames cropped and new English captions added. Intro/outro one second each. Source imagery is fitted without upscaling, color remapping or generated detail. Frames repeat to encode at 15 fps; native source motion is not interpolated. No continuous physical chronology is implied across cuts.')

def page(m):
 cards=[]
 for c in m['clips']:
  esc=lambda key:escape(str(c[key]))
  cards.append(f'''<article class="movie-card" id="{c['id']}"><div class="movie-heading"><p class="eyebrow">{esc('category')}</p><h3>{esc('title')}</h3></div><div class="movie-media"><video controls playsinline preload="none" poster="{c['poster']}" width="{c['width']}" height="{c['height']}" aria-label="{esc('title')}"><source src="{c['video']}" type="video/mp4"><a href="{c['video']}">Download MP4</a></video></div><div class="movie-copy"><p>{esc('description')}</p><p class="movie-spec">{c['width']} × {c['height']} px · {c['fps']:g} fps · {c['duration_s']:.1f} s</p><p class="quality"><strong>Viewing quality:</strong> {esc('quality')}</p><details><summary>Source and editing record</summary><p>{esc('caption_evidence')}</p><p>Original: <code>{escape(c['source_relative_path'])}</code>. Complete original frame and timing retained in this version; H.264/yuv420p conversion only. No audio was present. Exact paper figure and run parameters remain unverified.</p></details><p class="movie-links"><a href="{c['video']}">MP4</a> · <a href="{c['poster']}">Poster / PNG</a></p></div></article>''')
 html='''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Historical simulation gallery · The Art of Thermal Instability</title><meta name="description" content="Eight recovered astronomical simulation movies: thermal instability, two-phase gas, shocks and MHD, with an edited showreel, native-resolution versions and provenance."><link rel="stylesheet" href="../styles.css?v=20261108-gallery"><link rel="stylesheet" href="../field-1965/exhibit.css"><link rel="stylesheet" href="gallery.css"></head><body><header><a class="brand" href="../">The Art of Thermal Instability</a><nav aria-label="Gallery navigation"><a href="#showreel">Edited selection</a><a href="#movies">Full-frame movies</a><a href="#provenance">Provenance</a></nav></header><main class="exhibit"><section class="hero"><p class="eyebrow">From the original movie archive / recovered November 2026</p><h1>An archive in motion.</h1><p class="intro">Clouds emerge, shocks gather gas, and magnetized structures evolve. Eight selected recordings bring the historical calculations back to the browser.</p><p class="notice">Archival footage, selected for visual demonstration. These movies have not been matched to recovered run inputs or independently validated against a paper. Original frame captions are retained in the full-frame versions.</p></section><section id="showreel"><p class="eyebrow">An edited introduction</p><h2>Four views into the original simulations.</h2><video class="showreel" controls playsinline preload="none" width="800" height="700" poster="showreel.png" aria-label="Edited introduction to four archival 3D simulation clips"><source src="showreel.mp4" type="video/mp4"><a href="showreel.mp4">Download the edited selection</a></video><p class="small">18 seconds, four distinct 3D excerpts. New English headings replace the old top title frames; the 3D scenes and their source cadence are preserved. Cuts separate different calculations. Repeated frames encode the low-cadence sources at 15 fps; no intermediate simulation states have been invented. <a href="showreel.mp4">Download MP4</a>.</p></section><section id="movies"><p class="eyebrow">The full-frame collection</p><h2>Keep the original picture in view.</h2><p>Play one recording at a time. Each retains its original aspect ratio, annotations and temporal sequence. Small lettering is best inspected at native size using the MP4 link. Downloading or enlarging a recording does not recover detail absent from the source.</p><div class="movie-grid">'''+''.join(cards)+'''</div></section><section id="provenance"><h2>What these recordings can tell us.</h2><p>They are effective visual demonstrations at compact gallery sizes. Their native resolutions are 500–800 pixels wide and their cadences are 3–10 fps. Fine textures can look grainy, and low-frame-rate motion is visibly stepped. They are suitable as historical exhibits; they are not high-resolution, full-screen source material.</p><p>The 3D recordings carry the embedded credit “Koyama &amp; Inutsuka 2002”. This identifies the caption on the recording; it does not establish a match to a specific paper figure, simulation dimensionality in that paper, or verified numerical parameters. The 2D titles here describe the visible panels and evolution. A number in a filename is not treated as a confirmed computational mesh size.</p><p>The gallery was curated from 87 video files, representing 55 byte-distinct recordings after SHA-256 deduplication. Astronomy candidates were visually inspected at several times. The eight selected originals and their MP4 derivatives were decoded completely; dimensions and frame counts were checked. These are technical media checks, not numerical or physical validation of the simulations.</p><p>Materials molecular-dynamics movies and related materials CG are excluded, along with unrelated footage and journal-supplement movies from other collections. Source originals remain untouched. Full-frame versions use H.264/yuv420p, retain the original dimensions and cadence, and add no spatial or temporal interpolation. The showreel has its own explicit crop and excerpt log.</p><div class="downloads"><a href="manifest.json">Sources, edit log &amp; media hashes / JSON</a><a href="selection.json">Selection &amp; quality notes / JSON</a><a href="README.md">Editing and regeneration guide</a><a href="https://github.com/koyamahr/the-art-of-thermal-instability/blob/main/scripts/build_gallery.py">Conversion and gallery builder</a><a href="../ti-2d/">The separate reproducible cover experiment</a></div></section></main><footer><a href="../#atlas">Back to the paper atlas</a><a href="../progress.html">Scientific progress</a></footer><script src="gallery.js" defer></script></body></html>'''
 (DEST/'index.html').write_text(html+'\n')

if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--source-root',type=Path);ap.add_argument('--verify-only',action='store_true');a=ap.parse_args()
 if a.verify_only:verify()
 elif a.source_root:encode(a.source_root)
 else:ap.error('Choose --source-root for local regeneration, or --verify-only for committed media.')
