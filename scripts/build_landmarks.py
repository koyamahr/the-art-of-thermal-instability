"""Build a model-specific equilibrium landmark registry from validated outputs and explicitly marked source records."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/'docs'

def load(path):return json.loads((D/path).read_text())
def main():
    ki=load('ki-2002/validation.json');fgh=load('fgh-1969/validation-saturation.json')
    kim=load('kim-kim-2013/source-audit.json')
    sources=['ki-2002/validation.json','fgh-1969/validation-saturation.json','kim-kim-2013/source-audit.json']
    registry=dict(schema_version=1,meaning='Planar static ZP pressure and thermal-equilibrium intersections. Null means unavailable in the stated model/domain; never substitute a pressure extremum.',
        density_warning='KI/Iwasaki n=rho/mH for pure H; FGH prototype n0 counts H+He nuclei and pressure also includes electrons. Kim2013 uses n=rho/(1.37 mH) in its source convention. Do not equate these density columns.',
        reports_sha256={x:hashlib.sha256((D/x).read_bytes()).hexdigest() for x in sources},
        models=[dict(id='ki-2002-corrected',status='computed_corrected_fit',exhibit='ki-2002/#saturation',conductivity=ki['saturation_conductivity'],turning_points=ki['turning_points'],P_sat_over_kB_K_cm3=ki['saturation_equilibria'][0]['P_over_kB_K_cm3'],equilibria=ki['saturation_equilibria']),
                dict(id='iwasaki-2012',status='same_law_and_conventions',same_as='ki-2002-corrected',exhibit='iwasaki-2012/',note='An independent implementation was checked; this is the same equilibrium model, not a new distinct cooling curve.'),
                dict(id='kim-kim-2013',status='source_reported_not_recomputed',exhibit='kim-kim-2013/',P_sat_over_kB_K_cm3=None,source_reported_P_sat_over_kB_K_cm3=kim['model']['source_reported_P_sat_over_kB_K_cm3'],equilibria=[],conductivity='Constant kappa=1e5 erg cm^-1 s^-1 K^-1',density_convention=kim['model']['density_convention'],note='Source value only. Different conductivity from the computed KI/Iwasaki model; eigenfunction reconstruction is planned.'),
                *[dict(id='fgh-prototype-'+x['case'].replace(' ','-').lower(),exhibit='fgh-1969/two-phase.html#saturation',historical_reproduction=False,**x) for x in fgh['results']]],
        future_records='For every added equilibrium model, retain heating/loss normalization, composition and density convention, conductivity, geometry, thermal/chemical closure, domain, P_sat and all equilibrium roots. If unavailable, retain status and reason.')
    (D/'equilibrium-landmarks.json').write_text(json.dumps(registry,indent=2)+'\n')
if __name__=='__main__':main()
