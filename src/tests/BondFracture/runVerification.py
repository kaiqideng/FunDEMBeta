#!/usr/bin/env python3
"""Run native strength-initiated Euclidean-separation BK cases and retain scoped verification evidence."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import shutil
import subprocess
import time

CASES=("tension","shear","bending","torsion","bending-zero-torque")
PARITY_TOL=2e-9
REFINEMENT_TOL=5e-3
WORK_ABS_TOL=2e-7
GC=.005
INCREMENTS=(("base",1),("half",2),("quarter",4))
FIELDS=("damage","onset_effective_stress","onset_effective_displacement","maximum_effective_displacement",
        "failure_effective_displacement","effective_fracture_energy","normal_energy","shear_energy","bending_energy","torsional_energy",
        "force_x","force_y","force_z","torque_x","torque_y","torque_z","master_torque_x","master_torque_y","master_torque_z",
        "slave_torque_x","slave_torque_y","slave_torque_z","cumulative_external_work","cumulative_effective_work","cumulative_damage_dissipation")

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def states(path):
    with path.open(newline="") as stream:return [{key:float(value) for key,value in row.items()} for row in csv.DictReader(stream)]
def write_csv(path,rows):
    with path.open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)
def source_files(source):
    relatives=["src/interaction/bond.h","src/interaction/bond.cpp","src/execution/bondFunctions.h","src/execution/cuda/bondForceKernel.cu",
               "src/execution/cuda/bondForceKernel.cuh","src/tests/BondFractureTest.cpp","src/tests/PhysicsTest.cpp","src/tests/MaterialPairTest.cpp","src/tests/CMakeLists.txt","CMakeLists.txt"]
    return relatives+sorted(str(path.relative_to(source)) for path in (source/"src/tests/BondFracture").rglob("*")
                            if path.is_file() and "__pycache__" not in path.parts)
def source_identity(source,relatives):
    hashes={relative:sha(source/relative) for relative in relatives}
    def git(*arguments):
        try:
            result=subprocess.run(["git","-C",str(source),*arguments],text=True,capture_output=True,check=True)
            return result.stdout.strip()
        except (OSError,subprocess.CalledProcessError):return None
    status=git("status","--porcelain","--untracked-files=no","--",".")
    return {"git_head":git("rev-parse","HEAD"),"tracked_tree_dirty":None if status is None else bool(status),
            "tracked_status":status,"files":hashes,
            "source_sha256":hashlib.sha256(json.dumps(hashes,sort_keys=True).encode()).hexdigest()}
def compare(first,second,stride=1):
    if len(second)!=stride*(len(first)-1)+1:raise ValueError("The refinement histories do not share matching prescribed poses.")
    errors={}
    for field in FIELDS:
        scale=max(max(abs(row[field]) for row in first),max(abs(row[field]) for row in second),1e-8)
        absolute=max(abs(row[field]-second[stride*index][field]) for index,row in enumerate(first))
        errors[field]={"absolute":absolute,"normalized":absolute/scale}
    return errors
def independent_measures(row,case):
    normal=[row[f"master_{axis}"]-row[f"slave_{axis}"] for axis in ("x","y","z")]
    length=math.sqrt(sum(value*value for value in normal));normal=[value/length for value in normal]
    force=[row[f"force_{axis}"] for axis in ("x","y","z")]
    torque=[row[f"torque_{axis}"] for axis in ("x","y","z")]
    projection=sum(a*b for a,b in zip(force,normal));twist=sum(a*b for a,b in zip(torque,normal))
    shear=math.sqrt(sum((value-projection*n)**2 for value,n in zip(force,normal)))
    bending=math.sqrt(sum((value-twist*n)**2 for value,n in zip(torque,normal)))
    radius=1./math.sqrt(math.pi);inertia=radius*radius/4.
    axial=-projection
    bending_stress=bending*radius/inertia
    sigma=axial+bending_stress
    tau=shear+abs(twist)*radius/(2.*inertia)
    opening_stress=max(axial,0.)+bending_stress
    stress=math.hypot(opening_stress,tau)
    opening=max(row["normal_energy"],0.) if axial>0. else 0.
    normal_energy=opening+max(row["bending_energy"],0.)
    shear_energy=max(row["shear_energy"],0.)+max(row["torsional_energy"],0.)
    dn=2.*normal_energy/opening_stress if opening_stress>0. else 0.
    ds=2.*shear_energy/tau if tau>0. else 0.
    displacement=math.hypot(dn,ds)
    return max(sigma,0.)**2+tau*tau,stress,displacement

def history_audit(rows,case):
    index=next((i for i,row in enumerate(rows) if row["damage_started"]!=0.),None)
    if index is None:return {"started":False,"pass":True}
    first=rows[index]
    criterion,stress,separation=independent_measures(first,case)
    onset_energy=(max(first["normal_energy"],0.) if first["axial_displacement"]>0. else 0.)+sum(
        max(first[field],0.) for field in ("shear_energy","bending_energy","torsional_energy"))
    radial_traction=2.*onset_energy/separation
    previous=independent_measures(rows[index-1],case)[0] if index else 0.
    frozen_fields=("onset_effective_stress","onset_effective_displacement","failure_effective_displacement","effective_fracture_energy")
    frozen=all(row[field]==first[field] for row in rows[index:] for field in frozen_fields)
    maximum_formula_error=0.
    for row in rows[index:]:
        maximum=row["maximum_effective_displacement"];onset=row["onset_effective_displacement"];failure=row["failure_effective_displacement"]
        expected=1. if maximum>=failure else 0. if maximum<=onset else (1.-onset/maximum)/(1.-onset/failure)
        maximum_formula_error=max(maximum_formula_error,abs(expected-row["damage"]))
    passed=(criterion>=1.-1e-10 and previous<1.+1e-10 and first["damage"]==0.
            and abs(first["onset_effective_stress"]-stress)<=1e-9
            and abs(first["onset_effective_displacement"]-separation)<=1e-10
            and first["maximum_effective_displacement"]==first["onset_effective_displacement"]
            and abs(first["failure_effective_displacement"]-2.*first["effective_fracture_energy"]/stress)<=1e-10
            and frozen and maximum_formula_error<=2e-9)
    return {"started":True,"pass":passed,"sample":index,"amplitude":first["amplitude"],"criterion":criterion,
            "previous_criterion":previous,"first_damage":first["damage"],"onset_effective_stress":first["onset_effective_stress"],
            "onset_effective_displacement":first["onset_effective_displacement"],"failure_effective_displacement":first["failure_effective_displacement"],
            "stress_check_error":abs(first["onset_effective_stress"]-stress),"displacement_check_error":abs(first["onset_effective_displacement"]-separation),
            "onset_radial_traction":radial_traction,"onset_driving_energy":onset_energy,
            "onset_fields_frozen":frozen,"maximum_Eq13_damage_error":maximum_formula_error}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build",type=Path,required=True);parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--backends",choices=("cpu","gpu"),nargs="+",default=("cpu","gpu"),
                        help="Actual backends to run; CPU-only is permitted but does not establish CPU/GPU parity.")
    parser.add_argument("--binary-dir",type=Path,help="Optional exact executable directory, including a multi-configuration build subdirectory.")
    parser.add_argument("--base-refinement",type=int,default=1,
                        help="Positive base sampling factor; the three runs use 2000 times this value times 1, 2 and 4 steps.")
    args=parser.parse_args();source=Path(__file__).resolve().parents[3];output=args.output.resolve()
    if not 1<=args.base_refinement<=(2**31-1)//2000//4:
        parser.error("Base refinement must be a positive integer with representable step counts for all three runs.")
    if output.exists():parser.error("Choose a fresh output directory.")
    if output==source or source in output.parents:parser.error("Keep output outside the source tree.")
    output.mkdir(parents=True);(output/"logs").mkdir()
    binary_root=args.binary_dir.resolve() if args.binary_dir else args.build.resolve()/"src/tests/BondFracture"
    backends=tuple(dict.fromkeys(args.backends))
    relatives=source_files(source);provenance=source_identity(source,relatives)
    snapshot=output/"source"
    for relative in relatives:
        destination=snapshot/relative;destination.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source/relative,destination)
    (output/"source-identity.json").write_text(json.dumps(provenance,indent=2)+"\n")
    runs={};run_rows=[];history_checks={};native_checks={}
    for case in CASES:
        for backend in backends:
            for increment,factor in INCREMENTS:
                actual_factor=args.base_refinement*factor
                name=f"{case}-{backend}-{increment}"
                binary=binary_root/("bond_fracture_"+case.replace("-","_"))
                if platform.system()=="Windows":binary=binary.with_suffix(".exe")
                command=[str(binary),f"--{backend}","--refinement",str(actual_factor),"--output",str(output/name)]
                start=time.monotonic();completed=subprocess.run(command,text=True,capture_output=True);elapsed=time.monotonic()-start
                (output/"logs"/(name+".txt")).write_text(f"command={json.dumps(command)}\nreturncode={completed.returncode}\nelapsed_seconds={elapsed:.9f}\nstdout:\n{completed.stdout}\nstderr:\n{completed.stderr}")
                print(name,completed.returncode,completed.stdout.strip(),completed.stderr.strip(),flush=True)
                checks={"native_exit":completed.returncode==0}
                try:
                    rows=states(output/name/"states.csv");last=rows[-1]
                except (OSError,ValueError,IndexError) as error:
                    (output/"audit.json").write_text(json.dumps({"status":"FAIL","failed_run":name,"error":str(error),
                         "native_exit_code":completed.returncode,"source_identity":provenance,"completed_checks":native_checks},indent=2)+"\n")
                    raise SystemExit("Native history unavailable; failure evidence retained: "+name)
                runs[name]=rows
                checks["sample_count"]=len(rows)==2000*actual_factor+1
                checks["finite"]=all(math.isfinite(value) for row in rows for value in row.values())
                checks["irreversible_damage"]=all(rows[i]["damage"]>=rows[i-1]["damage"]-1e-12 for i in range(1,len(rows)))
                if case!="bending-zero-torque":
                    checks["loading_completion"]=last["damage"]==1. and any(0.<row["damage"]<1. for row in rows)
                else:
                    checks["elastic_history"]=all(row["damage"]==0. and row["damage_started"]==0. for row in rows)
                    checks["native_torque_cancellation"]=all(max(abs(row["torque_"+axis]) for axis in ("x","y","z"))<1e-12 for row in rows)
                    checks["nonzero_endpoint_torque"]=abs(last["master_torque_z"])>1.
                    elastic=sum(last[field] for field in ("normal_energy","shear_energy","bending_energy","torsional_energy"))
                    checks["elastic_work"]=abs(last["cumulative_external_work"]-elastic)<=WORK_ABS_TOL
                dissipation=0.;previous_energy=0.;previous_damage=0.;balance_error=0.;integral_error=0.
                for row in rows:
                    driving=(max(row["trial_normal_energy"],0.) if row["axial_displacement"]>0. else 0.)+sum(
                        max(row[field],0.) for field in ("trial_shear_energy","trial_bending_energy","trial_torsional_energy"))
                    dissipation+=.5*(previous_energy+driving)*(row["damage"]-previous_damage)
                    stored=sum(row[field] for field in ("normal_energy","shear_energy","bending_energy","torsional_energy"))
                    balance_error=max(balance_error,abs(row["cumulative_external_work"]-dissipation-stored))
                    integral_error=max(integral_error,abs(row["cumulative_damage_dissipation"]-dissipation))
                    previous_energy=driving;previous_damage=row["damage"]
                checks["native_energy_balance"]=balance_error<=WORK_ABS_TOL
                checks["independent_damage_integral"]=integral_error<=1e-12
                if case=="tension":checks["fracture_work"]=abs(last["cumulative_external_work"]-GC)<=WORK_ABS_TOL
                if case=="tension":checks["opening_work"]=checks["fracture_work"]
                try:history_checks[name]=history_audit(rows,case)
                except (ValueError,ZeroDivisionError,KeyError) as error:
                    history_checks[name]={"pass":False,"error":str(error)}
                checks["separation_history"]=history_checks[name]["pass"]
                native_checks[name]={"pass":all(checks.values()),"checks":checks,"exit_code":completed.returncode}
                effective_target=sum(last[field] for field in ("normal_energy","shear_energy","bending_energy","torsional_energy")) if case=="bending-zero-torque" else GC
                run_rows.append({"run":name,"case":case,"backend":backend,"increment":increment,
                    "refinement":actual_factor,"parameter_time_step_s":1./(2000*actual_factor),
                    "samples":len(rows),"elapsed_seconds":elapsed,"final_damage":last["damage"],
                    "true_particle_work_J":last["cumulative_external_work"],"effective_traction_area_J_per_m2":last["cumulative_effective_work"],
                    "damage_dissipation_J":dissipation,"maximum_energy_balance_error_J":balance_error,
                    "fracture_budget_work_deviation_J":None if case=="bending-zero-torque" else last["cumulative_external_work"]-GC,
                    "effective_area_target_J_per_m2":effective_target,
                    "effective_area_error":abs(last["cumulative_effective_work"]-effective_target),
                    "effective_area_required":case=="tension",
                    "executable_sha256":sha(binary),"states_sha256":sha(output/name/"states.csv")})
    area_pass=all(row["effective_area_error"]<=WORK_ABS_TOL for row in run_rows if row["effective_area_required"])
    native_pass=all(check["pass"] for check in native_checks.values())
    source_unchanged=provenance["files"]==source_identity(source,relatives)["files"]
    checks={};overall=area_pass and native_pass and source_unchanged;flat=[]
    for case in CASES:
        entry={"parity":{},"refinement":{}}
        if "cpu" in backends and "gpu" in backends:
            for increment,_ in INCREMENTS:
                errors=compare(runs[f"{case}-cpu-{increment}"],runs[f"{case}-gpu-{increment}"])
                passed=all(error["absolute"]<=1e-10 or error["normalized"]<=PARITY_TOL for error in errors.values())
                entry["parity"][increment]={"pass":passed,"errors":errors};overall &= passed
        for backend in backends:
            for coarse,fine,stride in (("base","half",2),("half","quarter",2),("base","quarter",4)):
                errors=compare(runs[f"{case}-{backend}-{coarse}"],runs[f"{case}-{backend}-{fine}"],stride)
                passed=all(error["absolute"]<=1e-10 or error["normalized"]<=REFINEMENT_TOL for error in errors.values())
                entry["refinement"][f"{backend}-{coarse}-{fine}"]={"pass":passed,"errors":errors};overall &= passed
        checks[case]=entry
        for group in ("parity","refinement"):
            tolerance=PARITY_TOL if group=="parity" else REFINEMENT_TOL
            for comparison,check in entry[group].items():
                for field,error in check["errors"].items():
                    flat.append({"case":case,"group":group,"comparison":comparison,"field":field,"absolute_error":error["absolute"],
                        "normalized_error":error["normalized"],"normalized_tolerance":tolerance,
                        "pass":error["absolute"]<=1e-10 or error["normalized"]<=tolerance})
    parity_ran="cpu" in backends and "gpu" in backends
    parity_pass=parity_ran and all(check["pass"] for entry in checks.values() for check in entry["parity"].values())
    refinement_pass=all(check["pass"] for entry in checks.values() for check in entry["refinement"].values())
    def verdict(passed):return "PASS" if passed else "FAIL"
    def area_verdict(case_rows):
        required=[row for row in case_rows if row["effective_area_required"]]
        return verdict(all(row["effective_area_error"]<=WORK_ABS_TOL for row in required)) if required else "DIAGNOSTIC (section-stress failure law)"
    parity_verdict=verdict(parity_pass) if parity_ran else "NOT RUN"
    def native_verdict(key):
        relevant=[check["checks"][key] for check in native_checks.values() if key in check["checks"]]
        return verdict(bool(relevant) and all(relevant))
    audit={"status":"PASS" if overall else "FAIL","scope":"native single-point strength-initiated Euclidean-separation BK adaptation, not paper geometric gap",
           "source_identity":provenance,"source_unchanged_during_campaign":source_unchanged,
           "backends":backends,"base_refinement":args.base_refinement,
           "increments":{name:1./(2000*args.base_refinement*factor) for name,factor in INCREMENTS},
           "native_run_completion":verdict(native_pass),"effective_traction_area":verdict(area_pass),
           "pure_opening_true_work":native_verdict("opening_work"),"zero_torque_elastic_work":native_verdict("elastic_work"),"pure_opening_true_work_equals_A_Gc":native_verdict("fracture_work"),
           "all_native_work_energy_balance":native_verdict("native_energy_balance"),"independent_damage_dissipation":native_verdict("independent_damage_integral"),
           "section_stress_fracture_budget":"DIAGNOSTIC outside pure axial opening: current df uses the section-stress norm; actual deviations remain in runs.csv",
           "stress_separation_history":native_verdict("separation_history"),"CPU_GPU_parity":parity_verdict,
           "loading_increment_refinement":"PASS" if refinement_pass else "FAIL","sample_count":sum(len(rows) for rows in runs.values()),
           "tolerances":{"CPU_GPU_normalized":PARITY_TOL,"loading_increment_normalized":REFINEMENT_TOL,"absolute_noise_floor":1e-10,
                         "work_absolute":WORK_ABS_TOL},"checks":checks,"native_checks":native_checks,"separation_history":history_checks,
           "physical_validation":"NOT PERFORMED"}
    (output/"audit.json").write_text(json.dumps(audit,indent=2)+"\n")
    write_csv(output/"runs.csv",run_rows);write_csv(output/"comparisons.csv",flat)
    onset_rows=[]
    for name,check in history_checks.items():
        onset_rows.append({key:({"run":name}|check).get(key,"") for key in ("run","started","pass","sample","amplitude","criterion","previous_criterion",
            "first_damage","onset_effective_stress","onset_effective_displacement","failure_effective_displacement","maximum_Eq13_damage_error")})
    write_csv(output/"onset.csv",onset_rows)
    lines=["# FunDEM Euclidean-separation BK adaptation verification","",
           f"Campaign verdict: **{audit['status']}**. {len(run_rows)} actual native runs retain {audit['sample_count']:,} states.","",
           f"Source HEAD: `{provenance['git_head'] or 'unavailable'}`; tracked tree dirty: {provenance['tracked_tree_dirty']}; exact tested source digest: `{provenance['source_sha256']}`. See [source identity](source-identity.json) and the pre-run source snapshot.","",
           f"Native completion: {audit['native_run_completion']}; separation history: {audit['stress_separation_history']}; pure-opening mechanical work: {audit['pure_opening_true_work']}. Selected effective-area check: {audit['effective_traction_area']}; CPU/GPU parity: {audit['CPU_GPU_parity']}; loading-increment refinement: {audit['loading_increment_refinement']}.","",
           "Numerical consistency and loading-increment comparisons are not physical validation. Any failed energy/refinement checks remain FAIL and make this strict campaign return nonzero; historical campaigns are not overwritten.","",
           "## Formulation and scope","",
           "The strength criterion uses the unchanged native B3/B4 torsion/bending torque and an equivalent circular section. "
           "The undamaged damage-driving energies are UN=Un_opening+Ub and US=Us+Ut; axial compression is excluded. "
           "The components are dn=2UN/(A*sigma_open) and ds=2US/(A*tau), with a zero component when its denominator is nonpositive; de=hypot(dn,ds) is their ordinary Euclidean norm. "
           "Here sigma_open=max(axial_stress,0)+bending_stress; the signed axial-plus-bending stress remains the strength-criterion input. "
           "Native force/torque expressions are retained; modal energies use squared-vector evaluation, raw bending stiffness and a nonnegative torsional-energy clamp. "
           "The adaptation is inspired by the cited strength/softening construction but is not an exact reproduction of its multi-point model.","",
           "At the first sampled Q>=1 state the opening/shear section-stress norm sigma0, d0 and BK Gc are frozen. "
           "The maximum Euclidean separation is irreversible; sigma0=hypot(sigma_open,tau) gives df=2*Gc/sigma0. "
           "D=(1-d0/dmax)/(1-d0/df) is bounded by [0,1] and previous damage. "
           "If df<=d0, the bond fails immediately and retains the requested Gc; no tolerance-based energy inflation is permitted.","",
           "BK with eta uses UN/US at onset only. Positive opening, shear, bending and torsional response use (1-D); partial damage leaves the axial compressive branch undamaged. "
           "A fully broken bond is removed and native contact subsequently carries compression. These prescribed fracture paths contain no axial compression. "
           "The paper's multi-point Voronoi interfaces, damping and specimen experiments are not reproduced.","",
           "The selected actual CPU bond::calculateForce and/or CUDA bondForceKernel are evaluated at prescribed poses. "
           f"Base sampling refinement is {args.base_refinement}; parameter-time increments are h={audit['increments']['base']:.12g}, h/2={audit['increments']['half']:.12g} and h/4={audit['increments']['quarter']:.12g} s; "
           "there is no dynamic integration. Host-reconstructed undamaged energy components and stress/separation measurements support postprocessing; native histories, forces and torques are separately compared.","",
           "## Parameters and torque cancellation","",
           "L=1 m; A=1 m^2; kn=1000 N/m; ks=400 N/m; kt=100 N m/rad; tensile strength=cohesion=1 Pa; friction=0; GIc=GIIc=0.005 J/m^2; eta=1.75. "
           "Only nondegenerate bending uses kb=200 N m/rad. Other cases preserve kb=100. The original kb=100 bending fixture has native TBT cancellation and remains elastic. "
           "These are verification parameters, not calibrated rock properties.","",
           "## Two distinct work measures","",
           "True particle external work integrates the original particle forces/torques against prescribed body motion. "
           "Radial scalar traction area integrates (1-D)*2(UN+US)/(A*de) against de, in J/m^2. "
           "Only pure axial opening is checked against true work=A*Gc and radial area=Gc. "
           "Other modes use a section-stress failure length rather than the radial work-conjugate traction: true work must balance independently integrated Ua*dD plus remaining stored energy; actual A*Gc and scalar-area deviations remain diagnostics. "
           "Zero-TBT bending remains elastic and true work must equal stored energy; its radial scalar area is likewise only diagnostic. "
           "C++ regression independently verifies the two-spring path against its current triangular work 0.5*k_path*d0*df and retains strict native external-work/energy-balance checks.","",
           "| Case / backend | True work h / h2 / h4 (J) | Radial area h / h2 / h4 (J/m^2) | True work minus A Gc (J) | Increment refinement |",
           "| --- | ---: | ---: | ---: | --- |"]
    for case in CASES:
        for backend in backends:
            final=[runs[f"{case}-{backend}-{name}"][-1] for name,_ in INCREMENTS]
            true_work=" / ".join(f"{row['cumulative_external_work']:.12g}" for row in final)
            selected_area=" / ".join(f"{row['cumulative_effective_work']:.12g}" for row in final)
            budget_deviation="not a fracture path" if case=="bending-zero-torque" else " / ".join(f"{row['cumulative_external_work']-GC:.12g}" for row in final)
            passed=all(check["pass"] for name,check in checks[case]["refinement"].items() if name.startswith(backend+"-"))
            lines.append(f"| {case} / {backend.upper()} | {true_work} | {selected_area} | {budget_deviation} | {verdict(passed)} |")
    lines += ["","## Compression and degenerate onset","",
              "Axial compression is excluded before summing the modal damage-driving energies, so it cannot cancel positive bending energy. "
              "Positive driving energies below the shared geometric length tolerance remain valid; modal energies use the current squared-vector expressions. "
              "Zero driving separation at a strength crossing or df<=d0 produces immediate brittle failure without division by zero or a fabricated Gc. "
              "Partial damage retains axial compressive force and energy; complete bond removal leaves subsequent compression to the separate native contact law. "
              "CPU/GPU, energy or refinement discrepancies remain FAIL; successful numerical checks alone do not establish physical validity. Historical evidence is not rewritten.","",
              "## Independent checks and limits","",
              "Every valid softening first-crossing state is checked independently from its actual force, native torque and pose: prior Q<1, current Q>=1, D=0, correct sigma0/d0/section-stress df, and frozen onset fields. "
              "Every recorded damaged state is compared against the irreversible Euclidean-separation softening law. Separate C++ tests cover pure and mixed native work, section-stress component recovery, different section areas, unloading, compression retention, insufficient budgets, "
              "all four energy contributions, raw stiffness lifecycle, exact native torque, rigid rotation and BK exponent/frozen mixity.","",
              f"Unchanged limits are {PARITY_TOL:g} normalized CPU/GPU error, {REFINEMENT_TOL:g} loading-increment error, 1e-10 absolute noise floor, and {WORK_ABS_TOL:g} effective-area/work tolerance. "
              "Both adjacent refinements and h versus h/4 are compared at matching poses. This records three discretizations, not an assumed convergence order. "
              "Discrete sampled onset can change sigma0/d0/df on halving the increment; any failed comparisons remain FAIL. "
              "No specimen strength, crack pattern, non-proportional mixed-path fracture calibration, or dynamic timestep convergence is validated.","",
               "Refinement errors below are maxima over matching coarse/fine poses. Near-zero components retain the absolute noise floor; component-normalized percentages do not describe the entire load vector.","",
               "| Case / comparison | Onset stress error (%) | Onset separation error (%) | Failure separation error (%) | Maximum absolute damage difference | Verdict |",
               "| --- | ---: | ---: | ---: | ---: | --- |"]
    for case in CASES:
        for name,check in checks[case]["refinement"].items():
            errors=check["errors"]
            lines.append(f"| {case} / {name} | {100*errors['onset_effective_stress']['normalized']:.6g} | {100*errors['onset_effective_displacement']['normalized']:.6g} | {100*errors['failure_effective_displacement']['normalized']:.6g} | {errors['damage']['absolute']:.9g} | {verdict(check['pass'])} |")
    lines += ["",
              "Raw evidence: [runs](runs.csv), [all comparisons](comparisons.csv), [onset states](onset.csv), [audit](audit.json), and [file hashes](SHA256SUMS.csv). "
              "Previous campaigns remain immutable historical evidence and are superseded for this Euclidean-separation formulation.","",
              "Reference: Suarez et al. (2025), *Multiscale modeling of rock fracture in comminution — A comparative study of FEM accuracy and DEM scalability*, Minerals Engineering 232, 109488. https://doi.org/10.1016/j.mineng.2025.109488."]
    (output/"REPORT.md").write_text("\n".join(lines)+"\n")
    for case in CASES:
        rows=[row for row in run_rows if row["case"]==case]
        parity=verdict(all(check["pass"] for check in checks[case]["parity"].values())) if parity_ran else "NOT RUN"
        refine=all(check["pass"] for check in checks[case]["refinement"].values())
        case_lines=[f"# Euclidean-separation BK case: {case}","",
                    "This is a native FunDEM single-point strength-initiated Euclidean-separation adaptation, not an official LAMMPS/PFC demo or exact paper specimen reproduction.","",
                    "Native CPU/GPU kernels are evaluated at prescribed poses. True work is checked against damage dissipation plus stored energy. The frozen fracture-budget equality is required only for pure axial opening; all other modes retain their actual budget deviation as diagnostic evidence under the section-stress failure law. "
                    "See [REPORT.md](REPORT.md) for definitions, parameters, stiffness distinction and limitations.","",
                    f"CPU/GPU parity: **{parity}**; loading-increment refinement: **{'PASS' if refine else 'FAIL'}**; "
                    f"radial-area budget check: **{area_verdict(rows)}**.","",
                    "| Run | States | True work (J) | Selected area (J/m^2) | Raw history |","| --- | ---: | ---: | ---: | --- |"]
        for row in rows:
            name=row["run"];case_lines.append(f"| {name} | {row['samples']} | {row['true_particle_work_J']:.12g} | {row['effective_traction_area_J_per_m2']:.12g} | [states.csv]({name}/states.csv) |")
        case_lines += ["","Radial-area interpretation: "+("target 0.005 J/m^2 for pure axial opening. The path passes through elastic, partially damaged and completely failed states." if case=="tension" else "diagnostic only under the section-stress failure law; true work retains its independent stored-energy/dissipation check."),
                       "All original numerical tolerances and any refinement failures are retained in [comparisons.csv](comparisons.csv).",""]
        (output/f"REPORT_{case}.md").write_text("\n".join(case_lines))
    manifest=[{"path":str(path.relative_to(output)),"bytes":path.stat().st_size,"sha256":sha(path)}
              for path in sorted(output.rglob("*")) if path.is_file() and path.name!="SHA256SUMS.csv"]
    write_csv(output/"SHA256SUMS.csv",manifest)
    print("Campaign",audit["status"],output,flush=True)
    return 0 if overall else 1
if __name__=="__main__":raise SystemExit(main())
