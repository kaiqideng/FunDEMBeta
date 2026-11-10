// Prescribed-pose native bond verification with independent radial/external work measurements.
#pragma once
#include "solver/solvers.h"
#include "execution/bondFunctions.h"
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
#include "execution/cuda/bondForceKernel.cuh"
#endif
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace cohesive_validation {
using namespace fundem;
using math::Real;
using math::Vec3;
using math::Quaternion;
inline void require(bool ok,const char* message) { if(!ok) throw std::runtime_error(message); }
inline const char* modeName(int mode) { return mode==0?"tension":mode==1?"shear":mode==2?"bending":mode==3?"torsion":"bending-zero-torque"; }
inline Real maximumAmplitude(int mode) { return mode==0?.012:mode==1?.06:(mode==2||mode==4)?.1:.05; }
inline Real totalElasticEnergy(const bond& b) {
    return b.normalElasticEnergy()+b.shearElasticEnergy()+b.bendingElasticEnergy()+b.torsionalElasticEnergy();
}
struct Measures {
    Real criterion{},sigma{},tau{},effectiveStress{},effectiveDisplacement{};
    Real axial{},un{},us{},ub{},ut{},da{},db{},ds{},dtor{},dn{},dt{};
};
struct Experiment {
    materialContainer materials;
    particleContainer particles;
    bondContainer bonds;
    bool gpu{};
    explicit Experiment(bool device,int mode):gpu(device) {
        materials.host().push_back(material{1.e6,.25,1.,1000.});
        for(int i=0;i<2;++i) {
            particle p; p.setRadius(.1); p.setMaterial(materials,0); p.setPosition({Real(i),0.,0.});
            particles.host().push_back(p);
        }
        bond b{1.};
        require(b.setConnection(particles,0,1),"Endpoint connection rejected.");
        require(b.setStiffness(1000.,400.,mode==2?200.:100.,100.),"Stiffness rejected.");
        require(b.setCrossSectionArea(1.),"Area rejected.");
        require(b.setModeICriticalEnergy(.005)&&b.setModeIICriticalEnergy(.005),"Fracture energies rejected.");
        require(b.setTensileStrength(1.)&&b.setCohesion(1.)&&b.setInternalFrictionCoefficient(0.),"Strengths rejected.");
        bonds.host().push_back(b);
        particles.setUseDevice(gpu); bonds.setUseDevice(gpu);
#if !defined(FUNDEM_HAS_CUDA) || !FUNDEM_HAS_CUDA
        require(!gpu,"GPU requested but this build has no CUDA backend.");
#endif
    }
    void setPose(int mode,Real a) {
        auto& m=particles.host()[0];auto& s=particles.host()[1];
        m.setPosition({0.,0.,0.});s.setPosition({1.,0.,0.});
        m.setOrientation(Quaternion::identity());s.setOrientation(Quaternion::identity());
        m.setForce(Vec3::zero());s.setForce(Vec3::zero());m.setTorque(Vec3::zero());s.setTorque(Vec3::zero());
        if(mode==0)s.setPosition({1.+a,0.,0.});
        if(mode==1)s.setPosition({std::cos(a),std::sin(a),0.});
        if(mode==2||mode==4) {
            m.setOrientation(Quaternion::fromUnitAxisAngle(Vec3::unitZ(),-.5*a));
            s.setOrientation(Quaternion::fromUnitAxisAngle(Vec3::unitZ(),.5*a));
        }
        if(mode==3)s.setOrientation(Quaternion::fromUnitAxisAngle(Vec3::unitX(),a));
    }
    Measures undamagedMeasures() const {
        // Host postprocessing from prescribed poses. Native CPU/GPU histories are audited separately.
        const auto& b=bonds.host()[0];const auto& m=particles.host()[0];const auto& s=particles.host()[1];
        const auto qm=m.orientation(),qs=s.orientation();
        Vec3 force,torque,cm,cs,normal;
        Real axial,un,us,ub,ut;
        execution::calculateBondForceAndTorque(force,torque,cm,cs,normal,axial,un,us,ub,ut,
            m.position()+math::rotateUnit(qm,b.masterEndpointLocalPosition()),
            s.position()+math::rotateUnit(qs,b.slaveEndpointLocalPosition()),
            math::rotateUnit(qm,b.masterEndpointLocalNormal()),
            math::rotateUnit(qm,b.masterEndpointLocalTangent1()),
            math::rotateUnit(qm,b.masterEndpointLocalTangent2()),
            math::rotateUnit(qs,b.slaveEndpointLocalNormal()),
            math::rotateUnit(qs,b.slaveEndpointLocalTangent1()),
            math::rotateUnit(qs,b.slaveEndpointLocalTangent2()),
            b.coefficientB1(),b.coefficientB2(),b.coefficientB3(),b.coefficientB4(),b.bendingStiffness(),b.equivalentLength());
        Measures result;
        Real axialStress=0.,bendingStress=0.;
        result.criterion=execution::calculateBondStrengthCriterion(axialStress,bendingStress,result.tau,force,torque,normal,
            b.crossSectionArea(),b.tensileStrength(),b.cohesion(),b.internalFrictionCoefficient());
        result.sigma=axialStress+bendingStress;
        const Real opening=(axialStress>0.?axialStress:0.)+bendingStress;
        result.axial=axial;result.un=un;result.us=us;result.ub=ub;result.ut=ut;
        const Real openingDenominator=b.crossSectionArea()*opening;
        const Real shearDenominator=b.crossSectionArea()*result.tau;
        result.da=axialStress>0.&&openingDenominator>0.?2.*std::max(un,0.)/openingDenominator:0.;
        result.db=openingDenominator>0.?2.*std::max(ub,0.)/openingDenominator:0.;
        result.ds=shearDenominator>0.?2.*std::max(us,0.)/shearDenominator:0.;
        result.dtor=shearDenominator>0.?2.*std::max(ut,0.)/shearDenominator:0.;
        result.dn=result.da+result.db;
        result.dt=result.ds+result.dtor;
        result.effectiveStress=std::sqrt(opening*opening+result.tau*result.tau);
        result.effectiveDisplacement=std::sqrt(result.dn*result.dn+result.dt*result.dt);
        return result;
    }
    void evaluate() {
        if(!gpu) {
            auto& b=bonds.host()[0];
            require(b.updateMasterSlaveParticle(particles),"Could not refresh body snapshots.");
            require(b.calculateForce(),"Native CPU bond evaluation failed.");
        } else {
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
            particles.copyHostToDevice();bonds.copyHostToDevice();
            cuda::launchBondForce(bonds,particles);
            bonds.copyDeviceToHost();particles.copyDeviceToHost();
#endif
        }
    }
};
inline void vec(std::ostream& out,const Vec3& v){out<<','<<v.x<<','<<v.y<<','<<v.z;}
inline void quat(std::ostream& out,const Quaternion& q){out<<','<<q.w<<','<<q.x<<','<<q.y<<','<<q.z;}
inline int run(int argc,char** argv,int mode) {
    try {
        bool gpu=false,smoke=false;int refinement=1;std::filesystem::path output;
        for(int i=1;i<argc;++i) {
            const std::string arg=argv[i];
            if(arg=="--gpu")gpu=true;
            else if(arg=="--cpu")gpu=false;
            else if(arg=="--half-step")refinement=2;
            else if(arg=="--quarter-step")refinement=4;
            else if(arg=="--refinement"&&i+1<argc) {
                const std::string value=argv[++i];std::size_t consumed=0;
                refinement=std::stoi(value,&consumed);
                require(consumed==value.size()&&refinement>0&&refinement<=std::numeric_limits<int>::max()/2000,
                    "Refinement must be a positive integer with a representable step count.");
            }
            else if(arg=="--smoke")smoke=true;
            else if(arg=="--output"&&i+1<argc)output=argv[++i];
            else if(arg=="--help"){std::cout<<"--cpu|--gpu [--refinement N|--half-step|--quarter-step] [--smoke] --output NEW_DIRECTORY\n";return 0;}
            else throw std::runtime_error("Unknown or incomplete option: "+arg);
        }
        require(!output.empty(),"Choose a new --output directory outside the source tree.");
        require(!std::filesystem::exists(output),"Refusing to overwrite an existing run directory.");
        std::filesystem::create_directories(output);
        Experiment experiment(gpu,mode);
        const int steps=2000*refinement;const Real dt=1./steps;
        std::ofstream configuration(output/"configuration.txt");
        configuration<<std::setprecision(17)
          <<"model=single-point strength-initiated Euclidean-separation BK adaptation\nexperiment="<<modeName(mode)
          <<"\nbackend="<<(gpu?"native CUDA bondForceKernel":"native CPU bond::calculateForce")
          <<"\ncontrol=prescribed pose; no dynamic time integration\n"
          <<"verification_scope="<<(smoke?"constitutive smoke; not increment convergence":"single increment; campaign comparison required")<<'\n'
          <<"duration=1 s\ntimeStep="<<dt<<"\nsteps="<<steps<<"\nrefinement="<<refinement
          <<"\nlength=1 m\narea=1 m^2\nnormalStiffness=1000 N/m\nshearStiffness=400 N/m\n"
          <<"bendingStiffness="<<(mode==2?200.:100.)<<" N m/rad\ntorsionalStiffness=100 N m/rad\n"
          <<"tensileStrength=1 Pa\ncohesion=1 Pa\ninternalFrictionCoefficient=0\n"
          <<"GIc=0.005 J/m^2\nGIIc=0.005 J/m^2\nBK_exponent=1.75\nmaximumAmplitude="<<maximumAmplitude(mode)
          <<"\nparticle_radius=0.1 m\nmaterial_E=1e6 Pa\nmaterial_nu=0.25\nmaterial_density=1000 kg/m^3\nmaterial_restitution=1\n"
          <<"onset=first sampled Q>=1; freeze section stress norm sigma_e0 and d0; df=2Gc/sigma_e0; df<=d0 fails without increasing Gc\n"
          <<"separation=dn=2(Un_opening+Ub)/(A*sigma_open); ds=2(Us+Ut)/(A*tau); zero stress gives zero component; de=hypot(dn,ds)\n"
          <<"Stress and separation measurements are undamaged host reconstructions from the prescribed pose; native histories/loads are recorded separately.\n"
          <<"True external work and reconstructed radial traction area are separate outputs; the current section-stress failure length does not generally impose work=A*Gc.\n"
          <<"Tensile/shear/rotational response uses (1-D); intact compressive response remains until bond removal; torque is the existing B3/B4 expression.\n";
        std::ofstream csv(output/"states.csv");
        csv<<std::setprecision(17)
           <<"step,time,amplitude,damage,damage_started,onset_effective_stress,onset_effective_displacement"
           <<",maximum_effective_displacement,failure_effective_displacement,effective_fracture_energy"
           <<",normal_energy,shear_energy,bending_energy,torsional_energy"
           <<",force_x,force_y,force_z,torque_x,torque_y,torque_z,master_torque_x,master_torque_y,master_torque_z"
           <<",slave_torque_x,slave_torque_y,slave_torque_z,master_x,master_y,master_z,slave_x,slave_y,slave_z"
           <<",master_qw,master_qx,master_qy,master_qz,slave_qw,slave_qx,slave_qy,slave_qz"
           <<",trial_normal_stress,trial_shear_stress,trial_criterion,effective_trial_stress,effective_displacement"
           <<",axial_displacement,trial_normal_energy,trial_shear_energy,trial_bending_energy,trial_torsional_energy"
           <<",axial_separation,bending_separation,shear_separation,torsional_separation,normal_separation,tangential_separation"
           <<",effective_traction,work_increment,cumulative_external_work,effective_work_increment,cumulative_effective_work"
           <<",damage_dissipation_increment,cumulative_damage_dissipation\n";
        Real previousAmplitude=0.,previousGeneralizedForce=0.,work=0.,previousDamage=0.;
        Real previousSeparation=0.,previousEffectiveTraction=0.,effectiveWork=0.,previousTrialEnergy=0.,dissipated=0.;
        bool started=false,partial=false;
        for(int i=0;i<=steps;++i) {
            const Real time=i*dt,amplitude=maximumAmplitude(mode)*time;
            experiment.setPose(mode,amplitude);
            const Measures measure=experiment.undamagedMeasures();
            experiment.evaluate();const auto& b=experiment.bonds.host()[0];
            const Real generalized=mode==0?b.force().x:
                mode==1?math::dot(b.force(),Vec3{-std::sin(amplitude),std::cos(amplitude),0.}):
                (mode==2||mode==4)?.5*(b.masterTorque().z-b.slaveTorque().z):-b.slaveTorque().x;
            const Real increment=.5*(previousGeneralizedForce+generalized)*(amplitude-previousAmplitude);
            const Real trialEnergy=(measure.axial>0.?std::max(measure.un,0.):0.)
                +std::max(measure.us,0.)+std::max(measure.ub,0.)+std::max(measure.ut,0.);
            const Real traction=measure.effectiveDisplacement>0.?
                (1.-b.damageFactor())*(2.*trialEnergy)/(b.crossSectionArea()*measure.effectiveDisplacement):0.;
            const Real effectiveIncrement=.5*(previousEffectiveTraction+traction)*(measure.effectiveDisplacement-previousSeparation);
            const Real dissipationIncrement=.5*(previousTrialEnergy+trialEnergy)*(b.damageFactor()-previousDamage);
            work+=increment;effectiveWork+=effectiveIncrement;dissipated+=dissipationIncrement;previousTrialEnergy=trialEnergy;
            require(std::isfinite(b.damageFactor())&&b.damageFactor()+1.e-12>=previousDamage,"Invalid or healing damage.");
            previousDamage=b.damageFactor();started|=b.getDamageStarted();partial|=b.damageFactor()>0.&&b.damageFactor()<1.;
            csv<<i<<','<<time<<','<<amplitude<<','<<b.damageFactor()<<','<<b.getDamageStarted()
               <<','<<b.getOnsetEffectiveStress()<<','<<b.getOnsetEffectiveDisplacement()<<','<<b.getMaximumEffectiveDisplacement()
               <<','<<b.getFailureEffectiveDisplacement()<<','<<b.getEffectiveFractureEnergy()
               <<','<<b.normalElasticEnergy()<<','<<b.shearElasticEnergy()<<','<<b.bendingElasticEnergy()<<','<<b.torsionalElasticEnergy();
            vec(csv,b.force());vec(csv,b.torque());vec(csv,b.masterTorque());vec(csv,b.slaveTorque());
            for(const auto& p:experiment.particles.host())vec(csv,p.position());
            for(const auto& p:experiment.particles.host())quat(csv,p.orientation());
            csv<<','<<measure.sigma<<','<<measure.tau<<','<<measure.criterion<<','<<measure.effectiveStress<<','<<measure.effectiveDisplacement
               <<','<<measure.axial<<','<<measure.un<<','<<measure.us<<','<<measure.ub<<','<<measure.ut
               <<','<<measure.da<<','<<measure.db<<','<<measure.ds<<','<<measure.dtor<<','<<measure.dn<<','<<measure.dt
               <<','<<traction<<','<<increment<<','<<work<<','<<effectiveIncrement<<','<<effectiveWork
               <<','<<dissipationIncrement<<','<<dissipated<<'\n';
            previousAmplitude=amplitude;previousGeneralizedForce=generalized;
            previousSeparation=measure.effectiveDisplacement;previousEffectiveTraction=traction;
        }
        const auto& last=experiment.bonds.host()[0];
        if(mode!=4) {
            require(started&&partial&&last.damageFactor()==1.,"Loading did not traverse elastic, softening, and failed states.");
            require(math::norm(last.force())==0.&&math::norm(last.masterTorque())==0.&&math::norm(last.torque())==0.,"Failed bond retained a load.");
            require(std::abs(work-dissipated)<2.e-7,"True native external work must balance the independently integrated damage dissipation.");
            if(mode==0) {
                require(std::abs(work-last.crossSectionArea()*last.getEffectiveFractureEnergy())<2.e-7,"Pure axial opening must retain its A Gc work budget.");
                require(std::abs(effectiveWork-last.getEffectiveFractureEnergy())<2.e-7,"Pure axial opening scalar work must retain its frozen Gc budget.");
            }
        } else {
            require(!started&&!partial&&last.damageFactor()==0.,"Zero-TBT bending must remain elastic.");
            require(math::norm(last.torque())<1.e-12&&math::norm(last.masterTorque())>1.,"Expected native torque cancellation absent.");
            const Real elastic=last.normalElasticEnergy()+last.shearElasticEnergy()+last.bendingElasticEnergy()+last.torsionalElasticEnergy();
            require(std::abs(work-elastic)<2.e-7,"Zero-TBT true external elastic work must equal its stored energy.");
        }
        require(std::abs(work-(dissipated+totalElasticEnergy(last)))<2.e-7,"Native work must balance damage dissipation plus final stored energy.");
        std::ofstream summary(output/"summary.txt");
        summary<<std::setprecision(17)<<"COMPLETED\ncase="<<modeName(mode)<<"\nbackend="<<(gpu?"GPU":"CPU")
               <<"\nverification_scope="<<(smoke?"constitutive smoke":"single increment")
               <<"\nsteps="<<steps<<"\ndt="<<dt<<"\nrefinement="<<refinement<<"\nfinal_damage="<<last.damageFactor()
               <<"\ntrue_particle_external_work="<<work<<"\nselected_effective_traction_area="<<effectiveWork
               <<"\ncumulative_damage_dissipation="<<dissipated
               <<"\nfracture_budget_work_deviation="<<(mode==4?0.:work-last.crossSectionArea()*last.getEffectiveFractureEnergy())
               <<"\neffective_area_tolerance=2e-7\neffective_area_check="<<(mode!=0?"DIAGNOSTIC_SECTION_STRESS_LAW":std::abs(effectiveWork-last.getEffectiveFractureEnergy())<2.e-7?"PASS":"FAIL")
               <<"\nTrue work is checked against damage dissipation plus stored energy; equality to A*Gc is checked only for pure axial opening.\n"
               <<"Completion does not assert effective-area agreement, increment convergence or physical validation; use runVerification.py for strict comparison.\n";
        std::cout<<modeName(mode)<<' '<<(gpu?"GPU":"CPU")<<" COMPLETED, true_work="<<std::setprecision(10)<<work
                 <<", selected_area="<<effectiveWork<<'\n';
        return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
} // namespace cohesive_validation
