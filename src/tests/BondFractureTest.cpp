// Focused constitutive checks for the single-point strength-initiated BK bond.
#include "solver/solvers.h"
#include "execution/bondFunctions.h"
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
#include "execution/cuda/bondForceKernel.cuh"
#endif

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
using namespace fundem;
using math::Real;
using math::Vec3;
using math::Quaternion;

void require(bool ok, const char* message) { if (!ok) throw std::runtime_error(message); }
bool close(Real a, Real b, Real tol=1.e-10) { return std::abs(a-b)<=tol*(1.+std::max(std::abs(a),std::abs(b))); }
struct Fixture {
    materialContainer materials;
    particleContainer particles;
    bondContainer bonds;
    Fixture() {
        materials.host().push_back(material{1.e6,.25,1.,1000.});
        for (int i=0;i<2;++i) {
            particle p; p.setRadius(.1); p.setMaterial(materials,0); p.setPosition({Real(i),0.,0.});
            particles.host().push_back(p);
        }
        bond b{1.};
        require(b.setConnection(particles,0,1),"Configure endpoints");
        require(b.setStiffness(1000.,400.,100.,100.),"Configure stiffness");
        require(b.setCrossSectionArea(1.),"Configure area");
        require(b.setModeICriticalEnergy(.005) && b.setModeIICriticalEnergy(.005),"Configure energies");
        require(b.setTensileStrength(1.) && b.setCohesion(1.) && b.setInternalFrictionCoefficient(0.),"Configure strengths");
        bonds.host().push_back(b);
    }
    bond& value() { return bonds.host().front(); }
    void evaluate() {
        require(value().updateMasterSlaveParticle(particles),"Refresh snapshots");
        require(value().calculateForce(),"Evaluate native bond");
    }
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
    void evaluateGPU() {
        particles.copyHostToDevice(); bonds.copyHostToDevice();
        cuda::launchBondForce(bonds,particles);
        bonds.copyDeviceToHost(); particles.copyDeviceToHost();
    }
#endif
    void pose(int mode, Real a) {
        auto& m=particles.host()[0]; auto& s=particles.host()[1];
        m.setPosition({0.,0.,0.}); s.setPosition({1.,0.,0.});
        m.setOrientation(Quaternion::identity()); s.setOrientation(Quaternion::identity());
        if(mode==0) s.setPosition({1.+a,0.,0.});
        if(mode==1) s.setPosition({std::cos(a),std::sin(a),0.});
        if(mode==2) { m.setOrientation(Quaternion::fromUnitAxisAngle(Vec3::unitZ(),-.5*a)); s.setOrientation(Quaternion::fromUnitAxisAngle(Vec3::unitZ(),.5*a)); }
        if(mode==3) s.setOrientation(Quaternion::fromUnitAxisAngle(Vec3::unitX(),a));
    }
};
Real total(const bond& b) { return b.normalElasticEnergy()+b.shearElasticEnergy()+b.bendingElasticEnergy()+b.torsionalElasticEnergy(); }
void requireFiniteEnergies(const bond& b) {
    for(Real energy : {b.normalElasticEnergy(),b.shearElasticEnergy(),b.bendingElasticEnergy(),b.torsionalElasticEnergy()})
        require(std::isfinite(energy),"Every stored modal energy must remain finite for a valid finite pose");
}

void testConfiguration() {
    bond b;
    require(!b.setTensileStrength(-1.) && !b.setCohesion(-1.) && !b.setInternalFrictionCoefficient(-1.),"Reject negative strengths");
    require(!b.setTensileStrength(std::numeric_limits<Real>::quiet_NaN()),"Reject nonfinite strength");
}
void testPureModes() {
    const Real maximum[]{.012,.06,.1,.05};
    for(int mode=0;mode<4;++mode) {
        Fixture f; bool started=false, partial=false;
        if(mode==2) require(f.value().setStiffness(1000.,400.,200.,100.),"Use nondegenerate native bend torque for fracture");
        Real previous=0.;
        for(int i=0;i<=4000;++i) {
            f.pose(mode,maximum[mode]*i/4000.); f.evaluate();
            const auto& b=f.value();
            require(b.damageFactor()+1.e-13>=previous,"Damage may not heal during loading");
            previous=b.damageFactor();
            started|=b.getDamageStarted()!=0;
            partial|=b.damageFactor()>0. && b.damageFactor()<1.;
        }
        require(started && partial,"Each mode must trigger onset and pass a softening interval");
        require(f.value().damageFactor()==1.,"Each pure loading mode must reach complete fracture");
        require(math::norm(f.value().force())==0. && math::norm(f.value().masterTorque())==0. && math::norm(f.value().torque())==0. && total(f.value())==0.,"Newly failed bond must clear loads and stored energy immediately");
    }
}
void testStrengthAndWork() {
    Fixture f;
    f.pose(0,.0009); f.evaluate();
    require(!f.value().getDamageStarted() && f.value().damageFactor()==0.,"Substrength opening must remain elastic");
    f.pose(0,.0014); f.evaluate();
    require(f.value().getDamageStarted(),"First sampled strength crossing must start damage history");
    require(close(f.value().getOnsetEffectiveStress(),1.4,1.e-10),"First crossing must store effective trial stress");
    require(close(f.value().getOnsetEffectiveDisplacement(),.0014,1.e-10),"Opening onset must store the four-mode energy-equivalent separation");
    require(close(f.value().getFailureEffectiveDisplacement(),.01/1.4,1.e-10),"Axial rupture separation must use the energy-conjugate onset traction");
    require(f.value().damageFactor()==0. && close(f.value().force().x,1.4),"The first crossing stores onset without instantaneous stiffness loss");
    require(close(f.value().getEffectiveFractureEnergy(),.005),"Equal mode energies must remain unchanged");
    f.pose(0,.002); f.evaluate();
    require(f.value().damageFactor()>0. && f.value().damageFactor()<1.,"A later larger sample must advance softening");
    require(close(f.value().getOnsetEffectiveDisplacement(),.0014,1.e-10),"The onset separation must stay frozen");
    // Restart from an undamaged fixture so the trapezoidal work includes the entire curve.
    Fixture work;
    Real integral=0.,lastForce=0.,lastDisplacement=0.;
    for(int i=0;i<=12000;++i) {
        const Real displacement=.012*i/12000.;
        work.pose(0,displacement); work.evaluate();
        const Real force=work.value().force().x;
        integral+=.5*(lastForce+force)*(displacement-lastDisplacement);
        lastForce=force; lastDisplacement=displacement;
    }
    require(close(integral,.005,1.e-8),"Monotonic tensile external work must equal A Gf after complete failure");
}
void testUnloadingAndCompressionRetention() {
    Fixture f;
    f.pose(0,.0014); f.evaluate();
    f.pose(0,.004); f.evaluate();
    const Real damage=f.value().damageFactor(), history=f.value().getMaximumEffectiveDisplacement();
    require(damage>0. && damage<1.,"Prepare a partially damaged bond");
    f.pose(0,.0003); f.evaluate();
    require(close(f.value().damageFactor(),damage) && close(f.value().getMaximumEffectiveDisplacement(),history),"Unloading must preserve damage and maximum effective separation");
    require(close(f.value().force().x,(1.-damage)*1000.*.0003),"Unloading stiffness must retain the original global factor");
    f.pose(0,-.0002); f.evaluate();
    require(close(f.value().damageFactor(),damage),"Pure compression must not advance damage");
    require(close(f.value().force().x,-.2),"Partial tensile damage must not weaken axial compression");
    require(close(f.value().normalElasticEnergy(),.5*1000.*.0002*.0002),"The stored compressive energy must remain undamaged");
    for(int mode=1;mode<4;++mode) {
        f.pose(mode,1.e-5); f.evaluate();
        Fixture elastic;
        elastic.value().setCrossSectionArea(0.); elastic.pose(mode,1.e-5); elastic.evaluate();
        require(math::norm(f.value().force()-(1.-damage)*elastic.value().force())<1.e-9,"All force branches must scale globally");
        require(math::norm(f.value().masterTorque()-(1.-damage)*elastic.value().masterTorque())<1.e-9,"All moment branches must scale globally");
        require(math::norm(f.value().torque()-(1.-damage)*elastic.value().torque())<1.e-9,"Stored native torsion-bending torque must scale globally");
    }
}
void testRigidRotation() {
    Fixture base, rotated;
    base.pose(0,.002); base.particles.host()[1].setOrientation(Quaternion::fromUnitAxisAngle(Vec3::unitX(),.003)); base.evaluate();
    rotated.pose(0,.002); rotated.particles.host()[1].setOrientation(Quaternion::fromUnitAxisAngle(Vec3::unitX(),.003));
    const Quaternion rotation=Quaternion::fromUnitAxisAngle(math::normalizedOrZero(Vec3{1.,2.,3.}),.8);
    const Vec3 translation{3.,-2.,5.};
    for(auto& p:rotated.particles.host()) {
        p.setPosition(translation+math::rotateUnit(rotation,p.position()));
        p.setOrientation(rotation*p.orientation());
    }
    rotated.evaluate();
    require(close(base.value().damageFactor(),rotated.value().damageFactor(),1.e-8),"Damage must be rigid-rotation objective");
    require(math::norm(rotated.value().force()-math::rotateUnit(rotation,base.value().force()))<1.e-8,"Force must rotate objectively");
    require(close(base.value().getOnsetEffectiveDisplacement(),rotated.value().getOnsetEffectiveDisplacement(),1.e-8),"Onset separation must be objective");
    require(math::norm(rotated.value().torque()-math::rotateUnit(rotation,base.value().torque()))<1.e-8,"Stored native torque must rotate objectively");
}

void testStrengthEnvelopeAndBudget() {
    Fixture disabled;
    disabled.value().setCrossSectionArea(0.);
    disabled.pose(0,.02); disabled.evaluate();
    require(!disabled.value().getDamageStarted() && disabled.value().damageFactor()==0. && close(disabled.value().force().x,20.),"Zero area must retain elastic response and disable fracture");

    require(close(execution::calculateBKCriticalEnergy(1.,3.,.004,.012,2.),.0085),"Unequal pure-mode energies and eta must follow standard BK");
    require(close(execution::calculateBKCriticalEnergy(1.,3.,.004,.012,1.),.010),"BK exponent must affect mixed-mode critical energy");
    const Vec3 axis=-Vec3::unitX(), zero=Vec3::zero();
    Real sigma=0.,bending=0.,tau=0.;
    const Real tensile=execution::calculateBondStrengthCriterion(sigma,bending,tau,{1.,0.,0.},zero,axis,1.,1.,1.,.5);
    require(close(sigma,1.)&&close(tau,0.),"Criterion must expose its normal and shear stresses");
    const Real compressedShear=execution::calculateBondStrengthCriterion(sigma,bending,tau,{-1.,1.,0.},zero,axis,1.,1.,1.,.5);
    const Real unconfinedShear=execution::calculateBondStrengthCriterion(sigma,bending,tau,{0.,1.,0.},zero,axis,1.,1.,1.,.5);
    require(close(tensile,1.),"Unit tensile stress must reach the tensile strength");
    require(close(compressedShear,1./2.25) && close(unconfinedShear,1.),"Compression must increase Mohr-Coulomb shear capacity");
    Fixture budget;
    budget.value().setModeICriticalEnergy(.0001); budget.value().setModeIICriticalEnergy(.0001);
    budget.pose(0,.0011); budget.evaluate();
    require(budget.value().getDamageStarted(),"Low-budget case must still detect onset");
    require(close(budget.value().getOnsetEffectiveStress(),1.1),"Low-budget onset must record effective stress");
    require(close(budget.value().getOnsetEffectiveDisplacement(),.0011),"Low-budget onset must record effective separation");
    require(budget.value().damageFactor()==1.,"An onset elastic energy exceeding A Gc must break immediately instead of inventing a larger budget");
    require(budget.value().hasInsufficientFractureEnergy(),"The insufficient positive fracture budget must be explicitly diagnosable");
    require(budget.value().getFailureEffectiveDisplacement()<=budget.value().getOnsetEffectiveDisplacement(),"An insufficient budget must remain visible in the rupture-separation diagnostic");
    require(close(budget.value().getEffectiveFractureEnergy(),.0001,1.e-15),"An insufficient budget must not modify the prescribed critical fracture energy");
    budget.pose(0,.0012); budget.evaluate();
    require(budget.value().damageFactor()==1.,"A later sample beyond the minimum softening interval must fail");
}

void testBKScaleInvariance() {
    const Real GIc=.004,GIIc=.012,eta=1.75;
    const Real expected=GIc+(GIIc-GIc)*std::pow(.75,eta);
    const Real denormal=std::numeric_limits<Real>::denorm_min();
    const Real maximum=std::numeric_limits<Real>::max();
    for(Real scale : {denormal,std::numeric_limits<Real>::min(),1.e-200,1.,1.e200,maximum/4.}) {
        const Real critical=execution::calculateBKCriticalEnergy(scale,3.*scale,GIc,GIIc,eta);
        require(std::isfinite(critical)&&close(critical,expected,1.e-13),"BK mixity must be invariant under representable positive energy scaling, including subnormals");
        require(close(execution::calculateBKCriticalEnergy(scale,0.,GIc,GIIc,eta),GIc,1.e-13),"Tiny and large pure opening energies must select GIc");
        require(close(execution::calculateBKCriticalEnergy(0.,scale,GIc,GIIc,eta),GIIc,1.e-13),"Tiny and large pure shear energies must select GIIc");
    }
    const Real equalExpected=GIc+(GIIc-GIc)*std::pow(.5,eta);
    require(close(execution::calculateBKCriticalEnergy(maximum,maximum,GIc,GIIc,eta),equalExpected,1.e-13),"BK must not overflow when the sum of two finite energies exceeds the floating-point maximum");
    require(execution::calculateBKCriticalEnergy(0.,0.,GIc,GIIc,eta)==0.,"Exactly zero total modal energy must retain its explicit no-mixity result");
}

void testLowEnergyStrengthCrossing() {
    Real damage=0.,stress0=0.,delta0=0.,maximum=0.,failure=0.,energy=0.;int started=0;
    const Real stiffness=1.e15,force=1.1,displacement=force/stiffness;
    const Real elasticEnergy=.5*stiffness*displacement*displacement;
    const bool failed=execution::updateBondDamage(damage,started,stress0,delta0,maximum,failure,energy,
        {force,0.,0.},Vec3::zero(),-Vec3::unitX(),elasticEnergy,0.,0.,0.,
        1.,.004,.012,1.75,1.,1.,0.);
    require(elasticEnergy>0.&&elasticEnergy<math::defaultTolerance,"Prepare a positive modal energy below the unrelated length tolerance");
    require(!failed&&started==1&&damage==0.,"A valid strength crossing must start softening even at a tiny positive modal energy");
    require(close(stress0,force)&&delta0>0.&&close(energy,.004,1.e-13),"Low-energy initiation must retain the effective stress, positive separation and correct pure-mode threshold");

    Fixture native;
    require(native.value().setStiffness(stiffness,stiffness,0.,0.),"Configure a native high-stiffness bond");
    native.pose(0,1.2e-15);native.evaluate();
    require(native.value().getDamageStarted()&&native.value().damageFactor()==0.,"Native CPU evaluation must preserve tiny-energy strength initiation");
    require(native.value().getOnsetEffectiveStress()>1.&&native.value().normalElasticEnergy()>0.
        &&native.value().normalElasticEnergy()<math::defaultTolerance,"Native high-stiffness state must cross strength without losing its positive stored energy");
}

void testPureModeElasticLoads() {
    for(Real angle : {1.e-9,1.e-7,.02}) {
        for(int mode=1;mode<4;++mode) {
            Fixture f;f.value().setCrossSectionArea(0.);
            require(f.value().setStiffness(1000.,400.,200.,100.),"Configure pure-mode elastic reference stiffnesses");
            f.pose(mode,angle);f.evaluate();
            const bond& b=f.value();
            if(mode==1) require(close(math::norm(b.force()),400.*std::sin(angle),1.e-10),"Pure shear force must match the existing angular response");
            if(mode==2) require(close(b.masterTorque().z,200.*std::sin(.5*angle)+100.*std::sin(angle),1.e-10),"Pure bending endpoint torque must preserve the existing B-coefficient response");
            if(mode==3) require(close(b.torque().x,100.*std::sin(angle),1.e-10),"Pure torsional torque must preserve the existing angular response");
        }
    }
}

void testSquaredElasticExpressions() {
    const Vec3 masterAxis=math::normalizedOrZero(Vec3{1.,2.,3.});
    const Vec3 slaveAxis=math::normalizedOrZero(Vec3{-1.,4.,2.});
    for(int sample=1;sample<=8;++sample) {
        Fixture f;f.value().setCrossSectionArea(0.);
        require(f.value().setStiffness(1000.,400.,200.,100.),"Configure algebraic elastic equivalence stiffnesses");
        f.particles.host()[0].setPosition({.03*sample,-.02*sample,.01*sample});
        f.particles.host()[1].setPosition({.95+.01*sample,.1+.02*sample,-.07});
        f.particles.host()[0].setOrientation(Quaternion::fromUnitAxisAngle(masterAxis,.04*sample));
        f.particles.host()[1].setOrientation(Quaternion::fromUnitAxisAngle(slaveAxis,-.03*sample));
        f.evaluate();
        const auto& b=f.value();
        const auto& master=f.particles.host()[0];const auto& slave=f.particles.host()[1];
        const Vec3 a=math::rotateUnit(master.orientation(),b.masterEndpointLocalNormal());
        const Vec3 bNormal=math::rotateUnit(slave.orientation(),b.slaveEndpointLocalNormal());
        const Vec3 mt1=math::rotateUnit(master.orientation(),b.masterEndpointLocalTangent1());
        const Vec3 mt2=math::rotateUnit(master.orientation(),b.masterEndpointLocalTangent2());
        const Vec3 st1=math::rotateUnit(slave.orientation(),b.slaveEndpointLocalTangent1());
        const Vec3 st2=math::rotateUnit(slave.orientation(),b.slaveEndpointLocalTangent2());
        const Vec3 mp=master.position()+math::rotateUnit(master.orientation(),b.masterEndpointLocalPosition());
        const Vec3 sp=slave.position()+math::rotateUnit(slave.orientation(),b.slaveEndpointLocalPosition());
        const Real length=math::norm(mp-sp);const Vec3 direction=(sp-mp)/length;
        const Real expectedNormal=.5*b.coefficientB1()*(length-b.equivalentLength())*(length-b.equivalentLength());
        const Real expectedShear=.5*b.coefficientB2()*math::normSquared(.5*(a-bNormal)-direction);
        const Real expectedBending=.5*b.bendingStiffness()*math::normSquared(a+bNormal);
        const Real expectedTorsion=.25*b.coefficientB4()*(math::normSquared(mt1-st1)+math::normSquared(mt2-st2)-math::normSquared(a+bNormal));
        require(close(b.normalElasticEnergy(),expectedNormal,1.e-12)
            &&close(b.shearElasticEnergy(),expectedShear,1.e-12)&&close(b.bendingElasticEnergy(),expectedBending,1.e-12)
            &&close(b.torsionalElasticEnergy(),expectedTorsion,1.e-12),"Modal energies must follow the enabled squared-vector expressions using the stored raw bending stiffness");
        requireFiniteEnergies(b);
        const Vec3 expectedForce=b.coefficientB1()*(length-b.equivalentLength())*direction
            +.5*b.coefficientB2()/length*((bNormal-a)-math::dot(bNormal-a,direction)*direction);
        const Vec3 expectedTBT=b.coefficientB3()*math::cross(bNormal,a)
            -.5*b.coefficientB4()*(math::cross(st1,mt1)+math::cross(st2,mt2));
        const Vec3 expectedMaster=-.5*b.coefficientB2()*math::cross(direction,a)+expectedTBT
            +math::cross(mp-master.position(),expectedForce);
        const Vec3 expectedSlave=.5*b.coefficientB2()*math::cross(direction,bNormal)-expectedTBT
            +math::cross(sp-slave.position(),-expectedForce);
        require(math::norm(b.force()-expectedForce)<1.e-10&&math::norm(b.torque()-expectedTBT)<1.e-10
            &&math::norm(b.masterTorque()-expectedMaster)<1.e-10&&math::norm(b.slaveTorque()-expectedSlave)<1.e-10,"Elastic evaluation must retain the established forces and endpoint torques");
    }
}

void testSquaredModalResolution() {
    for(Real angle : {1.e-9,1.e-7,.02}) for(int mode=1;mode<4;++mode) {
        Fixture f;f.value().setCrossSectionArea(0.);
        require(f.value().setStiffness(1000.,400.,200.,100.),"Configure squared-energy small-angle resolution");
        f.pose(mode,angle);f.evaluate();requireFiniteEnergies(f.value());
        const Real modal[]{f.value().normalElasticEnergy(),f.value().shearElasticEnergy(),f.value().bendingElasticEnergy(),f.value().torsionalElasticEnergy()};
        const Real stiffness[]{1000.,400.,200.,100.};
        const Real halfSine=std::sin(.5*angle),expected=2.*stiffness[mode]*halfSine*halfSine;
        require(modal[mode]>0.&&std::abs(modal[mode]-expected)<1.e-10*expected,"Squared-vector angular energies must resolve positive small-angle energy without subtracting nearly equal dot products");
    }
    for(int mode=0;mode<4;++mode) {
        Fixture f;f.value().setCrossSectionArea(0.);
        require(f.value().setStiffness(1.e-20,1.e-20,2.e-20,1.e-20),"Configure tiny positive modal energies");
        f.pose(mode,.02);f.evaluate();requireFiniteEnergies(f.value());
        const Real modal[]{f.value().normalElasticEnergy(),f.value().shearElasticEnergy(),f.value().bendingElasticEnergy(),f.value().torsionalElasticEnergy()};
        require(modal[mode]>0.&&modal[mode]<math::defaultTolerance,"Each positive modal energy below tolerance must remain positive");
    }
    Fixture tinyBending;tinyBending.value().setCrossSectionArea(0.);
    require(tinyBending.value().setStiffness(1000.,400.,1.e-30,100.),"Configure raw bending stiffness below B-coefficient reconstruction precision");
    tinyBending.pose(2,1.e-4);tinyBending.evaluate();
    const Real halfSine=std::sin(.5e-4),expected=2.e-30*halfSine*halfSine;
    require(tinyBending.value().bendingElasticEnergy()>0.
        &&std::abs(tinyBending.value().bendingElasticEnergy()-expected)<1.e-10*expected,"Squared bending energy must use the original raw stiffness even when its B-coefficient contribution cancels numerically");
    Fixture stressFree;stressFree.evaluate();requireFiniteEnergies(stressFree.value());
}

void testBrittleFracture() {
    Fixture f;
    require(f.value().setModeICriticalEnergy(0.) && f.value().setModeIICriticalEnergy(0.),"Configure strength-only brittle fracture");
    for(Real opening : {-.002,.0009}) {
        f.pose(0,opening); f.evaluate();
        require(f.value().damageFactor()==0. && !f.value().getDamageStarted(),"Compression and substrength tension must not break a zero-energy bond");
        require(close(f.value().force().x,1000.*opening),"An intact brittle bond must retain its elastic force");
    }
    f.pose(0,.0011); f.evaluate();
    require(f.value().damageFactor()==1.,"A zero-energy bond must fail at its first sampled strength crossing");
    require(f.value().getDamageStarted(),"Brittle failure must record that the strength threshold was reached");
    require(math::norm(f.value().normal())==0. && math::norm(f.value().force())==0. && math::norm(f.value().masterTorque())==0.
        && math::norm(f.value().slaveTorque())==0. && math::norm(f.value().torque())==0. && total(f.value())==0.,"Brittle failure must clear every load and stored energy in the same evaluation");
    f.pose(0,.002); f.evaluate();
    require(f.value().damageFactor()==1. && math::norm(f.value().force())==0. && total(f.value())==0.,"Repeated evaluation cannot reactivate a broken brittle bond");

    Real d=0.,s0=0.,d0=0.,maximum=0.,failure=0.,energy=0.;int started=0;
    const auto advance=[&](Real force) {
        return execution::updateBondDamage(d,started,s0,d0,maximum,failure,energy,
            {force,0.,0.},Vec3::zero(),-Vec3::unitX(),0.,0.,0.,0.,
            1.,0.,0.,1.75,1.,1.,0.);
    };
    require(!advance(.999) && d==0.,"Q below one must not fail even when separation stiffnesses are zero");
    require(advance(1.) && d==1.,"Q equal to one must report brittle failure without computing separation");
    require(started==1,"The shared brittle helper must set the onset flag at failure");
    require(advance(0.) && d==1.,"An already failed helper must report failure again on unloading");
    d=.25; started=1;
    require(!advance(.999) && d==.25,"Switching a started bond to zero energies below strength must preserve damage");
    require(advance(1.) && d==1.,"Switching a started bond to zero energies at strength must break it immediately");

    for(int mode : {0,3}) {
        Fixture softening;
        require(softening.value().setModeICriticalEnergy(mode==0?.005:0.)
            && softening.value().setModeIICriticalEnergy(mode==3?.005:0.),"Configure one positive fracture energy");
        const Real onset=mode==0?.0014:.003;
        softening.pose(mode,onset); softening.evaluate();
        require(softening.value().getDamageStarted() && softening.value().damageFactor()==0.,"One positive energy must retain first-crossing softening initialization");
        softening.pose(mode,1.4*onset); softening.evaluate();
        require(softening.value().damageFactor()>0. && softening.value().damageFactor()<1.,"The matching positive mode energy must retain a finite softening interval");
    }
}

#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
void testSquaredNativeGPUParity() {
    const Vec3 axis=math::normalizedOrZero(Vec3{1.,2.,3.});
    for(int sample=0;sample<=8;++sample) {
        Fixture cpu,gpu;
        for(Fixture* f : {&cpu,&gpu}) {
            f->value().setCrossSectionArea(0.);
            require(f->value().setStiffness(1000.,400.,200.,100.),"Configure squared CPU/GPU modal energy parity");
            if(sample>0) {
                f->particles.host()[0].setPosition({.03*sample,-.02*sample,.01*sample});
                f->particles.host()[1].setPosition({.95+.01*sample,.1+.02*sample,-.07});
                f->particles.host()[0].setOrientation(Quaternion::fromUnitAxisAngle(axis,.04*sample));
                f->particles.host()[1].setOrientation(Quaternion::fromUnitAxisAngle(axis,-.03*sample));
            }
        }
        gpu.particles.setUseDevice(true);gpu.bonds.setUseDevice(true);
        cpu.evaluate();gpu.evaluateGPU();
        const bond& a=cpu.value();const bond& b=gpu.value();
        requireFiniteEnergies(a);requireFiniteEnergies(b);
        require(close(a.normalElasticEnergy(),b.normalElasticEnergy(),1.e-10)&&close(a.shearElasticEnergy(),b.shearElasticEnergy(),1.e-10)
            &&close(a.bendingElasticEnergy(),b.bendingElasticEnergy(),1.e-10)&&close(a.torsionalElasticEnergy(),b.torsionalElasticEnergy(),1.e-10),"Native CPU/GPU modal energies must agree at stress-free and general deformed poses");
        require(math::norm(a.force()-b.force())<1.e-9&&math::norm(a.torque()-b.torque())<1.e-9
            &&math::norm(a.masterTorque()-b.masterTorque())<1.e-9&&math::norm(a.slaveTorque()-b.slaveTorque())<1.e-9,"Squared modal energies must preserve native CPU/GPU force and torque parity");
    }
}

void testBrittleNativeGPUParity() {
    const Real crossing[]{.0011,.003,.003,.003};
    for(int mode=0;mode<4;++mode) for(bool zeroSeparationStiffness : {false,true}) {
        if(zeroSeparationStiffness && mode<2) continue;
        Fixture cpu,gpu;
        for(Fixture* f : {&cpu,&gpu}) {
            f->value().setModeICriticalEnergy(0.); f->value().setModeIICriticalEnergy(0.);
            f->value().setStiffness(zeroSeparationStiffness?0.:1000.,zeroSeparationStiffness?0.:400.,mode==2?200.:100.,100.);
            require(f->value().isValid(),"Brittle bending and twist must accept zero axial/shear stiffness");
        }
        gpu.particles.setUseDevice(true); gpu.bonds.setUseDevice(true);
        for(Real factor : {.1,1.,1.5}) {
            cpu.pose(mode,factor*crossing[mode]); gpu.pose(mode,factor*crossing[mode]);
            for(Fixture* f : {&cpu,&gpu}) for(auto& p : f->particles.host()) { p.setForce(Vec3::zero()); p.setTorque(Vec3::zero()); }
            cpu.evaluate(); gpu.evaluateGPU();
            const bond& expected=cpu.value(); const bond& actual=gpu.value();
            require(actual.damageFactor()==expected.damageFactor() && actual.damageFactor()==(factor<1.?0.:1.),"Native CPU/GPU must agree on intact, first brittle failure and repeated failure states");
            require(actual.getDamageStarted()==expected.getDamageStarted()&&actual.getDamageStarted()==(factor>=1.),"Native CPU/GPU brittle onset flags must agree before and after the strength crossing");
            require(math::norm(actual.force()-expected.force())<1.e-10 && math::norm(actual.masterTorque()-expected.masterTorque())<1.e-10
                && math::norm(actual.slaveTorque()-expected.slaveTorque())<1.e-10 && math::norm(actual.torque()-expected.torque())<1.e-10,"Native CPU/GPU brittle loads must agree in every mode");
            if(factor>=1.) require(math::norm(actual.normal())==0. && math::norm(actual.force())==0. && math::norm(actual.masterTorque())==0.
                && math::norm(actual.slaveTorque())==0. && math::norm(actual.torque())==0. && total(actual)==0.,"Native GPU brittle failure must clear normal, all loads and all energies immediately");
        }
    }
}

void testLowEnergyNativeGPUParity() {
    Fixture cpu,gpu;
    for(Fixture* f : {&cpu,&gpu}) {
        require(f->value().setStiffness(1.e15,1.e15,0.,0.),"Configure native CPU/GPU tiny-energy stiffnesses");
        f->pose(0,1.2e-15);
    }
    gpu.particles.setUseDevice(true);gpu.bonds.setUseDevice(true);
    cpu.evaluate();gpu.evaluateGPU();
    require(cpu.value().getDamageStarted()&&gpu.value().getDamageStarted(),"Native CPU/GPU must initiate damage at a high-stiffness tiny-energy strength crossing");
    require(cpu.value().damageFactor()==0.&&gpu.value().damageFactor()==0.,"Tiny-energy first crossing must retain the established no-instantaneous-softening rule");
    require(close(cpu.value().getOnsetEffectiveStress(),gpu.value().getOnsetEffectiveStress(),1.e-10)
        &&close(cpu.value().getEffectiveFractureEnergy(),gpu.value().getEffectiveFractureEnergy(),1.e-13),"Native CPU/GPU low-energy onset stress and BK threshold must agree");
}
#endif

void testExactNativeTorque() {
    for(int mode=0;mode<4;++mode) {
        Fixture f;
        f.value().setCrossSectionArea(0.);
        f.pose(mode,.02); f.evaluate();
        const bond& b=f.value();
        const Quaternion qm=f.particles.host()[0].orientation(), qs=f.particles.host()[1].orientation();
        const Vec3 nm=math::rotateUnit(qm,b.masterEndpointLocalNormal()), ns=math::rotateUnit(qs,b.slaveEndpointLocalNormal());
        const Vec3 tm1=math::rotateUnit(qm,b.masterEndpointLocalTangent1()), ts1=math::rotateUnit(qs,b.slaveEndpointLocalTangent1());
        const Vec3 tm2=math::rotateUnit(qm,b.masterEndpointLocalTangent2()), ts2=math::rotateUnit(qs,b.slaveEndpointLocalTangent2());
        const Vec3 expected=b.coefficientB3()*math::cross(ns,nm)
            -.5*b.coefficientB4()*(math::cross(ts1,tm1)+math::cross(ts2,tm2));
        require(math::norm(b.torque()-expected)<1.e-12,"Stored bond torque must equal the exact existing TBT expression");
    }
    Fixture shear;
    shear.value().setCrossSectionArea(0.);
    shear.pose(1,.02); shear.evaluate();
    require(math::norm(shear.value().torque())==0.,"Parallel endpoint frames under shear must have zero native TBT");
    require(math::norm(shear.value().masterTorque())>1.,"Zero native TBT must not be replaced by nonzero endpoint couples");
    Fixture degenerateBend;
    degenerateBend.pose(2,.1); degenerateBend.evaluate();
    require(math::norm(degenerateBend.value().torque())<1.e-12,"kb=ks L^2/4 must retain its native TBT cancellation");
    require(math::norm(degenerateBend.value().masterTorque())>1.,"The cancellation case must retain nonzero endpoint load couples");
    require(!degenerateBend.value().getDamageStarted(),"Strength must use native TBT, not a translated or averaged endpoint couple");
}

void testEnergyEquivalentHistory() {
    Real damage=0.,stress0=0.,delta0=0.,maximum=0.,failure=0.,energy=0.;
    int started=0;
    const auto advance=[&](const Vec3& force,Real un,Real us) {
        return execution::updateBondDamage(damage,started,stress0,delta0,maximum,failure,energy,
            force,Vec3::zero(),-Vec3::unitX(),un,us,0.,0.,1.,.005,.005,1.75,1.,1.,0.);
    };
    require(!advance({1.,0.,0.},.0005,0.),"First aligned onset must not fail");
    require(close(stress0,1.)&&close(delta0,.001)&&close(maximum,.001)&&close(failure,.01)&&damage==0.,"Pure opening must retain its original conjugate traction-separation law");
    advance({2.,0.,0.},.002,0.);
    require(close(damage,5./9.),"Pure-opening linear softening must equal 5/9");
    advance({1.5,0.,0.},.001125,0.);
    require(close(damage,5./9.)&&close(maximum,.002),"Unloading cannot reduce the maximum equivalent separation or damage");
    damage=stress0=delta0=maximum=failure=energy=0.;started=0;
    advance({0.,1.,0.},0.,.00125);
    const Real shearOnset=2.*.00125;
    require(close(stress0,1.)&&close(delta0,shearOnset)&&close(failure,.01)&&damage==0.,"Pure shear must recover its energy-over-section-stress separation and frozen stress-based failure length");
    advance({0.,2.,0.},0.,.005);
    require(close(damage,2./3.),"Pure-shear Eq13 value must equal 2/3");
    advance({0.,1.2,0.},0.,.0018);
    require(close(damage,2./3.)&&close(maximum,2.*shearOnset),"Shear unloading must be irreversible");
    damage=stress0=delta0=maximum=failure=energy=0.;started=0;
    execution::updateBondDamage(damage,started,stress0,delta0,maximum,failure,energy,
        {-3.,4.,0.},Vec3::zero(),-Vec3::unitX(),1.,1.,0.,0.,1.,100.,100.,1.75,1.,1.,0.);
    require(close(stress0,4.)&&close(delta0,.5)&&close(failure,50.),"Compression must be excluded; shear separation is 2Us/(A tau) and failure uses the onset shear stress");
}
void testStiffnessPreconditions() {
    for(Real GI : {0.,.005}) for(Real GII : {0.,.005}) {
        Fixture f;
        require(f.value().setModeICriticalEnergy(GI) && f.value().setModeIICriticalEnergy(GII),"Configure fracture energy validation matrix");
        require(f.value().isValid(),"Positive strengths and stiffnesses must support brittle and softening bonds");
        f.value().setTensileStrength(0.);
        require(!f.value().isValid(),"Every positive-area bond requires a positive tensile strength");
        f.value().setTensileStrength(1.); f.value().setCohesion(0.);
        require(!f.value().isValid(),"Every positive-area bond requires a positive cohesion");
        f.value().setCohesion(1.);
        const bool brittle=GI==0. && GII==0.;
        f.value().setStiffness(0.,400.,100.,100.);
        require(f.value().isValid()==brittle,"Positive normal stiffness is required only when softening recovers separation");
        f.value().setStiffness(1000.,0.,100.,100.);
        require(f.value().isValid()==brittle,"Euclidean softening requires positive normal and shear stiffnesses");
        f.value().setCrossSectionArea(0.); f.value().setTensileStrength(0.); f.value().setCohesion(0.);
        f.value().setStiffness(0.,0.,100.,100.);
        require(f.value().isValid(),"Zero area must still disable fracture without requiring strengths or separation stiffnesses");
    }
    Fixture torsion;
    require(torsion.value().setStiffness(1000.,0.,100.,100.)&&!torsion.value().isValid(),"Positive-toughness torsion cannot recover a shear-group separation from zero shear stiffness");
}

void testFourEnergySeparation() {
    const auto sample=[&](Real axialStress,Real un,Real us,Real ub,Real ut) {
        Real d=0.,s0=0.,d0=0.,maximum=0.,failure=0.,energy=0.;int started=0;
        // A=1 gives R/I=4*sqrt(pi); the prescribed moment produces sigma_b=5.
        const Vec3 torque{0.,0.,5./(4.*std::sqrt(math::pi))};
        execution::updateBondDamage(d,started,s0,d0,maximum,failure,energy,
            {axialStress,3.,0.},torque,-Vec3::unitX(),un,us,ub,ut,
            1.,100.,100.,1.75,1.,1.,0.);
        require(started==1&&d==0.,"Four-energy first crossing must store history with D=0");
        return d0;
    };
    const Real un=.5*2.*.1*.1,us=.5*3.*.2*.2;
    const Real ub=.5*5.*.3*.3,ut=.5*7.*.4*.4;
    require(close(sample(2.,un,us,ub,ut),std::hypot(2.*(un+ub)/7.,2.*(us+ut)/3.)),"Grouped energies must use opening/shear section stresses before an unweighted Euclidean norm");
    require(close(sample(-2.,un,0.,ub,0.),2.*ub/5.),"Axial compression cannot cancel positive bending energy");
    require(close(sample(-100.,100.,0.,ub,0.),2.*ub/5.),"Arbitrarily large compression must remain outside the damage-driving energy");
    require(close(sample(2.,un,us,0.,0.),std::hypot(2.*un/7.,2.*us/3.)),"The component separation must not divide by absent rotational stiffnesses");
    require(close(sample(2.,un,-1.e-15,-1.e-15,-1.e-15),2.*un/7.),"Negative roundoff energies must be excluded from the damage-driving components");
}
void testRawStiffnessStorage() {
    Fixture f;const bond before=f.value();
    require(before.normalStiffness()==1000.&&before.shearStiffness()==400.&&before.bendingStiffness()==100.&&before.torsionalStiffness()==100.,"Store the original four stiffness inputs exactly");
    require(!f.value().setStiffness(2.,-1.,3.,4.),"Reject negative stiffness transactionally");
    require(f.value().normalStiffness()==before.normalStiffness()&&f.value().shearStiffness()==before.shearStiffness()
        &&f.value().bendingStiffness()==before.bendingStiffness()&&f.value().torsionalStiffness()==before.torsionalStiffness(),"Rejected inputs cannot partially replace raw stiffnesses");
    const bond copied=before;bond assigned;assigned=before;
    require(copied.normalStiffness()==1000.&&copied.shearStiffness()==400.&&copied.bendingStiffness()==100.&&copied.torsionalStiffness()==100.
        &&assigned.normalStiffness()==1000.&&assigned.shearStiffness()==400.&&assigned.bendingStiffness()==100.&&assigned.torsionalStiffness()==100.,"Copy and assignment must preserve all four raw inputs");
    bond large{1.e200};
    require(large.setStiffness(2.,1.e-300,3.,4.),"Prepare finite coefficients at a large length");
    require(!large.setStiffness(7.,1.,8.,9.),"Reject derived B2 overflow");
    require(large.normalStiffness()==2.&&large.shearStiffness()==1.e-300&&large.bendingStiffness()==3.&&large.torsionalStiffness()==4.,"Derived-overflow rejection must preserve raw inputs");
    require(f.value().setEquivalentLength(2.),"Change reference length");
    require(f.value().normalStiffness()==0.&&f.value().shearStiffness()==0.&&f.value().bendingStiffness()==0.&&f.value().torsionalStiffness()==0.,"Length invalidation must clear all raw stiffnesses");
    Fixture tiny;tiny.value().setCrossSectionArea(0.);
    require(tiny.value().setStiffness(1000.,400.,1.e-30,100.),"Retain a representable tiny raw bending stiffness");
    require(tiny.value().bendingStiffness()==1.e-30,"Raw stiffness storage must preserve a tiny original input independently of the existing energy formulas");
}
void testCompressionCannotCancelBending() {
    Fixture f;require(f.value().setStiffness(1.,400.,2000.,100.),"Configure a compression/bending edge");
    f.pose(2,1.e-4);f.particles.host()[1].setPosition({.999,0.,0.});f.evaluate();
    const auto& b=f.value();
    require(b.getDamageStarted()&&b.getOnsetEffectiveStress()>1.,"Section stress must cross strength despite compression");
    require(b.getOnsetEffectiveDisplacement()>0.&&b.damageFactor()==0.,"Bending-driven strength onset must have a finite energy coordinate even under axial compression");
    const Real frozenEnergy=b.getEffectiveFractureEnergy();
    f.pose(2,1.e-4);f.evaluate();
    require(f.value().damageFactor()<1.&&close(f.value().getEffectiveFractureEnergy(),frozenEnergy),"Removing compression cannot produce artificial zero-onset rupture or change frozen Gc");
}

void testFrozenBKMixity() {
    Fixture lowEta, highEta;
    for(Fixture* f : {&lowEta,&highEta}) {
        f->value().setModeICriticalEnergy(.004);
        f->value().setModeIICriticalEnergy(.012);
        f->pose(0,.0012);
        f->particles.host()[1].setOrientation(Quaternion::fromUnitAxisAngle(Vec3::unitX(),.002));
    }
    lowEta.value().setModeMixityExponent(1.);
    highEta.value().setModeMixityExponent(2.);
    lowEta.evaluate(); highEta.evaluate();
    require(lowEta.value().getDamageStarted()&&highEta.value().getDamageStarted(),"Mixed loading must initiate both bonds");
    const Real un=.5*1000.*.0012*.0012, ut=100.*(1.-std::cos(.002));
    const Real fraction=ut/(un+ut);
    const Real expectedLow=.004+.008*fraction;
    const Real expectedHigh=.004+.008*fraction*fraction;
    require(close(lowEta.value().getEffectiveFractureEnergy(),expectedLow,1.e-9),"Mixed initiation must use standard BK mode-II fraction");
    require(close(highEta.value().getEffectiveFractureEnergy(),expectedHigh,1.e-9),"Stored onset energy must retain eta dependence");
    const Real frozen=lowEta.value().getEffectiveFractureEnergy();
    lowEta.pose(0,.003); lowEta.evaluate();
    require(close(lowEta.value().getEffectiveFractureEnergy(),frozen),"Changing the loading path must not recompute onset Gc");
    Fixture tensile, twist;
    for(Fixture* f : {&tensile,&twist}) { f->value().setModeICriticalEnergy(.004); f->value().setModeIICriticalEnergy(.012); }
    tensile.pose(0,.0015); tensile.evaluate();
    twist.pose(3,.005); twist.evaluate();
    require(close(tensile.value().getEffectiveFractureEnergy(),.004),"Pure tension must select GIc");
    require(close(twist.value().getEffectiveFractureEnergy(),.012),"Pure twist must select GIIc");
}

void testEuclideanMixedCoordinate() {
    for(Real area : {1.,1.5,2.}) {
        Real damage=0.,stress0=0.,onset=0.,maximum=0.,failure=0.,energy=0.;int started=0;
        const auto advance=[&](Real scale) {
            return execution::updateBondDamage(damage,started,stress0,onset,maximum,failure,energy,
                {2.,1.,0.},Vec3::zero(),-Vec3::unitX(),
                .0002*scale,.0003*scale,.0001*scale,.0004*scale,
                area,.004,.012,1.75,1.,1.,0.);
        };
        require(!advance(1.)&&started&&damage==0.,"Prepare mixed-mode onset for the Euclidean-coordinate check");
        const Real expected=std::hypot(.0003,.0014);
        const Real onsetStress=std::sqrt(5.)/area;
        require(close(onset,expected,1.e-12),"Normal and shear component separations must use an ordinary unweighted norm");
        require(close(stress0,onsetStress,1.e-12),"The frozen onset stress must follow the independent circular-section force projection");
        require(close(failure,2.*energy/onsetStress,1.e-12),"The failure separation must use the onset opening/shear stress norm");
        advance(3.);
        const Real expectedDamage=(1.-1./3.)/(1.-expected/failure);
        require(close(damage,expectedDamage,1.e-12),"At fixed section stress, tripled energy triples the separation and preserves irreversible softening");
    }
}

void testZeroEnergyStrengthCrossing() {
    Real damage=0.,stress0=0.,onset=0.,maximum=0.,failure=0.,energy=0.;int started=0;
    require(execution::updateBondDamage(damage,started,stress0,onset,maximum,failure,energy,
        {2.,0.,0.},Vec3::zero(),-Vec3::unitX(),0.,0.,0.,0.,
        1.,.004,.012,1.75,1.,1.,0.),"A strength crossing with exactly zero driving energy must fail safely");
    require(started&&damage==1.&&std::isfinite(onset)&&std::isfinite(failure),"A zero-energy crossing cannot create NaN history or a spurious finite softening branch");
}

void setEnergyWorkPose(Fixture& f,int mode,Real a) {
    if(mode<4) { f.pose(mode,a);return; }
    f.pose(0,a);
    if(mode==4) f.particles.host()[1].setOrientation(Quaternion::fromUnitAxisAngle(Vec3::unitX(),2.*a));
    else {
        const Real c=std::sqrt(.5*f.value().normalStiffness()/(4.*f.value().coefficientB2()));
        const Real angle=mode==6?2.*std::asin(c*a):2.*a;
        f.particles.host()[1].setPosition({(1.+a)*std::cos(angle),(1.+a)*std::sin(angle),0.});
    }
}

Real generalizedExternalForce(const Fixture& f,int mode,Real a) {
    const bond& b=f.bonds.host().front();
    if(mode==0) return b.force().x;
    if(mode==1) return math::dot(b.force(),Vec3{-std::sin(a),std::cos(a),0.});
    if(mode==2) return .5*(b.masterTorque().z-b.slaveTorque().z);
    if(mode==3) return -b.slaveTorque().x;
    if(mode==4) return b.force().x-2.*b.slaveTorque().x;
    const Real c=std::sqrt(.5*b.normalStiffness()/(4.*b.coefficientB2()));
    const Real angle=mode==6?2.*std::asin(c*a):2.*a;
    const Real angularRate=mode==6?2.*c/std::sqrt(1.-c*c*a*a):2.;
    const Vec3 derivative{std::cos(angle)-(1.+a)*std::sin(angle)*angularRate,
        std::sin(angle)+(1.+a)*std::cos(angle)*angularRate,0.};
    return math::dot(b.force(),derivative);
}

struct WorkMeasurement { Real work{},critical{},dissipated{}; };
WorkMeasurement measureFractureWork(int mode,int steps) {
    Fixture f,elastic;
    require(f.value().setStiffness(1000.,400.,200.,100.),"Configure positive unequal physical stiffnesses for external-work checks");
    require(elastic.value().setStiffness(1000.,400.,200.,100.)&&elastic.value().setCrossSectionArea(0.),"Configure independent undamaged native work reference");
    require(f.value().setModeICriticalEnergy(.004)&&f.value().setModeIICriticalEnergy(.008),"Configure unequal pure-mode fracture budgets");
    const Real maximum[]{.025,.16,.16,.16,.16,.16,.16};
    Real work=0.,previousForce=0.,previousAmplitude=0.,frozen=0.,previousTrialEnergy=0.,previousDamage=0.,dissipated=0.;
    bool partial=false;
    for(int i=0;i<=steps;++i) {
        const Real amplitude=maximum[mode]*Real(i)/steps;
        setEnergyWorkPose(f,mode,amplitude);setEnergyWorkPose(elastic,mode,amplitude);
        f.evaluate();elastic.evaluate();
        const auto& b=f.value();
        const Real force=generalizedExternalForce(f,mode,amplitude);
        work+=.5*(previousForce+force)*(amplitude-previousAmplitude);
        const Real trialEnergy=total(elastic.value());
        dissipated+=.5*(previousTrialEnergy+trialEnergy)*(b.damageFactor()-previousDamage);
        previousTrialEnergy=trialEnergy;previousDamage=b.damageFactor();
        previousForce=force;previousAmplitude=amplitude;
        if(b.getDamageStarted()) {
            if(frozen==0.) frozen=b.getEffectiveFractureEnergy();
            require(close(b.getEffectiveFractureEnergy(),frozen,1.e-12),"The physical-work path cannot change Gc after onset");
        }
        partial|=b.damageFactor()>0.&&b.damageFactor()<1.;
    }
    require(partial&&f.value().damageFactor()==1.,"Every pure and proportional mixed work path must finish its softening interval");
    require(total(f.value())==0.,"Complete opening/shearing rupture must release all remaining stored energy");
    return {work,f.value().crossSectionArea()*frozen,dissipated};
}

void testNativeFractureWorkBudget() {
    for(int mode=0;mode<7;++mode) {
        // The independent Ua*dD trapezoidal oracle needs this density on the
        // largest-amplitude mixed path; keep every path at the same resolution.
        const auto coarse=measureFractureWork(mode,8000);
        const auto fine=measureFractureWork(mode,16000);
        const Real coarseError=std::abs(coarse.work-coarse.dissipated);
        const Real fineError=std::abs(fine.work-fine.dissipated);
        if(std::abs(fine.work-fine.dissipated)>=2.e-7||fineError>=2.e-7||fineError>coarseError+1.e-8)
            std::cerr<<std::setprecision(17)<<"Fracture work mode "<<mode
                <<": coarse_work="<<coarse.work<<", coarse_dissipated="<<coarse.dissipated
                <<", coarse_balance_error="<<std::abs(coarse.work-coarse.dissipated)
                <<", fine_work="<<fine.work<<", fine_dissipated="<<fine.dissipated
                <<", fine_balance_error="<<std::abs(fine.work-fine.dissipated)
                <<", A_Gc="<<fine.critical<<", fine_budget_deviation="<<fine.work-fine.critical<<'\n';
        require(fineError<2.e-7,"Native paths must satisfy actual energy balance; equality to A Gc is not a current mixed-mode guarantee");
        require(fineError<=coarseError+1.e-8,"Halving the prescribed loading increment cannot worsen fracture-work error beyond quadrature noise");
    }
}

void testFixedDirectionTwoSpringWork() {
    // Analytic damage-law verification, not a native finite-angle beam experiment.
    constexpr Real kn=1000.,ks=400.,normalDirection=.8,shearDirection=.6;
    const Real stiffness=kn*normalDirection*normalDirection+ks*shearDirection*shearDirection;
    Real damage=0.,stress0=0.,onset=0.,maximum=0.,failure=0.,energy=0.,work=0.,previousQ=0.,previousTraction=0.;int started=0;
    for(int i=0;i<=12000;++i) {
        const Real q=.03*Real(i)/12000.;
        const Real dn=normalDirection*q,ds=shearDirection*q;
        execution::updateBondDamage(damage,started,stress0,onset,maximum,failure,energy,
            {kn*dn,ks*ds,0.},Vec3::zero(),-Vec3::unitX(),.5*kn*dn*dn,.5*ks*ds*ds,0.,0.,
            1.,.004,.012,1.75,1.,1.,0.);
        const Real traction=(1.-damage)*stiffness*q;
        work+=.5*(previousTraction+traction)*(q-previousQ);previousQ=q;previousTraction=traction;
    }
    require(started&&damage==1.,"The fixed-direction two-spring path must finish fracture");
    const Real stressPerDisplacement=std::hypot(kn*normalDirection,ks*shearDirection);
    require(close(stress0,stressPerDisplacement*onset,1.e-12),"The two-spring onset must use the independently calculated section-stress norm");
    require(close(failure,2.*energy/stress0,1.e-12),"The two-spring failure separation must use the frozen section stress");
    const Real expectedWork=.5*stiffness*onset*failure;
    require(std::abs(work-expectedWork)<2.e-7,"Fixed-direction two-spring work must match the current triangular softening law, not an assumed A Gc budget");
}

#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
void testEnergyLawNativeGPUParity() {
    for(int mode=0;mode<7;++mode) {
        Fixture cpu,gpu;
        for(Fixture* f : {&cpu,&gpu}) {
            f->value().setStiffness(1000.,400.,200.,100.);
            f->value().setModeICriticalEnergy(.004);f->value().setModeIICriticalEnergy(.008);
        }
        gpu.particles.setUseDevice(true);gpu.bonds.setUseDevice(true);
        for(Real amplitude : {0.,.0005,.0015,.003,.006,.0002,.012,.04,.16}) {
            setEnergyWorkPose(cpu,mode,amplitude);setEnergyWorkPose(gpu,mode,amplitude);
            cpu.evaluate();gpu.evaluateGPU();
            const bond& a=cpu.value();const bond& b=gpu.value();
            require(close(a.damageFactor(),b.damageFactor(),1.e-8)&&a.getDamageStarted()==b.getDamageStarted(),"Native CPU/GPU energy damage and onset must agree on loading, unloading and rupture");
            require(close(a.getEffectiveFractureEnergy(),b.getEffectiveFractureEnergy(),1.e-8)
                &&close(a.getMaximumEffectiveDisplacement(),b.getMaximumEffectiveDisplacement(),1.e-8),"Native CPU/GPU must preserve the same frozen budget and irreversible energy history");
            require(math::norm(a.force()-b.force())<1.e-7&&math::norm(a.masterTorque()-b.masterTorque())<1.e-7
                &&math::norm(a.slaveTorque()-b.slaveTorque())<1.e-7&&close(total(a),total(b),1.e-8),"Native CPU/GPU force, couple and stored-energy degradation must agree");
        }
    }
    Fixture cpu,gpu;
    for(Fixture* f : {&cpu,&gpu}) { f->pose(0,.0014); }
    gpu.particles.setUseDevice(true);gpu.bonds.setUseDevice(true);
    cpu.evaluate();gpu.evaluateGPU();
    for(Fixture* f : {&cpu,&gpu}) f->pose(0,.004);
    cpu.evaluate();gpu.evaluateGPU();
    const Real damage=cpu.value().damageFactor();
    for(Fixture* f : {&cpu,&gpu}) f->pose(0,-.0002);
    cpu.evaluate();gpu.evaluateGPU();
    require(damage>0.&&damage<1.&&close(cpu.value().force().x,-.2)&&close(gpu.value().force().x,-.2),"CPU/GPU partially damaged bonds must retain the full compressive branch");
    require(close(cpu.value().damageFactor(),damage)&&close(gpu.value().damageFactor(),damage),"CPU/GPU compression cannot heal or advance tensile damage");
}
#endif

}
int main(int argc,char** argv) {
    try {
        if(argc==2&&std::string(argv[1])=="--work-probe") {
            std::cout<<"Diagnostic quadrature probe only; not a verification verdict.\n"
                <<"mode,steps,work_J,damage_dissipation_J,balance_error_J,A_Gc_J,budget_deviation_J\n"
                <<std::setprecision(17);
            for(int mode=0;mode<7;++mode) for(int steps : {4000,8000,16000}) {
                const auto measured=measureFractureWork(mode,steps);
                std::cout<<mode<<','<<steps<<','<<measured.work<<','<<measured.dissipated<<','
                    <<std::abs(measured.work-measured.dissipated)<<','<<measured.critical<<','
                    <<measured.work-measured.critical<<'\n';
            }
            return 0;
        }
        require(argc==1,"Unknown test option; use --work-probe only for diagnostic quadrature measurements.");
        testConfiguration(); testPureModes(); testStrengthAndWork(); testStrengthEnvelopeAndBudget();
        testBKScaleInvariance(); testLowEnergyStrengthCrossing(); testPureModeElasticLoads(); testSquaredElasticExpressions(); testSquaredModalResolution(); testBrittleFracture();
        testUnloadingAndCompressionRetention(); testRigidRotation(); testExactNativeTorque(); testEnergyEquivalentHistory(); testStiffnessPreconditions(); testFourEnergySeparation(); testRawStiffnessStorage(); testCompressionCannotCancelBending(); testFrozenBKMixity();
        testEuclideanMixedCoordinate();testZeroEnergyStrengthCrossing();testNativeFractureWorkBudget();testFixedDirectionTwoSpringWork();
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
        testSquaredNativeGPUParity();testBrittleNativeGPUParity();testLowEnergyNativeGPUParity();testEnergyLawNativeGPUParity();
#endif
        std::cout<<"Strength-initiated BK bond: all focused constitutive checks passed.\n";
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
    return 0;
}
