//---------------------------------------------------------------
// Model is for combustion of Ethylene air mixture in presence of shocks (detonations)
// The implemented model is from "Two-step chemical-kinetic descriptions
// for hydrocarbon–oxygen-diluent ignition and detonation applications" (2005)
// B. Varatharajan et al

// Enthalpy values are from Nasa Polynomials (gri-mech)
// http://combustion.berkeley.edu/gri-mech/data/nasa_plnm.html

// Written by James Karr Feb 2026

//--------------------------------------------------------------

#include <CCA/Components/Models/FluidsBased/ethyleneDetonation.h>

#include <CCA/Ports/Scheduler.h>
#include <CCA/Components/ICE/Core/ICELabel.h>

#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Material.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/GeometryPiece/GeometryPieceFactory.h>

#include <numeric>


using namespace Uintah;
using namespace SpeciesIndex;

//------------------------------------------------------------------
// Constructor / Destructor
//------------------------------------------------------------------
ethyleneDetonation::ethyleneDetonation(const ProcessorGroup* myworld,
                                       const MaterialManagerP& materialManager,
                                       const ProblemSpecP& params)
  : FluidsBasedModel(myworld, materialManager),
    d_params(params)
{
  Ilb = scinew ICELabel();
}

ethyleneDetonation::~ethyleneDetonation()
{
  if (d_matl_set && d_matl_set->removeReference()) {
    delete d_matl_set;
  }

  for (auto* lbl : d_Y_labels) {
    VarLabel::destroy(lbl);
  }
  for (auto* lbl : d_Y_src_labels) {
    VarLabel::destroy(lbl);
  }

  for (auto* r : d_regions) {
    delete r;
  }

  delete Ilb;
}

//------------------------------------------------------------------
// Output UPS
//------------------------------------------------------------------
void ethyleneDetonation::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP model_ps = ps->appendChild("Model");
  model_ps->setAttribute("type", "ethyleneDetonation");

  d_matl->outputProblemSpec(model_ps);

  ProblemSpecP ed_ps = model_ps->appendChild("ethyleneDetonation");
  ed_ps->appendElement("debug", d_debug);

}

void ethyleneDetonation::scheduleRestartInitialize(SchedulerP&, const LevelP&) {}
void ethyleneDetonation::scheduleTestConservation(SchedulerP&, const PatchSet*) {}

//------------------------------------------------------------------
// problemSetup
//------------------------------------------------------------------
void ethyleneDetonation::problemSetup(GridP&, const bool)
{
  ProblemSpecP ps = d_params->findBlock("ethyleneDetonation");
  if (!ps) {
    throw ProblemSetupException("Missing <ethyleneDetonation> block", __FILE__, __LINE__);
  }

  d_matl = m_materialManager->parseAndLookupMaterial(ps, "material");

  std::vector<int> m(1);
  m[0] = d_matl->getDWIndex();
  d_matl_set = scinew MaterialSet();
  d_matl_set->addAll(m);
  d_matl_set->addReference();

  ps->getWithDefault("debug", d_debug, false);

  //----------------------------------------------------------------
  // Create 6 passive scalars
  //----------------------------------------------------------------
  static const char* names[N_SPECIES] = {
    "YC2H4", "YO2", "YCO2", "YCO", "YH2O", "YH"
  };

  for (int i = 0; i < N_SPECIES; i++) {
    std::string yname = std::string("scalar-") + names[i];
    std::string sname = std::string("scalar_") + names[i] + "_src";

    VarLabel* Y = VarLabel::create(yname, CCVariable<double>::getTypeDescription());
    VarLabel* S = VarLabel::create(sname, CCVariable<double>::getTypeDescription());

    d_Y_labels.push_back(Y);
    d_Y_src_labels.push_back(S);

    registerTransportedVariable(d_matl_set, Y, S);
  }

  //----------------------------------------------------------------
  // Geometry-based initialization
  //----------------------------------------------------------------
  for (ProblemSpecP geom_ps = ps->findBlock("geom_object");
       geom_ps != nullptr;
       geom_ps = geom_ps->findNextBlock("geom_object")) {

    std::vector<GeometryPieceP> pieces;
    GeometryPieceFactory::create(geom_ps, pieces);

    std::vector<double> Yinit;
    geom_ps->require("Y", Yinit);

    if (Yinit.size() != N_SPECIES) {
    throw ProblemSetupException(
        "Initial Y vector must have length 6",
        __FILE__, __LINE__);
    }

    for (auto& piece : pieces) {
    d_regions.push_back(scinew Region(piece, Yinit));
    }

  }
}

//------------------------------------------------------------------
// Initialization
//------------------------------------------------------------------
void ethyleneDetonation::scheduleInitialize(SchedulerP& sched,
                                            const LevelP& level)
{
  Task* t = scinew Task("ethyleneDetonation::initialize",
                        this, &ethyleneDetonation::initialize);

  for (auto* lbl : d_Y_labels) {
    t->computesVar(lbl);
  }

  sched->addTask(t, level->eachPatch(), d_matl_set);
}

void ethyleneDetonation::initialize(const ProcessorGroup*,
                                    const PatchSubset* patches,
                                    const MaterialSubset* matls,
                                    DataWarehouse*,
                                    DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    for (int m = 0; m < matls->size(); m++) {
      int indx = matls->get(m);

      std::vector<CCVariable<double>> Y(N_SPECIES);
      //                  Y0 = [YC2H4,              YO2,                 YCO2, YCO, YH2O, YH]
      std::vector<double> Y0 = {0.0637524374728892, 0.21814541536937784, 0.0, 0.0, 0.0, 0.0};
      for (int k = 0; k < N_SPECIES; k++) {
        new_dw->allocateAndPut(Y[k], d_Y_labels[k], indx, patch);
        Y[k].initialize(Y0[k]);
      }

      for (CellIterator iter(patch->getCellIterator()); !iter.done(); iter++) {
        Point pt = patch->cellPosition(*iter);

        for (auto* r : d_regions) {
          if (r->piece->inside(pt)) {
            for (int k = 0; k < N_SPECIES; k++) {
              Y[k][*iter] = r->Yinit[k];
            }
          }
        }
      }
    }
  }
}

//------------------------------------------------------------------
// Source terms
//------------------------------------------------------------------
void ethyleneDetonation::scheduleComputeModelSources(SchedulerP& sched,
                                                     const LevelP& level)
{
  Task* t = scinew Task("ethyleneDetonation::computeModelSources",
                        this, &ethyleneDetonation::computeModelSources);

  t->requiresVar(Task::OldDW, Ilb->delTLabel);
  t->requiresVar(Task::OldDW, Ilb->temp_CCLabel, Ghost::None, 0);
  t->requiresVar(Task::OldDW, Ilb->rho_CCLabel, Ghost::None, 0);

  t->modifiesVar(Ilb->modelEng_srcLabel);

  for (int k = 0; k < N_SPECIES; k++) {
    t->requiresVar(Task::OldDW, d_Y_labels[k], Ghost::None, 0);
    t->modifiesVar(d_Y_src_labels[k]);
  }

  sched->addTask(t, level->eachPatch(), d_matl_set);
}

// ----------------------------------------------------------------
// Combustion Functions
// ----------------------------------------------------------------
// Step 1 calculate global rates (omegaI, omegaII)
double ethyleneDetonation::globalRates(double T,
                                       const std::vector<double>& Y,
                                       double rho,
                                       int rate)
{
    // Y is species mass fraction of all species, including nitrogen

    //Convert Density
    double rho_g = rho / 1000; // g /cm^3
    
    // 1.1) Molar Concentration mol / cm^3
    std::vector<double> conc(Y.size());
    for(size_t j = 0; j< Y.size(); j++){
        conc[j] = rho_g * Y[j] / MM_N2[j];
    }

    if(rate == 1){
        // if statements for speed. Rate won't be calculated in cells with no oxygen or fuel
        if(Y[O2] < massTol || Y[C2H4] < massTol){
          return 0.0;
        }

        if(Y[H] < massTol){ //avoid calculating 5 expensive rates if they will be multiplied by 0
          double k10 = arrhenius(T, 10, 0.0);
          return k10 * conc[C2H4] * conc[O2];
        }
        // 1.3 Third Body Reaction terms
        double iM13 = std::accumulate(conc.begin(), conc.end(), 0.0);
        double iM2  = std::inner_product(conc.begin(), conc.end(), chap2.begin(), 0.0);

        // 1.4 Calculate Elmentary rates
        double k1  = arrhenius(T, 1, 0.0);
        double k2  = arrhenius(T, 2, 0.0);
        double k10 = arrhenius(T, 10, 0.0);
        double k11 = arrhenius(T, 11, iM13);
        double k12 = arrhenius(T, 12, 0.0);
        double k13 = arrhenius(T, 13, iM13);
        
        

        // 1.5 Rate Correction Factors
        double t = k12 * conc[O2] / (k12 * conc[O2] + k13 * iM13);
        double q = k1 / (k1 + k2 * iM2 + t * k11);

        //Heaviside Function return 1 if fuel is still present, 0 if not
        double heavy;
        if(conc[C2H4] > 1e-12){
            heavy = 1.0;
        }else{
            heavy = 0.0;
        }

        // 1.6 Global reaction rate 1 (omegaI)
        return k10 * conc[C2H4] * conc[O2] + q * k1 * conc[H] * conc[O2] * heavy; // mol / cm^3 s
    } else if(rate == 2){
        double f = 2.0 / 3.0; // constant
        // if statement for speed. Rate won't be calculated in cells with no H or OH and H2O
        if((Y[H] < massTol || Y[OH] < massTol) && Y[H2O] < massTol){
          return 0.0;
        }
   
        if(Y[H] < massTol || Y[OH] < massTol){//avoid calculating expensive rates if they will be multiplied by 0
          // 1.3 Third Body Reaction term
          double iM3  = std::inner_product(conc.begin(), conc.end(), chap3.begin(), 0.0);
          // 1.4 Elementary reaction rate
          double  k4 = arrhenius(T, 4, 0.0);
          return -k4 * iM3 * conc[H2O] / f;
        }
        if(Y[H2O] < massTol){//avoid calculating expensive rates if they will be multiplied by 0
          // 1.3 Third Body Reaction term
          double iM3  = std::inner_product(conc.begin(), conc.end(), chap3.begin(), 0.0);
          // 1.4 Elementary reaction rate
          double k3 = arrhenius(T, 3, 0.0);
          return k3 * iM3 * conc[H] * conc[OH] / f;
        }

        // 1.3 Third Body Reaction term
        double iM3  = std::inner_product(conc.begin(), conc.end(), chap3.begin(), 0.0);
        // 1.4 Elementary reaction rate
        double k3 = arrhenius(T, 3, 0.0);
        double k4 = arrhenius(T, 4, 0.0);

        // 1.6 Global reaction rate 2 (omega II)
        return (k3 * iM3 * conc[H] * conc[OH] - k4 * iM3 * conc[H2O]) / f; // mol / cm^3 s
    }
    return 0.0;
}

// Step 1.4
double ethyleneDetonation::arrhenius(double T, 
                                     int idx, 
                                     double M)
{
    if (idx == 11 || idx == 13) {

        double k0   = 0.0;
        double kinf = 0.0;

        if (idx == 11) {
            k0   = 1.9e35 * std::pow(T, -5.57)
                 * std::exp(-21.1 / (R * T));

            kinf = 1.08e12 * std::pow(T, 0.45)
                 * std::exp(-7.6 / (R * T));
        }
        else { // idx == 13
            k0   = 3.99e33 * std::pow(T, -4.99)
                 * std::exp(-167.4 / (R * T));

            kinf = 1.11e10 * std::pow(T, 1.04)
                 * std::exp(-153.8 / (R * T));
        }

        double Pr = std::max(k0 * M / kinf, 1e-300);
        double Fc = 0.832 * std::exp(-T / 1203.0);
        
        double logFc = std::log10(Fc);
        double N     = 0.75 - 1.27 * logFc;
        double base  = std::log10(Pr) / N; 
        double logF  = logFc / (1.0 + base * base);

        double F  = std::pow(10.0, logF);
        return kinf * Pr * F /(1.0 + Pr);
    }
    int i = idx - 1;

    return A[i] * std::pow(T, n[i]) * std::exp(-E[i] / (R * T));
}

std::vector<double> ethyleneDetonation::molarEnthalpy(double T)
{
  // Solve for molar enthalpies of C2H4, O2, CO2, CO, H2O, H, OH using nasa polynomials (mech-gri). See top for reference
  std::vector<double> h_species;
  double h;
  for(int j = 0; j < 7; j++){
    h = 1000 * R * (a6[j] + a1[j] * T + a2[j] * std::pow(T, 2) + a3[j] * std::pow(T, 3) + a4[j] * std::pow(T, 4) + a5[j] * std::pow(T, 5)); // J / mol
    h_species.push_back(h);
  }
  return h_species;
}

// Step 2 Calculate Heat release (q)
double ethyleneDetonation::heatRelease(double T,
                                       double rho,
                                       const std::vector<double>& Y,
                                       const std::vector<double>& nu1,
                                       const std::vector<double>& nu2)
{
    // 2.1 Call globalRates and convert from mol / cm^3 s to mol / m^3 s
    double omegaI  = 1e6 * globalRates(T, Y, rho, 1); // mol / m^3 s
    double omegaII = 1e6 * globalRates(T, Y, rho, 2); // mol / m^3 s

    // 2.2 Species Enthalpies   [C2H4, O2, CO2, CO, H2O, H, OH]
    std::vector<double> h = molarEnthalpy(T); //J /mol

    // 2.3 Calculate heat release per reaction
    double qI  = omegaI  * std::inner_product(nu1.begin(), nu1.end(), h.begin(), 0.0); // W/m^3
    double qII = omegaII * std::inner_product(nu2.begin(), nu2.end(), h.begin(), 0.0); // W/m^3

    // 2.4 Total heat release of reaction
    return -(qI + qII);

}

std::vector<double> ethyleneDetonation::massSource(double T,
                                                   double rho,
                                                   const std::vector<double>& Y,
                                                   const std::vector<double>& nu1,
                                                   const std::vector<double>& nu2)
{
    // 3.1 Call globalRates
    double omegaI  = globalRates(T, Y, rho, 1); //mol /cm^3 s
    double omegaII = globalRates(T, Y, rho, 2); //mol /cm^3 s

    std::vector<double> S;
    double tmp;
    for(size_t k = 0; k < MM.size(); k++){
      // 3.2 Calculate species creation/destruction for both global reactions
        tmp = omegaI * nu1[k] + omegaII * nu2[k]; // mol /cm^3 s
      // 3.3 Convert to appropiate units 
        tmp = MM[k] * tmp * 1e3; // kg / m^3 s
        S.push_back(tmp);
    }
    return S;
}


void ethyleneDetonation::computeModelSources(const ProcessorGroup*,
                                             const PatchSubset* patches,
                                             const MaterialSubset* matls,
                                             DataWarehouse* old_dw,
                                             DataWarehouse* new_dw)
{
  delt_vartype delT;
  old_dw->get(delT, Ilb->delTLabel);

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    Vector dx = patch->dCell();
    double cellVol = dx.x() * dx.y() * dx.z();

    for (int m = 0; m < matls->size(); m++) {
      int indx = matls->get(m);

      CCVariable<double> eng_src;
      new_dw->getModifiable(eng_src, Ilb->modelEng_srcLabel, indx, patch);

      std::vector<constCCVariable<double>> Yold(N_SPECIES);
      std::vector<CCVariable<double>>       Ysrc(N_SPECIES);

      for (int k = 0; k < N_SPECIES; k++) {
        old_dw->get(Yold[k], d_Y_labels[k], indx, patch, Ghost::None, 0);
        new_dw->getModifiable(Ysrc[k], d_Y_src_labels[k], indx, patch);
      }
      // Pull in Temperature and density from data warehouse
      constCCVariable<double> temp;
      constCCVariable<double> rho;

      old_dw->get(temp, Ilb->temp_CCLabel, indx, patch, Ghost::None, 0);
      old_dw->get(rho,  Ilb->rho_CCLabel,  indx, patch, Ghost::None, 0);

      for (CellIterator iter(patch->getCellIterator()); !iter.done(); iter++) {

        // Current Properties for cell
        double T      = temp[*iter];
        double rho_kg = rho[*iter];
        double YN2    = 0.718102147157733; //stoichmetric air ethylene (will be specified in input file)

        // Build the mass fraction vector for all species [C2H4, O2, N2, CO2, CO, H2O, H, OH]
        std::vector<double> Y;
        double Ytmp;

        for (int j = 0; j< N_SPECIES; j++){
            Ytmp = Yold[j][*iter];
            Y.push_back(Ytmp);
        }
        Y.insert(Y.begin() + 2, YN2);
        double YOH = std::max(1 - std::accumulate(Y.begin(), Y.end(), 0.0), 0.0);
        Y.push_back(YOH);

        // Compute Scalar sources
        std::vector<double> S = massSource(T, rho_kg, Y, nu1, nu2);
        for (int j = 0; j< N_SPECIES; j++){
            Ysrc[j][*iter] += S[j] * delT / rho_kg;
        }

        // Compute Energy source
        eng_src[*iter] += heatRelease(T, rho_kg, Y, nu1, nu2) * cellVol * delT;
      }
    }
  }
}
