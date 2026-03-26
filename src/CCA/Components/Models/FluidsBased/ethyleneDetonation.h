#ifndef Uintah_Component_Models_FluidsBased_ethlyeneDetonation_h
#define Uintah_Component_Models_FluidsBased_ethlyeneDetonation_h

#include <CCA/Components/Models/FluidsBased/FluidsBasedModel.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/GeometryPiece/GeometryPiece.h>

#include <vector>
#include <string>

#include <cmath>

namespace Uintah {


class ICELabel;

namespace SpeciesIndex {
  constexpr int C2H4 = 0;
  constexpr int O2   = 1;
  constexpr int N2   = 2;
  constexpr int CO2  = 3;
  constexpr int CO   = 4;
  constexpr int H2O  = 5;
  constexpr int H    = 6;
  constexpr int OH   = 7;
}

/**
 * ethyleneDetonation
 *
 * Single-material ICE combustion model with multiple passive scalars
 * representing pseudospecies mass fractions.
 */
class ethyleneDetonation : public FluidsBasedModel {

public:
  ethyleneDetonation(const ProcessorGroup* myworld,
                     const MaterialManagerP& materialManager,
                     const ProblemSpecP& params);

  virtual ~ethyleneDetonation();

  virtual void problemSetup(GridP& grid, const bool isRestart);

  virtual void scheduleInitialize(SchedulerP& sched,
                                  const LevelP& level);

  virtual void initialize(const ProcessorGroup*,
                          const PatchSubset* patches,
                          const MaterialSubset* matls,
                          DataWarehouse* old_dw,
                          DataWarehouse* new_dw);

  virtual void scheduleComputeModelSources(SchedulerP& sched,
                                           const LevelP& level);

  virtual void computeModelSources(const ProcessorGroup*,
                                   const PatchSubset* patches,
                                   const MaterialSubset* matls,
                                   DataWarehouse* old_dw,
                                   DataWarehouse* new_dw);


  virtual void outputProblemSpec(ProblemSpecP& ps);
  virtual void scheduleRestartInitialize(SchedulerP&, const LevelP&);
  virtual void scheduleTestConservation(SchedulerP&, const PatchSet*);

  virtual void scheduleComputeStableTimeStep(SchedulerP&, const LevelP&) {}
  virtual void scheduleModifyThermoTransportProperties(SchedulerP&, const LevelP&, const MaterialSet*) {}
  virtual void computeSpecificHeat(CCVariable<double>&, const Patch*, DataWarehouse*, const int) {}
  virtual void scheduleErrorEstimate(const LevelP&, SchedulerP&) {}

private:
  ethyleneDetonation(const ethyleneDetonation&) = delete;
  ethyleneDetonation& operator=(const ethyleneDetonation&) = delete;

  //------------------------------------------------------------------
  // Geometry-based initialization
  //------------------------------------------------------------------
  struct Region {
    GeometryPieceP piece;
    std::vector<double> Yinit; // size = 6 (tracked species)

    Region(GeometryPieceP p, const std::vector<double>& Y)
      : piece(p), Yinit(Y) {}
  };
  //------------------------------------------------------------------
  // Combustion Function declarations
  //------------------------------------------------------------------
  double arrhenius(double T, int idx, double M);
  double globalRates(double T, const std::vector<double>& Y, double rho, int rate);
  std::vector<double> molarEnthalpy(double T);
  double heatRelease(double T, double rho, const std::vector<double>& Y, const std::vector<double>& nu1, const std::vector<double>& nu2);
  std::vector<double> massSource(double T, double rho, const std::vector<double>& Y, const std::vector<double>& nu1, const std::vector<double>& nu2);
  //------------------------------------------------------------------
  // Constants
  //------------------------------------------------------------------
  // Mass tolerance
  static constexpr double massTol = 1e-12;
  
  // Molar Masses                              [C2H4,  O2,     CO2,    CO,    H2O,      H,     OH] g/mol
  inline static const std::vector<double> MM = {28.05, 31.999, 44.009, 28.01, 18.01528, 1.008, 17.008};

  // Molar Masses with N2                         [C2H4,  N2,    O2,     CO2,    CO,    H2O,      H,     OH] g/mol
  inline static const std::vector<double> MM_N2 = {28.05, 28.02, 31.999, 44.009, 28.01, 18.01528, 1.008, 17.008};

  // 1.2 Chaperon Efficiencies
  inline static const std::vector<double> chap2 = {1.0, 0.3, 1.0, 1.5, 0.75, 7.0, 1.0, 1.0};
  inline static const std::vector<double> chap3 = {1.0, 1.0, 1.0, 3.8, 1.9, 12.0, 1.0, 1.0};

  // Universial gas constant
  static constexpr double R = 0.008314462618; //(kJ / mol K)

  // Arrhenius Parameters
  inline static const std::vector<double> A = {
        3.52e16, 2.60e19, 2.20e22, 2.18e23,
        4.40e6,  4.97e8,  4.6e15,  6.26e16,
        0.0,     4.22e13, 0.0,     2.00e12
  };

  inline static const std::vector<double> n = {
        -0.70, -1.20, -2.00, -1.93,
         1.50,  1.50, -0.54,  0.00,
         0.0,   0.00,  0.0,   0.00
  };

  inline static const std::vector<double> E = {
        71.4, 0.0, 0.0, 499.0,
       -3.1, 89.7, 188.0, 326.0,
        0.0, 241.0, 0.0, 20.9
  };

  // NASA polynomial coefficients for molar enthalpy [C2H4, O2, CO2, CO, H2O, H, OH]
  inline static const std::vector<double> a1 = {
    2.03611116, 3.28253784, 3.85746029, 
    2.71518561, 3.03399249, 2.50000001, 3.09288767};

  inline static const std::vector<double> a2 = {
    7.32270755e-03, 7.41543770e-04, 2.20718513e-03, 
    1.03126372e-03, 1.08845902e-03, -1.15421486e-11, 2.74214858e-04};
  
  inline static const std::vector<double> a3 = {
    -2.23692638e-06, -2.52655556e-07, -7.38271347e-07, 
    -3.32941924e-07, -5.46908393e-08, 5.38539827e-15, 4.21684093e-08};
  
  inline static const std::vector<double> a4 = {
    3.68057308e-10, 5.23676387e-11, 1.30872547e-10,
     5.75132520e-11, -2.42604967e-11, -1.18378809e-18, -2.19865389e-11};
  
  inline static const std::vector<double> a5 = {
    -2.51412122e-14, -4.33435588e-15, -9.44168328e-15,
     -4.07295432e-15, 3.36401984e-15, 9.96394714e-23, 2.34824752e-15};
  
  inline static const std::vector<double> a6 = {
    4939.88614, -1088.45772, -48759.166,
    -14151.8724, -30004.2971, 25473.6599, 3858.657};
  
  // Stoichmetric Coefficients
  inline static const std::vector<double> nu1 = {-1.0, -2.0, 0.0, 2.0, 4.0/3.0, 2.0/3.0, 2.0/3.0};
  inline static const std::vector<double> nu2 = {0.0, -0.30268766157832044, 0.6084054799653067, -0.6084054799653067, 0.6381929192026506, -0.6351627623939848, -0.6412230760113165};

  //------------------------------------------------------------------
  // Data members
  //------------------------------------------------------------------
  ICELabel*    Ilb{nullptr};
  Material*    d_matl{nullptr};
  MaterialSet* d_matl_set{nullptr};
  ProblemSpecP d_params;

  // Species bookkeeping
  static constexpr int N_SPECIES = 6;

  // VarLabels for passive scalars and their sources
  std::vector<VarLabel*> d_Y_labels;      // scalar-YC2H4, scalar-YO2, ...
  std::vector<VarLabel*> d_Y_src_labels;  // scalar_YC2H4_src, ...

  // Geometry regions for initialization
  std::vector<Region*> d_regions;

  // Chemistry parameters from UPS
//   std::vector<double> d_h;    // size 7
//   std::vector<double> d_nu1;  // size 7
//   std::vector<double> d_nu2;  // size 7

  bool d_debug{false};
};

} // namespace Uintah

#endif
