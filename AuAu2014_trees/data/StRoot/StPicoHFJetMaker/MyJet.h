#pragma once

#ifndef __CINT__
#include "fastjet/PseudoJet.hh"
#endif

#include <vector>

class MyJet {

public:
  float pt;
  float pt_corr;
  float eta;
  float phi;
  float area;
  float rho;
  float pt_lead;
  int n_constituents;
  float neutral_fraction;
  bool trigger_match;

  MyJet()
      : pt(-999), pt_corr(-999), eta(-999), phi(-999), area(-999), rho(-999), pt_lead(-999),
        n_constituents(-999), neutral_fraction(-999), trigger_match(false) {}

  MyJet(float pt, float pt_corr, float eta, float phi, float area, float rho,
        float pt_lead, int n_constituents, float neutral_fraction,
        bool trigger_match)
      : pt(pt), pt_corr(pt_corr), eta(eta), phi(phi), area(area), rho(rho),
        pt_lead(pt_lead), n_constituents(n_constituents),
        neutral_fraction(neutral_fraction), trigger_match(trigger_match) {}

  float deltaR(const MyJet &other) const {
    if (pt < 0 || other.pt < 0) {
      return 10000; // Return 10000 to indicate invalid jets
    }
    float deta = eta - other.eta;
    float dphi = TVector2::Phi_mpi_pi(phi - other.phi);
    return sqrt(deta * deta + dphi * dphi);
  }

#ifndef __CINT__
MyJet(fastjet::PseudoJet jet, float rho)
  : rho(rho) {
  pt   = jet.perp();
  eta  = jet.eta();
  phi  = jet.phi();
  area = jet.area_4vector().pt(); 
  pt_corr = pt - area * rho;

  vector<fastjet::PseudoJet> constituents = sorted_by_pt(jet.constituents());
  n_constituents = constituents.size();

  if (n_constituents == 0 || pt <= 0) {
    // mark as invalid jet
    pt          = -999;
    pt_corr     = -999;
    pt_lead     = -999;
    neutral_fraction = -999;
    trigger_match    = false;
    return;
  }

  pt_lead = constituents[0].perp();

  float neutral_sum = 0.0;
  trigger_match = false;
  for (const auto &constituent : constituents) {
    int uidx = constituent.user_index();
    if (uidx == 9999) trigger_match = true;
    if (uidx == 0 || uidx == 9999)
      neutral_sum += constituent.perp();
  }
  neutral_fraction = neutral_sum / pt;
}

#endif
};

typedef pair<MyJet, MyJet> MatchedJetPair;
