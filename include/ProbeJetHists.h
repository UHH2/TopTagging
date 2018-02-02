#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/TopJet.h"
#include "TH2D.h"

namespace uhh2examples {

class ProbeJetHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    ProbeJetHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;
    virtual void fill_probe(const uhh2::Event & ev, const TopJet & jet);
    virtual ~ProbeJetHists();

  private:

    TH2D *h_ratio_mLB_vs_mB;
    TH2D *h_tau32_vs_pt;
    bool fill_PDF;
    std::string mass_scale;
};

}
