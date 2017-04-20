#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/TopJet.h"

namespace uhh2examples {

class ProbeJetHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    ProbeJetHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;
    virtual void fill_probe(const uhh2::Event & ev, const TopJet & jet);
    virtual ~ProbeJetHists();
};

}
