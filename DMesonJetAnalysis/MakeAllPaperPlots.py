#!/usr/bin/env python
# python script to do extract B feed down correction factors

import IPython
import BFeedDown_JetZSpectrum_JetPt_5_15_Paper
import BFeedDown_JetZSpectrum_JetPt_15_30_Paper
import BFeedDown_JetPtSpectrum_DPt_30_Paper
import Efficiency_Paper
import JetPtSpectrum_Paper
import JetZSpectrum_JetPt_5_15_Paper
import JetZSpectrum_JetPt_15_30_Paper
import Uncertainties_JetPtSpectrum_Paper
import SideBandInvMass_Paper


def main():
    BFeedDown_JetZSpectrum_JetPt_5_15_Paper.main()
    BFeedDown_JetZSpectrum_JetPt_15_30_Paper.main()
    BFeedDown_JetPtSpectrum_DPt_30_Paper.main()
    Efficiency_Paper.main()
    JetPtSpectrum_Paper.main()
    JetZSpectrum_JetPt_5_15_Paper.main()
    JetZSpectrum_JetPt_15_30_Paper.main()
    Uncertainties_JetPtSpectrum_Paper.main()
    SideBandInvMass_Paper.main()


if __name__ == '__main__':
    main()

    IPython.embed()
