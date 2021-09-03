data <- TCGAanalyze_DMR(data, groupCol = "methylation_subtype",
                        group1 = "CIMP.H",
                        group2="CIMP.L",
                        p.cut = 10^-5,
                        diffmean.cut = 0.25,
                        legend = "State",
                        plot.filename = "coad_CIMPHvsCIMPL_metvolcano.png")