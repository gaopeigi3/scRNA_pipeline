hierarchical_markers = {
    "CD4": {
        # "general": ["CD4", "IL7R", "MAL", "LEF1", "LTB", "LDHB", "TPT1", "TRAC", "TCF7", "TMSB10", "LEPROTL1","CD3D", "CD3E","CD3G", "CD27"],   # broad CD4 T cell markers
         "general": ["CD4", "IL7R", "MAL", "LTB"],   # broad CD4 T cell markers
        "subtypes": {
            "CD4 Memory T cells": ["CXCR3", "IL7R", "LTB"],
            "CD4 Naive T cells": ["LEF1","TCF7","SELL","CD44"], #Naive "IL7R", "LTB"
            "CD4 CTL T cells": ["GZMB","PRF1","GNLY","NKG7"],
            "CD4 Exhausted T cells": ["CTLA4","LAG3","TIGIT","HAVCR2","PDCD1"],
            "CD4 Th1 T cells": ["TBX21","STAT4","IFNG","IL12A"],
            "CD4 Th2 T cells": ["STAT6","GATA6","IL4"],
            "CD4 Th17 T cells": ["IL7RA","STAT3","RORC"],
            "CD4 Tfh T cells": ["BCL6","CXCR5"],
            "Treg": ["TGFB1","FOXP3","IL2RA","IKZF2"]
        }
    },
    "CD8": {
        # "general": ["CD8A","CD8B","NELL2", "CD3D", "CD3E", "S100B", "GZMH", "CD3G", "TRGC2", "CCL5"],
        "general": ["CD8A","CD8B", "CD3D" ],
        "subtypes": {
            "CD8 Naive T cells": ["NELL2","CD8B","CCR7","LEF1","SELL","IL7R"],
            "CD8 Effector Memory T cells": ["PRF1","FGFBP2","FCGR3A","KLRD1","NKG7","KLRG1","CX3CR1"],
            "CD8 Exhausted T cells": ["CTLA4","LAG3","TIGIT","PDCD1","TOX"],
            "CD8 CTL T cells": ["GZMB","PRF1","GNLY","NKG7","GZMA"],
            "MAIT": ["KLRB1","IL7R","SLC4A10","RORC"]
        }
    },
    # "B": {
    #     "general": ["FCRL1", "FCRL2", "CD22", "ARHGAP24", #Azimuth
    #                                                            "BANK1", "MS4A1", "RALGPS2", "CD37",
    #                                                            "SWAP70", "CD79A",
    #                                                            "CD19", "B220"],#Abbas Lab
    #     "subtypes": {}
    # },
    "B": {
        "general": ["MS4A1", "CD79A", "CD79B", "BANK1", "RALGPS2"],
        "subtypes": {
            "B intermediate": ["MS4A1", "TNFRSF13B", "IGHM", "IGHD", "AIM2", "CD79A", "LINC01857", "RALGPS2", "BANK1", "CD79B"],
            "B memory": ["MS4A1", "COCH", "AIM2", "BANK1", "SSPN", "CD79A", "TEX9", "RALGPS2", "TNFRSF13C", "LINC01781"],
            "B naive": ["IGHM", "IGHD", "CD79A", "IL4R", "MS4A1", "CXCR4", "BTG1", "TCL1A", "CD79B", "YBX3"],
            "Plasmablast": ["IGHA2", "MZB1", "TNFRSF17", "DERL3", "TXNDC5", "TNFRSF13B", "POU2AF1", "CPNE5", "HRASLS2", "NT5DC2"]
        }
    },

    "pre B": {
        "general": ["NPY", "LCN6", "RAG2", "HMHB1",
                                                               "ARPP21", "AKAP12", "RAG1", "C10orf10",
                                                               "CYGB", "SLC8A1-AS1"],#Abbas Lab
        "subtypes": {}
    },

    # ---- Myeloid / Mono / Macrophage ----
    "Mono": {
        "general": ["LYPD2", "FOLR3", "CLEC4E", "LILRA1",
                    "CDA", "RBP7", "CD300LF", "FPR1", "CD93", "MTMR11"],
        "subtypes": {
            "CD14 Mono": ["FOLR3", "CLEC4E", "MCEMP1", "RBP7",
                          "CDA", "FPR1", "CD300E", "C5AR1",
                          "CD93", "APOBEC3A"],
            "CD16 Mono": ["LYPD2", "VMO1", "TPPP3", "C1QA",
                          "C5AR1", "CD300E", "GPBAR1", "LILRA1",
                          "HES4", "APOBEC3A"]
        }
    },
    # ===== NK lineage =====
    "NK": {
        "general": ["TRDC", "FCER1G", 
                    "KLRF1" 
                    ],
        # "general": ["NKG7", "KLRC1", "KLRD1", "KLRB1", 
        #                                                        "KIR2DL1", "KIR3DL1", "KIR2DL3", "KIR3DL2",
        #                                                        "KIR3DL3", "GNLY", "NCAM1", "FCGR3A"],
        "subtypes": {
            "NK CD56-dim": ["GNLY", "TYROBP", "NKG7", "GZMB",  "PRF1", "FGFBP2", "SPON2"],
            "NK Proliferating": ["MKI67" , "TYMS",  "TOP2A", "PCLAF", "CD247", "CLSPN", "ASPM"],
            "NK CD56-bright": ["XCL2",  "SPINK2", "KLRC1", "XCL1", "SPTSSB", "PPP1R9A", "NCAM1", "TNFRSF11A"]
        }
    },
    "EMPs": {
        "general": ["MYCT1", "CRHBP", "NPR3", "AVP", "HPGDS", "CRYGD","IGSF10", "PBX1", 
                    "CYTL1","GATA2"
                    ],
        
        "subtypes": {
            "Megakaryocyte": [
                "GFI1B", "SELP", "GP1BA", "CD9", "ITGA2B", "GATA2", "FLI1",
                "GP1BB", "VWF", "THPO", "ELF1", "THBS1", "MPIG6B", "GP9",
                "F2R", "FOG1", "NFE2", "SPI1", "PF4"
            ],
            "Erythroid progenitor": [
                "GATA1", "KLF1", "FCER1A", "ITAG2B", "EPOR", "HBD", "ZFPM1",
                "GATA2", "GYPA", "TFRC", "TFR2", "CSF2RB", "APOE", "APOC1",
                "CNRIP1", "FOXO3", "ETS1", "BRD1", "TAL1"
            ]
        }
    },
    # ---- DC ----
    "DC": {
        "general": ["CLEC4C", "PROC", "SCT", "SCN9A",
                    "SHD", "PPM1J", "ENHO", "CLEC10A",
                    "LILRA4", "DNASE1L3"],
        "subtypes": {
            'ASDC':['PPP1R14A', 'LILRA4', 'AXL', 'IL3RA', 'SCT', 'SCN9A', 'LGMN', 'DNASE1L3', 'CLEC4C', 'GAS6'],
            'cDC1':['CLEC9A', 'DNASE1L3', 'C1orf54', 'IDO1', 'CLNK', 'CADM1', 'FLT3', 'ENPP1', 'XCR1', 'NDRG2'],
            'cDC2':['FCER1A', 'HLA-DQA1', 'CLEC10A', 'CD1C', 'ENHO', 'PLD4', 'GSN', 'SLC38A1', 'NDRG2', 'AFF3'],
            'pDC':['ITM2C', 'PLD4', 'SERPINF1', 'LILRA4', 'IL3RA', 'TPM2', 'MZB1', 'SPIB', 'IRF4', 'SMPD3']
        }
    },

    # ---- Progenitors / Stem cells ----
    "HSC": {
        "general": ["CD34", "AVP", "CRHBP"],
        "subtypes": {
            'MonocyticLineage' : ["LYZ", "S100A8", "S100A9", "ITGAM"],
            
            'Antigen' : ["CD74", "HLA-DRA", "HLA-DRB1", "CIITA"],
        }
    },
    "CLP": {
        "general": ["ACY3", "PRSS2", "C1QTNF4", "SPINK2",
                    "SMIM24", "NREP", "CD34", "DNTT", "FLT3", "SPNS3"],
        "subtypes": {}
    },
    "GMP": {
        "general": ["SERPINB10", "RNASE3", "MS4A3", "PRTN3",
                    "ELANE", "AZU1", "CTSG", "RNASE2", "RETN", "NPW"],
        "subtypes": {}
    },

    # ---- Erythroid / Megakaryocyte ----
    # "Early Erythroid": {
    #     "general": [
    #         "CNRIP1",  "ITGA2B", "TFR2", "MAP7", "FSCN1", "APOC1" "GATA1", "KLF1", 
    #                 "CYTL1","GATA2",
    #                 ],
    #     "subtypes": {}
    # },
    # "Late Erythroid": {
    #     "general": ["CTSE", "TSPO2", "IFIT1B", "TMEM56",
    #                            "RHCE", "RHAG", "SPTA1", "ADD2",
    #                            "EPCAM", "HBG1"],
    #     "subtypes": {}
    # },
    "Macrophage": {
        "general": ["SPIC", "FABP3", "CD5L", "CCL18",
                    "C1QC", "C1QB", "FABP4", "C1QA",
                    "APOE", "SELENOP"],
        "subtypes": {}
    },
    "Erythroid": {
        "general": ['HBD', 'HBA1', 'HBA2', 'AHSP', 'ALAS2', 'SLC4A1'],
        "subtypes": {
            "Early Erythroid": ["CNRIP1", "GATA2", "ITGA2B", "TFR2",
                                "GATA1", "KLF1", "CYTL1", "MAP7",
                                "FSCN1", "APOC1"],
            "Late Erythroid": ["CTSE", "TSPO2", "IFIT1B", "TMEM56",
                               "RHCE", "RHAG", "SPTA1", "ADD2",
                               "EPCAM", "HBG1"]
        }
    # # },"Erythroid": {

    # "subtypes": {
    #     # 1) BFU-E (earliest erythroid progenitor-like)
    #     # CD34 high, CD36 low/-, TFRC low, still hematopoietic background
    #     "BFU_E": [
    #         'CD34', 'PTPRC', 'TFRC', 'CD44'
    #     ],

    #     # 2) CFU-E (committed erythroid progenitor)
    #     # CD34 down, CD36 up, TFRC high
    #     "CFU_E": [
    #         'CD36', 'TFRC', 'CD44'
    #     ],

    #     # 3) Proerythroblast (start terminal differentiation)
    #     # GYPA turns on, TFRC high, EPOR may appear
    #     "Proerythroblast_ProE": [
    #         'GYPA', 'TFRC', 'EPOR', 'CD44'
    #     ],

    #     # 4) Basophilic erythroblast (early terminal)
    #     # TFRC high, GYPA+, ITGA4 often higher earlier, SLC4A1 starting
    #     "Basophilic_Erythroblast": [
    #         'GYPA', 'TFRC', 'ITGA4', 'SLC4A1'
    #     ],

    #     # 5) Polychromatic erythroblast (mid terminal)
    #     # TFRC starts trending down, Band3 (SLC4A1) up
    #     "Polychromatic_Erythroblast": [
    #         'GYPA', 'SLC4A1', 'TFRC'
    #     ],

    #     # 6) Orthochromatic erythroblast (late terminal; pre-enucleation)
    #     # TFRC low, Band3 high, ITGA4 down
    #     "Orthochromatic_Erythroblast": [
    #         'GYPA', 'SLC4A1', 'ITGA4'
    #     ],

    #     # 7) Reticulocyte (enucleated but still immature)
    #     # still can have TFRC, strong GYPA/SLC4A1
    #     "Reticulocyte": [
    #         'GYPA', 'SLC4A1', 'TFRC'
    #     ],

    #     # 8) Mature RBC
    #     # strong GYPA/SLC4A1, TFRC absent/very low
    #     "RBC_Mature_Erythrocyte": [
    #         'GYPA', 'SLC4A1'
    #     ], }
    },	


    "Platelets": {
        "general": ["RGS18", "C2orf88", "TMEM40", "GP9",
                    "PF4", "PPBP", "DAB2", "SPARC",
                    "RUFY1", "F13A1"],
        "subtypes": {}
    }
}