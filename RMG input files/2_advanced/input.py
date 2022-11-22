# Adding libraries, pressure dependence, and species constraints
# for your projects, you need to reconsider toleranceMoveToCore and the termination criteria

# Data sources
database(
    thermoLibraries=['primaryThermoLibrary',
    'BurkeH2O2',
    'thermo_DFT_CCSDTF12_BAC',
    'DFT_QCI_thermo',
    'Spiekermann_refining_elementary_reactions',
    'CurranPentane',
    'C3',
    'CBS_QB3_1dHR',
    'FFCM1(-)',
    'JetSurF2.0',
    's3_5_7_ane',
    'naphthalene_H',
    'USC-Mech-ii',
    'heavy_oil_ccsdtf12_1dHR',
    'bio_oil',
    'vinylCPD_H',
    'Klippenstein_Glarborg2016',
    'Fulvene_H',
    'Chernov',
    'C10H11',
    'CH'],
    reactionLibraries=['primaryH2O2',
    'C2H2_init',
    'C2H4+O_Klipp2017',
    'CurranPentane',
    'FFCM1(-)',
    'NOx2018',
    'JetSurF2.0',
    'Klippenstein_Glarborg2016',
    'C10H11',
    'C12H11_pdep',
    'Lai_Hexylbenzene',
    '1989_Stewart_2CH3_to_C2H5_H',
    '2001_Tokmakov_H_Toluene_to_CH3_Benzene',
    '2003_Miller_Propargyl_Recomb_High_P',
    '2005_Senosiain_OH_C2H2',
    'kislovB',
    'c-C5H5_CH3_Sharma',
    'fascella',
    '2006_Joshi_OH_CO',
    '2009_Sharma_C5H5_CH3_highP',
    '2015_Buras_C2H3_C4H6_highP',
    'C3',
    'Methylformate',
    'C6H5_C4H4_Mebel',
    'vinylCPD_H',
    'Mebel_Naphthyl',
    'Mebel_C6H5_C2H2',
    'Fulvene_H'],
    transportLibraries=['OneDMinN2', 'PrimaryTransportLibrary', 'NOx2018', 'GRI-Mech'],
    seedMechanisms=[],
    kineticsDepositories='default',
    kineticsFamilies='default',
    kineticsEstimator='rate rules',
)

# List of species
species(
    label='ethane',
    structure=SMILES("CC"),
)
species(
    label='O2',
    structure=SMILES('[O][O]'),
)
species(
    label='N2',
    reactive=False,
    structure=SMILES('N#N'),
)

# Reaction systems
simpleReactor(
    temperature=(1350,'K'),
    pressure=(1.0,'bar'),
    initialMoleFractions={
        "ethane": 1.0,
        "O2": 3.5,
        "N2": 12.845,
    },
    terminationTime=(1,'ms'),
    terminationRateRatio=0.1,
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

pressureDependence(
    method='modified strong collision',
    maximumGrainSize=(2.0, 'kJ/mol'),
    minimumNumberOfGrains=250,
    temperatures=(300, 2500, 'K', 10),
    pressures=(0.1, 100, 'bar', 10),
    interpolation=('Chebyshev', 6, 4),
    maximumAtoms=16,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.20,
    toleranceInterruptSimulation=0.20,
    maximumEdgeSpecies=100000,
    filterReactions=True,
)

options(
    units='si',
    generateOutputHTML=True,
    generatePlots=False,
    saveEdgeSpecies=True,
    saveSimulationProfiles=True,
)

generatedSpeciesConstraints(
    allowed=['input species', 'seed mechanisms', 'reaction libraries'],
    maximumCarbonAtoms=4,
    maximumOxygenAtoms=2,
    maximumNitrogenAtoms=0,
    maximumSiliconAtoms=0,
    maximumSulfurAtoms=0,
    maximumHeavyAtoms=5,
    maximumRadicalElectrons=2,
    maximumSingletCarbenes=1,
    maximumCarbeneRadicals=0,
    allowSingletO2=True,
)

