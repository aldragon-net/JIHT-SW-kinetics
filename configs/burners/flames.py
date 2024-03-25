from configs.burners.tprofile import read_temperature_profile
from configs.burners.impinging_jet_data import ImpingingJetData


ethylene_flame = ImpingingJetData(
    label="ethylene_base",
    fuel='C2H4:1',
    fuel_flow_lph=84,
    oxydizer='O2:0.21, N2:0.78, AR:0.01',
    oxydizer_flow_lph=573,
    height=0.023,
    T_room=293,
    T_burner=400,
    T_body=600,
    T_profile=read_temperature_profile('mckenna-constructed-C2H4.dat')
)

acetylene_flame = ImpingingJetData(
    label="acetylene_base",
    fuel='C2H2:1',
    fuel_flow_lph=84,
    oxydizer='O2:0.21, N2:0.78, AR:0.01',
    oxydizer_flow_lph=545,
    height=0.023,
    T_room=293,
    T_burner=400,
    T_body=600,
    T_profile=None
)

dme_fractions = [0, 15, 30, 45, 60, 75, 90, 100]
dme_flames = []
for dme_fraction in dme_fractions:
    dme_flames.append(
        ImpingingJetData(
            label=f'C2H4_{100-dme_fraction}_DME_{dme_fraction}',
            fuel=f'C2H4:{100-dme_fraction} CH3OCH3:{dme_fraction}',
            fuel_flow_lph=84,
            oxydizer='O2:0.21, N2:0.78, AR:0.01',
            oxydizer_flow_lph=545,
            height=0.023,
            T_room=293,
            T_burner=400,
            T_body=600,
            T_profile=None
        )
    )


ethylene_est = ImpingingJetData(
    label="ethylene_base",
    fuel='C2H4:1',
    fuel_flow_lph=84,
    oxydizer='O2:0.21, N2:0.78, AR:0.01',
    oxydizer_flow_lph=573,
    height=0.023,
    T_room=293,
    T_burner=400,
    T_body=600,
    T_profile=None
)

ethylene_he25 = ImpingingJetData(
    label="ethylene_he",
    fuel='C2H4:0.75 AR:0.25',
    fuel_flow_lph=84,
    oxydizer='O2:0.21, N2:0.78, AR:0.01',
    oxydizer_flow_lph=573,
    height=0.023,
    T_room=293,
    T_burner=400,
    T_body=600,
    T_profile=None
)

ethylene_h2 = ImpingingJetData(
    label="ethylene_h2",
    fuel='C2H4:0.75 H2:0.25',
    fuel_flow_lph=84,
    oxydizer='O2:0.21, N2:0.78, AR:0.01',
    oxydizer_flow_lph=573,
    height=0.023,
    T_room=293,
    T_burner=400,
    T_body=600,
    T_profile=None
)

ethylene_c3h8 = ImpingingJetData(
    label="ethylene_c3h8",
    fuel='C2H4:0.75 C3H8:0.25',
    fuel_flow_lph=84,
    oxydizer='O2:0.21, N2:0.78, AR:0.01',
    oxydizer_flow_lph=573,
    height=0.023,
    T_room=293,
    T_burner=400,
    T_body=600,
    T_profile=None
)
