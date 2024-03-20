from configs.burners.tprofile import read_temperature_profile
from configs.burners.impinging_jet_data import ImpingingJetData

ethylene_flame = ImpingingJetData(
    label="ethylene_base",
    mixture='C2H4:1 O2:1.43 N2:5.32 AR:0.068',
    inlet_velocity=0.0645,
    height=0.023,
    T_room=293,
    T_burner=400,
    T_body=600,
    T_profile=read_temperature_profile('mckenna-constructed-C2H4.dat')
)

acetylene_flame = ImpingingJetData(
    label="acetylene_base",
    mixture='C2H2:1 O2:1.363 N2:5.06 AR:0.064',
    inlet_velocity=0.0617,
    height=0.023,
    T_room=293,
    T_burner=400,
    T_body=600,
    T_profile=None
)
