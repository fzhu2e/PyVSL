def convert_unit_P(P_in, unit_in='kg/m2/s', unit_out='mm/month'):
    ''' Convert the unit for precipitation
    '''
    if unit_in == 'kg/m2/s' and unit_out == 'mm/month':
        P_out = P_in * 3600*24*30
    elif unit_in == 'mm/month' and unit_out == 'kg/m2/s':
        P_out = P_in / (3600*24*30)

    return P_out

def convert_unit_T(T_in, unit_in='K', unit_out='C'):
    ''' Convert the unit for temperature
    '''
    if unit_in == 'K' and unit_out == 'C':
        T_out = T_in - 273.15
    elif unit_in == 'C' and unit_out == 'K':
        T_out = T_in + 273.15

    return T_out