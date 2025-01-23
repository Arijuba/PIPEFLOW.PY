import numpy as np

def pres_bub_Standing(R_sb, rho_g_sc, rho_o_sc, T):
    """
    Calcula a pressão de bolha usando a correlação de Standing.
    
    Parâmetros:
    R_sb: Razão gás-óleo na pressão de bolha, m³/m³
    rho_g_sc: Densidade do gás em condições padrão, kg/m³
    rho_o_sc: Densidade do óleo em condições padrão, kg/m³
    T: Temperatura, °C
    
    Retorna:
    p_b: Pressão de bolha, Pa
    """
    if rho_g_sc == 0:
        return 1.0e5  # Pressão atmosférica
    
    help01 = (10 ** (0.00164 * T)) / (10 ** (1768 / rho_o_sc))
    p_b = 125e3 * ((716 * R_sb / rho_g_sc) ** 0.83 * help01 - 1.4)
    
    return max(p_b, 1.0e5)  # Garantir que a pressão não seja menor que a atmosférica
