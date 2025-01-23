import numpy as np

def from_kg_per_m3_to_gas_grav(rho_g_sc):
    return rho_g_sc / 1.225  # Gravidade específica do gás em relação ao ar

def from_kg_per_m3_to_deg_API(rho_o_sc):
    return (141.5 / (rho_o_sc / 999.1)) - 131.5  # Grau API

def from_m3_per_m3_to_ft3_per_bbl(R_sb):
    return R_sb * 5.615  # Conversão de m³/m³ para ft³/bbl

def from_deg_C_to_deg_F(T):
    return (T * 9/5) + 32  # Conversão de graus Celsius para Fahrenheit

def from_psi_to_Pa(p_psi):
    return p_psi * 6894.76  # Conversão de psi para Pa

def pres_bub_Glaso(R_sb, rho_g_sc, rho_o_sc, T):
    """
    Calcula a pressão de bolha usando a correlação de Glaso.
    
    Parâmetros:
    R_sb: Razão gás-óleo na pressão de bolha, m³/m³
    rho_g_sc: Densidade do gás em condições padrão, kg/m³
    rho_o_sc: Densidade do óleo em condições padrão, kg/m³
    T: Temperatura, °C
    
    Retorna:
    p_b: Pressão de bolha, Pa
    """
    gamma_g = from_kg_per_m3_to_gas_grav(rho_g_sc)
    gamma_API = from_kg_per_m3_to_deg_API(rho_o_sc)
    R_sb_FU = from_m3_per_m3_to_ft3_per_bbl(R_sb)
    T_FU = from_deg_C_to_deg_F(T)
    
    a, b, c = 0.816, 0.172, 0.989
    p_b_star = (R_sb_FU / gamma_g) ** a * T_FU ** b / gamma_API ** c
    p_b_FU = 10 ** (1.7669 + 1.7447 * np.log10(p_b_star) - 0.30218 * (np.log10(p_b_star)) ** 2)
    
    return from_psi_to_Pa(p_b_FU)
