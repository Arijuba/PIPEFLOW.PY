import numpy as np

def res_oil_simp(av, h, k_o, oil, p_in, q_sc, r_in, r_out, rho_sc, semi, S, T_R):
    """
    Calcula a queda de pressão para fluxo monofásico de óleo em um reservatório perfurado
    por um poço vertical sob condições de estado estacionário ou semi-estacionário.
    
    Parâmetros:
    av: Seleção da pressão do reservatório (0 = pressão na fronteira, 1 = pressão média)
    h: Altura do reservatório, m
    k_o: Permeabilidade efetiva ao óleo, m^2
    oil: Modelo de óleo (1 = Standing, 2 = Glaso, 3 = óleo volátil)
    p_in: Pressão de entrada, Pa
    q_sc: Vazões padronizadas [q_g_sc, q_o_sc, q_w_sc], m³/s
    r_in: Raio inicial, m (r_e ou r_w)
    r_out: Raio final, m (r_w se r_in = r_e ou vice-versa)
    rho_sc: Densidades padronizadas [rho_g_sc, rho_o_sc, rho_w_sc], kg/m³
    semi: Seleção das condições de escoamento (0 = estacionário, 1 = semi-estacionário)
    S: Fator de dano (skin)
    T_R: Temperatura do reservatório, °C
    
    Retorna:
    p_out: Pressão de saída, Pa
    r: Distâncias radiais, m
    p: Pressão ao longo de r, Pa
    J: Índice de produtividade, m³/(s Pa)
    """
    rho_g_sc, rho_o_sc, _ = rho_sc
    q_g_sc, q_o_sc, q_w_sc = q_sc
    
    if q_w_sc != 0:
        raise ValueError("Taxa de água não nula. Uso de res_oil_simp inválido.")
    
    R_go = q_g_sc / q_o_sc if q_o_sc != 0 else 0  
    R_sb = R_go  
    
    # Cálculo das propriedades do óleo
    mu_o = oil_viscosity(p_in, R_sb, rho_g_sc, rho_o_sc, T_R)  # Viscosidade do óleo, Pa.s
    q_o = q_o_sc  # Taxa de óleo local, m³/s
    
    r_e = max(r_in, r_out)
    r_w = min(r_in, r_out)
    
    if semi == 0 and av == 0:
        f_R = 0
    elif semi == 0 and av == 1:
        f_R = 0.5
    elif semi == 1 and av == 0:
        f_R = 0.5
    else:
        f_R = 0.75
    
    Delta_p_skin = -mu_o * q_o * S / (2 * np.pi * k_o * h)
    help01 = -mu_o * q_o * (np.log(r_e / r_w) - f_R + S) / (2 * np.pi * k_o * h)
    
    if r_out > r_in:
        p_out = p_in + help01
        p_wf = p_in
    else:
        p_out = p_in - help01
        p_wf = p_out
    
    if oil == 1:
        p_b = pres_bub_Standing(R_go, rho_g_sc, rho_o_sc, T_R)
    elif oil == 2:
        p_b = pres_bub_Glaso(R_go, rho_g_sc, rho_o_sc, T_R)
    else:
        p_b = pres_bub_volatile_oil(T_R)
    
    if p_wf < p_b:
        raise ValueError("Pressão abaixo da pressão de bolha. Uso de res_oil_simp inválido.")
    
    n_p = 100
    r_inc = (r_e - r_w) / (n_p - 1)
    r = np.linspace(r_w, r_e, n_p)
    p = p_wf + Delta_p_skin + (p_out - (p_wf - Delta_p_skin)) * np.log(r / r_w) / np.log(r_e / r_w)
    
    if av == 1:
        J = -q_o_sc / (p_out - p_wf)
    else:
        J = -q_o_sc / (p_out - p_wf)
    
    return p_out, r, p, J