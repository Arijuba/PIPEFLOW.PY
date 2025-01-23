import numpy as np
from scipy.integrate import solve_ivp

def res(beta, fluid, h, k, oil, p_in, q_sc, r_in, r_out, rel, rho_sc, T_R):
    """
    Calcula a queda de pressão para o fluxo em um reservatório drenado por um poço vertical
    produzindo de uma área cilíndrica de drenagem sob condições de estado estacionário.
    
    Parâmetros:
    beta: Coeficiente de Forchheimer, m^-1
    fluid: Tipo de fluido (1 = óleo, 2 = gás, >=3 = multifásico)
    h: Altura do reservatório, m
    k: Permeabilidade absoluta, m^2
    oil: Modelo de óleo (1 = Standing, 2 = Glaso, 3 = óleo volátil)
    p_in: Pressão de entrada, Pa
    q_sc: Vazões padronizadas [q_g_sc, q_o_sc, q_w_sc], m³/s
    r_in: Raio inicial, m (r_e ou r_w)
    r_out: Raio final, m (r_w se r_in = r_e ou vice-versa)
    rel: Parâmetros de permeabilidade relativa
    rho_sc: Densidades padronizadas [rho_g_sc, rho_o_sc, rho_w_sc], kg/m³
    T_R: Temperatura do reservatório, °C
    
    Retorna:
    p_out: Pressão de saída, Pa
    r: Coordenadas radiais, m
    p: Pressão ao longo de r, Pa
    """
    intervalo = [r_in, r_out]  # Intervalo de integração, m
    boundcon = [p_in]  # Condição de contorno, Pa
    
    if fluid == 1:  # Óleo
        k_o = k  # Permeabilidade efetiva para óleo, m^2
        resultado = solve_ivp(res_oil_dpdr, intervalo, boundcon, args=(h, k_o, oil, q_sc, rho_sc, T_R))
    elif fluid == 2:  # Gás
        k_g = k  # Permeabilidade efetiva para gás, m^2
        resultado = solve_ivp(res_gas_dpdr, intervalo, boundcon, args=(beta, h, k_g, q_sc, rho_sc, T_R))
    else:  # Multifásico
        resultado = solve_ivp(res_dpdr, intervalo, boundcon, args=(beta, h, k, oil, q_sc, rel, rho_sc, T_R))
    
    r = resultado.t
    p = resultado.y[0]
    p_out = p[-1]
    
    return p_out, r, p
