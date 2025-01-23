import numpy as np

def res_gas_simp(av, h, k_g, p_in, q_sc, r_in, r_out, rho_sc, semi, S, T_R):
    """
    Calcula a queda de pressão para fluxo monofásico de gás em um reservatório perfurado
    por um poço vertical sob condições de estado estacionário ou semi-estacionário.
    
    Parâmetros:
    av: Seleção da pressão do reservatório (0 = pressão na fronteira, 1 = pressão média)
    h: Altura do reservatório, m
    k_g: Permeabilidade efetiva ao gás, m^2
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
    J_pseu: Índice de produtividade pseudo, m³/(Pa s)
    """
    rho_g_sc, _, _ = rho_sc
    q_g_sc, q_o_sc, q_w_sc = q_sc
    
    if q_o_sc != 0 or q_w_sc != 0:
        raise ValueError("Taxa de óleo ou água não nula. Uso de res_gas_simp inválido para fluxo monofásico de gás.")
    
    r_e = max(r_in, r_out)
    r_w = min(r_in, r_out)
    
    p_pc = pres_pseu_crit_Sutton(rho_g_sc)
    T_pc = temp_pseu_crit_Sutton(rho_g_sc)
    T_R_abs = T_R + 273.15
    
    p_R = p_in
    p_R_av = p_in
    
    f_damp = 0.5
    tol_abs = 1.e2
    tol_rel = 1.e-3
    max_iter = 1000
    iter_count = 0
    repeat = True
    
    while repeat and iter_count < max_iter:
        iter_count += 1
        p_R_old, p_R_av_old = p_R, p_R_av
        
        p_pr_av = p_R_av / p_pc
        T_pr_av = T_R_abs / T_pc
        Z_av = Z_factor_DAK(p_pr_av, T_pr_av)
        
        p_sc = 1.e5
        T_sc_abs = 15 + 273.15
        Z_sc = 1
        B_g_av = (p_sc * T_R_abs * Z_av) / (p_R_av * T_sc_abs * Z_sc)
        mu_g_av = gas_viscosity(p_R_av, rho_g_sc, T_R)
        
        f_R = 0 if semi == 0 and av == 0 else 0.5 if semi == 0 and av == 1 else 0.5 if semi == 1 and av == 0 else 0.75
        help01 = mu_g_av * B_g_av * q_g_sc * (np.log(r_e / r_w) - f_R + S) / (np.pi * k_g * h)
        
        if r_out > r_in:
            if av == 1:
                p_wf = p_in
                p_R_av = np.sqrt(abs(p_wf**2 - help01 * p_R_av))
                p_out = p_R_av
                p_R_new = p_R_av - mu_g_av * B_g_av * q_g_sc / (4 * np.pi * k_g * h)
                p_R = f_damp * p_R_old + (1 - f_damp) * p_R_new
            else:
                p_wf = p_in
                p_R = np.sqrt(abs(p_wf**2 - help01 * p_R_av))
                p_out = p_R
                p_R_av_new = p_R + mu_g_av * B_g_av * q_g_sc / (4 * np.pi * k_g * h)
                p_R_av = f_damp * p_R_av_old + (1 - f_damp) * p_R_av_new
        else:
            if av == 1:
                p_R_av = p_in
                p_wf = np.sqrt(abs(p_R_av**2 + help01 * p_R_av))
                p_out = p_wf
                p_R_new = p_wf + (p_R_av - p_wf) * 3 / 2
                p_R = f_damp * p_R_old + (1 - f_damp) * p_R_new
            else:
                p_R = p_in
                p_wf = np.sqrt(abs(p_R**2 + help01 * p_R_av))
                p_out = p_wf
                p_R_av_new = p_wf + (p_R - p_wf) * 2 / 3
                p_R_av = f_damp * p_R_av_old + (1 - f_damp) * p_R_av_new
        
        diff = p_R - p_R_old if av == 1 else p_R_av - p_R_av_old
        rel_diff = diff / p_R_old if av == 1 else diff / p_R_av_old
        repeat = abs(diff) > tol_abs or abs(rel_diff) > tol_rel
    
    if iter_count >= max_iter:
        raise ValueError("Número máximo de iterações excedido.")
    
    n_p = 100
    r = np.linspace(r_w, r_e, n_p)
    p = np.sqrt(p_wf**2 + (p_R**2 - p_wf**2) * np.log(r / r_w) / np.log(r_e / r_w))
    J_pseu = mu_g_av * B_g_av * p_R_av * (np.log(r_w / r_e) - f_R + S) / (np.pi * k_g * h * p_R)
    
    return p_out, r, p, J_pseu
