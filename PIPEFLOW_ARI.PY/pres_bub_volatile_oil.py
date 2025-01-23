import numpy as np

def pres_bub_volatile_oil(T, vol_oil):
    """
    Determina a pressão de bolha de um óleo volátil tabulado para uma dada temperatura.
    
    Parâmetros:
    T: Temperatura, °C
    vol_oil: Tabela de dados do óleo volátil (numpy array de dimensões apropriadas)
    
    Retorna:
    p_b: Pressão de bolha, Pa
    """
    # Encontrar a temperatura mais próxima na tabela
    i = 0
    T_lo = vol_oil[0, 0, 0]  # Menor valor de temperatura na tabela, °C
    T_tab = T_lo
    
    while T_tab < T and i < len(vol_oil) - 1:
        i += 1
        T_tab = vol_oil[i, 0, 0]
    
    # Percorrer os valores de pressão até que não haja mais aumento em R_s
    j = 0
    R_s_lo = vol_oil[i, 0, 5]  # Menor valor de R_s na tabela, m³/m³
    R_s_tab = R_s_lo
    R_s_old = R_s_lo - np.finfo(float).eps  # Pequeno decremento para iniciar o loop
    
    while R_s_tab > R_s_old and j < vol_oil.shape[1] - 1:
        j += 1
        R_s_old = R_s_tab
        R_s_tab = vol_oil[i, j, 5]
    
    p_b = vol_oil[i, j, 1]  # Pressão de bolha correspondente, Pa
    
    return p_b